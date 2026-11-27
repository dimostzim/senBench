#!/usr/bin/env python3
"""Build the sample table and the label inputs from the downloaded records of GEO and ArrayExpress, in four parts:

1. The sample and series fields of every GEO series are copied from its SOFT file into tables, unchanged.
2. The sample sheets (SDRF) of the ArrayExpress experiments are copied into one long table, unchanged.
3. The sample table. Every sample gets organism, assay, tissue, cell type, cell line, sex, age, genotype, treatment,
   time and passage, taken from the repository fields by name. Samples of datasets with relevance `direct` also get a
   senescence state, an inducer and a kind of control, derived from the wording of their own fields. A dataset without
   a sample sheet gets one row that describes the dataset.
4. The label inputs: per labelled dataset, the series and sample fields that a label reading reads, unchanged, from
   GEO, or for an ArrayExpress experiment from its sample sheet and its investigation description (IDF).
"""

import argparse
import csv
import gzip
import itertools
import json
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# The folders of download.sh that hold the GEO series and the ArrayExpress experiments.
GEO_FOLDER = "geo"
ARRAYEXPRESS_FOLDER = "biostudies_arrayexpress"
# Joins the values of a field that occurs several times.
SEPARATOR = "; "
GEO_SAMPLES = "geo_samples.tsv.gz"
GEO_CHARACTERISTICS = "geo_sample_characteristics.tsv.gz"
GEO_SERIES = "geo_series.tsv"
SHEETS_NAME = "sample_sheets_downloads.tsv.gz"
SAMPLES_NAME = "samples.tsv.gz"
COVERAGE_NAME = "sample_table_coverage.tsv"
STATES_NAME = "sample_table_senescence_states.tsv"
LABEL_INPUTS_NAME = "label_inputs.jsonl.gz"
GEO_UNITS_NAME = "units_geo_records.tsv"
SHEET_UNITS_NAME = "units_sample_sheets.tsv"
# Assay values of the sample table that are bulk RNA-seq: GEO's library strategy, and the library strategy of an
# ArrayExpress sample sheet joined to its technology type.
BULK_RNA_SEQ_ASSAYS = ("RNA-Seq", "RNA-Seq; sequencing assay")

csv.field_size_limit(sys.maxsize)


def _one_line(text):
    return " ".join(str(text).split())


# 1. The fields of the GEO series.

CHANNELS = ("ch1", "ch2")
# Table column -> SOFT field. A field that occurs several times in a record is joined with SEPARATOR.
SOFT_SAMPLE_FIELDS = {
    "title": "Sample_title",
    "type": "Sample_type",
    "channel_count": "Sample_channel_count",
    "source_name_ch1": "Sample_source_name_ch1",
    "organism_ch1": "Sample_organism_ch1",
    "taxid_ch1": "Sample_taxid_ch1",
    "characteristics_ch1": "Sample_characteristics_ch1",
    "molecule_ch1": "Sample_molecule_ch1",
    "treatment_protocol_ch1": "Sample_treatment_protocol_ch1",
    "growth_protocol_ch1": "Sample_growth_protocol_ch1",
    "source_name_ch2": "Sample_source_name_ch2",
    "organism_ch2": "Sample_organism_ch2",
    "characteristics_ch2": "Sample_characteristics_ch2",
    "description": "Sample_description",
    "platform": "Sample_platform_id",
    "instrument_model": "Sample_instrument_model",
    "library_strategy": "Sample_library_strategy",
    "library_source": "Sample_library_source",
    "library_selection": "Sample_library_selection",
    "relations": "Sample_relation",
    "supplementary_files": "Sample_supplementary_file",
    "series_of_sample": "Sample_series_id",
    "submission_date": "Sample_submission_date",
}
SOFT_SERIES_FIELDS = {
    "title": "Series_title",
    "types": "Series_type",
    "summary": "Series_summary",
    "overall_design": "Series_overall_design",
    "pubmed_ids": "Series_pubmed_id",
    "relations": "Series_relation",
    "platforms": "Series_platform_id",
    "submission_date": "Series_submission_date",
    "last_update_date": "Series_last_update_date",
}
GEO_SAMPLE_COLUMNS = ("series", "sample", *SOFT_SAMPLE_FIELDS)
GEO_SERIES_COLUMNS = ("series", *SOFT_SERIES_FIELDS, "samples")
GEO_CHARACTERISTICS_COLUMNS = ("series", "sample", "channel", "key", "value")


def _records(path):
    """Yield (kind, accession, fields) for every record of a SOFT file; fields maps a field name to its values."""
    kind = accession = None
    fields = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[0] == "^":
                if kind:
                    yield kind, accession, fields
                kind, _, accession = line[1:].partition("=")
                kind, accession, fields = kind.strip(), accession.strip(), {}
            # Lines without "!" are rows of the data tables, which are not read.
            elif line[0] == "!" and " = " in line:
                name, _, value = line[1:].partition(" = ")
                # A supplementary file field carries a number, as in Sample_supplementary_file_1.
                if name.startswith("Sample_supplementary_file"):
                    name = "Sample_supplementary_file"
                fields.setdefault(name, []).append(_one_line(value))
    if kind:
        yield kind, accession, fields


def _row(accession_columns, columns, fields):
    return {**accession_columns, **{column: SEPARATOR.join(fields.get(name, [])) for column, name in columns.items()}}


def _harvest(series_dir):
    """Return (series, status, note, series row, sample rows, characteristics rows) for one series folder."""
    series = series_dir.name
    soft = series_dir / "soft" / f"{series}_family.soft.gz"
    if not soft.is_file():
        return series, "failed", "no SOFT file", None, [], []
    series_row = None
    samples = []
    characteristics = []
    try:
        for kind, accession, fields in _records(soft):
            if kind == "SERIES":
                series_row = _row({"series": series}, SOFT_SERIES_FIELDS, fields)
            elif kind == "SAMPLE":
                samples.append(_row({"series": series, "sample": accession}, SOFT_SAMPLE_FIELDS, fields))
                for channel in CHANNELS:
                    for entry in fields.get(f"Sample_characteristics_{channel}", []):
                        key, separator, value = entry.partition(": ")
                        characteristics.append((series, accession, channel, key if separator else "",
                                                value if separator else entry))
    except (OSError, EOFError) as problem:
        return series, "failed", f"SOFT file cannot be read: {problem}", None, [], []
    if series_row is None:
        return series, "failed", "no series record in the SOFT file", None, [], []
    series_row["samples"] = len(samples)
    return series, "ok", f"{len(samples)} samples", series_row, samples, characteristics


def harvest_geo_records(geo_dir, out_dir, workers):
    """Write the tables of the sample fields, the characteristics and the series fields of every GEO series."""
    series_dirs = sorted(path for path in geo_dir.iterdir() if path.is_dir())
    failed = []
    with gzip.open(out_dir / GEO_SAMPLES, "wt", newline="") as samples_file, \
            gzip.open(out_dir / GEO_CHARACTERISTICS, "wt", newline="") as characteristics_file, \
            (out_dir / GEO_SERIES).open("w", newline="") as series_file, \
            (out_dir / GEO_UNITS_NAME).open("w", newline="") as units_file, \
            ProcessPoolExecutor(max_workers=workers) as pool:
        samples_writer = csv.DictWriter(samples_file, GEO_SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n")
        series_writer = csv.DictWriter(series_file, GEO_SERIES_COLUMNS, delimiter="\t", lineterminator="\n")
        characteristics_writer = csv.writer(characteristics_file, delimiter="\t", lineterminator="\n")
        units_writer = csv.writer(units_file, delimiter="\t", lineterminator="\n")
        samples_writer.writeheader()
        series_writer.writeheader()
        characteristics_writer.writerow(GEO_CHARACTERISTICS_COLUMNS)
        results = pool.map(_harvest, series_dirs, chunksize=8)
        for done, (series, status, note, series_row, samples, characteristics) in enumerate(results, start=1):
            units_writer.writerow((series, status, note))
            if status == "ok":
                series_writer.writerow(series_row)
                samples_writer.writerows(samples)
                characteristics_writer.writerows(characteristics)
            else:
                failed.append(series)
            if done % 500 == 0 or done == len(series_dirs):
                print(f"PROGRESS {done}/{len(series_dirs)} {series} {status}: {note}", flush=True)
    if failed:
        print(f"failed series ({len(failed)}): {' '.join(failed)}", flush=True)
        sys.exit(1)


# 2. The sample sheets of the ArrayExpress experiments.

SHEET_COLUMNS = ("dataset", "data_folder", "sheet", "row", "column", "value")
SHEET_PATTERN = "*.sdrf.txt"
PROBE_BYTES = 4096


class _NotATable(Exception):
    """A file that carries the name of a sample sheet and is not a text table."""


def _sheet_rows(path):
    """Yield (row number, column, value) for every filled cell of a tab-separated sheet."""
    with path.open("rb") as handle:
        if b"\x00" in handle.read(PROBE_BYTES):
            raise _NotATable
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [_one_line(name) for name in next(reader, [])]
        for number, row in enumerate(reader, start=1):
            for column, value in zip(header, row):
                if value.strip():
                    yield number, column, _one_line(value)


def harvest_sample_sheets(arrayexpress_dir, out_dir):
    """Write every filled cell of the sample sheets of the ArrayExpress experiments into one long table."""
    failed = 0
    sheets = 0
    with gzip.open(out_dir / SHEETS_NAME, "wt", newline="") as out_file, \
            (out_dir / SHEET_UNITS_NAME).open("w", newline="") as units_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        units = csv.writer(units_file, delimiter="\t", lineterminator="\n")
        writer.writerow(SHEET_COLUMNS)
        for dataset in sorted(path for path in arrayexpress_dir.glob("*") if path.is_dir()):
            for sheet in sorted(dataset.rglob(SHEET_PATTERN)):
                name = str(sheet.relative_to(dataset))
                sheets += 1
                try:
                    rows = list(_sheet_rows(sheet))
                except _NotATable:
                    units.writerow((f"{dataset.name}/{name}", "skipped", "not a text table"))
                    continue
                except (OSError, csv.Error) as problem:
                    units.writerow((f"{dataset.name}/{name}", "failed", f"cannot be read: {problem!r}"))
                    failed += 1
                    continue
                writer.writerows((dataset.name, ARRAYEXPRESS_FOLDER, name, number, column, value)
                                 for number, column, value in rows)
                units.writerow((f"{dataset.name}/{name}", "ok", f"{len({number for number, _, _ in rows})} rows"))
                if sheets % 500 == 0:
                    print(f"PROGRESS {sheets} sheets read, the last: {dataset.name}/{name}", flush=True)
    print(f"{sheets} sheets read, {failed} failed", flush=True)
    if failed:
        sys.exit(1)


# 3. The sample table.

DIRECT = "direct"
COLUMNS = ("dataset", "repository", "relevance", "level", "sample", "title", "organism", "assay", "tissue", "cell_type",
           "cell_line", "sex", "age", "genotype", "treatment", "time", "passage", "senescence_state", "senescence_rule",
           "senescence_inducer", "control_kind", "senescence_evidence")
COVERAGE_COLUMNS = ("repository", "relevance", "level", "rows", "datasets", *COLUMNS[6:17])
STATE_COLUMNS = ("repository", "senescence_state", "senescence_rule", "samples", "datasets")
# Derived column -> names of the repository fields it is taken from, after _field_name().
FIELD_NAMES = {
    "organism": ("organism",),
    "tissue": ("tissue", "organ", "tissue type", "tissue source", "tissue of origin", "tissue origin", "source tissue",
               "body site", "anatomical site", "organism part", "organ hierarchy", "tissue/organ", "organ/tissue"),
    "cell_type": ("cell type", "cell types", "celltype", "cell population", "cell subtype", "cell lineage",
                  "cell subset"),
    "cell_line": ("cell line", "cell lines", "cellline", "cell line name", "parental cell line",
                  "cell line background"),
    "sex": ("sex", "gender", "donor sex", "sex of donor"),
    "age": ("age", "donor age", "age at collection", "age at sampling", "age (years)", "age (yrs)", "age (months)",
            "age (weeks)", "age (days)", "age in years", "age in months", "age in weeks"),
    "genotype": ("genotype", "genotype/variation", "genotype variation", "genetic modification",
                 "genetic background"),
    "treatment": ("treatment", "agent", "treated with", "drug", "compound", "stimulus", "stimulation",
                  "treatment group", "drug treatment", "exposure", "perturbation"),
    "time": ("time", "time point", "timepoint", "time post treatment", "days post treatment", "treatment duration",
             "treatment time", "duration", "day"),
    "passage": ("passage", "passages", "passage number", "population doubling", "population doublings",
                "population doubling level", "pd", "pdl"),
}
COLUMN_OF_FIELD = {name: column for column, names in FIELD_NAMES.items() for name in names}
SEX_VALUES = {"female": "female", "f": "female", "woman": "female", "women": "female", "male": "male", "m": "male",
              "man": "male", "men": "male"}
GEO_ARRAY_ASSAYS = {"RNA": "array, RNA", "genomic": "array, genomic DNA", "protein": "array, protein"}
# Fields of a sample sheet that describe the sample; its other columns name files and protocols.
_SHEET_SAMPLE_FIELD = re.compile(r"(?:characteristics|factor ?value)\s*\[", re.I)
# A longer field value is a description of the experiment, not of the sample.
LABEL_VALUE_MAX = 80
RULE_STATED = "stated in the fields that differ between the samples of the dataset"
RULE_INDUCER = "an inducer is named, the state is not stated"
RULE_CONTROL = "control wording in a dataset with samples stated to be senescent"
RULE_SHARED = "stated only in wording shared by all samples of the dataset (weak)"
RULE_MIXED = "wording for both states, or for a bypass, reversal or senolytic"
RULE_TWO_CHANNEL = "two-channel array: the sample compares two conditions"
# Wording rules for the senescence state. A negated form is turned into one word first, so that it does not match
# the plain form: "non-senescent" becomes "nonsenescent".
_NEGATED = re.compile(r"\b(?:non|not|pre|un)[- ]?(senescen\w*)|\bno (senescen\w*)", re.I)
_SPECIAL = re.compile(r"senolytic|bypass|escap|revers|post[- ]?senescen|senescence[- ]resist|rejuvenat|reprogram", re.I)
_SENESCENT = re.compile(r"\bsenescen|\bOIS\b|\bTIS\b|\bSIPS\b", re.I)
# In wording shared by all samples of a dataset, "senescence" is often the name of the experiment.
_SENESCENT_SHARED = re.compile(r"\bsenescent\b", re.I)
_NOT_SENESCENT = re.compile(r"nonsenescen", re.I)
_YES = re.compile(r"^(yes|true|positive|pos|high|\+|1|senescen\w*|induced|rs|ois|tis)$", re.I)
_NO = re.compile(r"^(no|false|negative|neg|low|-|0|none|nonsenescen\w*|control|ctrl|normal|proliferat\w*|young|"
                 r"growing|cycling)$", re.I)
# Wording that states that the cells divide or rest; it labels a sample on its own.
STATED_CONTROLS = (
    ("proliferating", re.compile(r"proliferat|\bcycling\b|\bgrowing\b|exponential", re.I)),
    ("quiescent", re.compile(r"quiescen|serum[- ]starv|contact[- ]inhibit", re.I)),
)
# A control of a reagent (siRNA, antibody) says nothing about the state of the cells.
_REAGENT_CONTROL = re.compile(r"\b(?:si|sh|sg)[- ]?(?:RNA[- ]?)?(?:controls?|ctrl|con|scr\w*)\b|"
                              r"\b(?:controls?|ctrl)[- ]?(?:si|sh|sg)(?:RNA)?\b|non[- ]?targeting(?: control)?|"
                              r"scrambled?(?: control)?|\bIgG(?: control)?\b|isotype(?: control)?|"
                              r"\binput(?: control)?\b|negative control", re.I)
# Wording of a control; it labels a sample only in a dataset where another sample is stated to be senescent.
CONTROLS = STATED_CONTROLS + (
    ("empty vector", re.compile(r"empty[- ]vector|vector[- ]control|\bEV\b", re.I)),
    ("early passage or young", re.compile(r"early[- ]passage|low[- ]passage|\byoung\b|\bearly PD", re.I)),
    ("untreated or vehicle", re.compile(r"untreated|non[- ]?treated|vehicle|\bDMSO\b|\bmock\b|\bcontrols?\b|\bctrl\b|"
                                        r"\bsham\b|non[- ]?irradiated|\b0 ?Gy\b", re.I)),
)
INDUCERS = (
    ("replicative", re.compile(r"replicative|late[- ]passage|high[- ]passage|\blate PD", re.I)),
    ("irradiation", re.compile(r"irradiat|\bIR\b|x[- ]?ray|gamma|γ|\b\d+(?:\.\d+)? ?Gy\b|radiation|\bUV[ABC]?\b", re.I)),
    ("oncogene", re.compile(r"\bOIS\b|oncogene[- ]induced|RAS[- ]?V12|RAS[- ]?G12V|\bG12V\b|ER:RAS|BRAF[- ]?V600E|"
                            r"RAS[- ]induced|\bH-?RAS", re.I)),
    ("DNA-damaging drug", re.compile(r"doxorubicin|\bdoxo\b|etoposide|bleomycin|cisplatin|camptothecin|mitomycin|"
                                     r"temozolomide|aphidicolin|hydroxyurea|busulfan", re.I)),
    ("CDK4/6 inhibitor", re.compile(r"palbociclib|abemaciclib|ribociclib|PD-?0332991", re.I)),
    ("oxidative stress", re.compile(r"\bH2O2\b|hydrogen peroxide|oxidative stress|\btBHP\b|paraquat", re.I)),
    ("therapy", re.compile(r"\bTIS\b|therapy[- ]induced|chemotherap", re.I)),
)


def _field_name(name):
    """Return the name of a repository field in one form: lower case, without the wrapper of a sample sheet."""
    name = name.strip()
    wrapped = re.fullmatch(r"(?:characteristics|factor ?value|comment|parameter ?value)\s*\[(.+)\]", name, re.I)
    name = wrapped.group(1) if wrapped else name.rsplit(".", 1)[-1]
    name = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", name).replace("_", " ")
    return " ".join(name.lower().split())


def _sex(values):
    known = {SEX_VALUES.get(value.strip().lower(), "") for value in values} - {""}
    if len(known) == 1:
        return known.pop()
    return "mixed" if known else ""


def _derived(fields):
    """Return the derived columns of a sample from its repository fields, a list of (field name, value)."""
    found = defaultdict(list)
    for name, value in fields:
        name = _field_name(name)
        column = COLUMN_OF_FIELD.get(name)
        if column and value not in found[column]:
            # An age field other than "age" says its unit or its time in its name.
            found[column].append(value if column != "age" or name == "age" else f"{name}: {value}")
    row = {column: SEPARATOR.join(values) for column, values in found.items()}
    row["sex"] = _sex(found.get("sex", []))
    return row


def _join_negations(text):
    return _NEGATED.sub(lambda match: "non" + (match.group(1) or match.group(2)), text)


def _snippet(text, match):
    return text[max(match.start() - 30, 0):match.end() + 30].strip()


def _first(patterns, text):
    """Return (name, wording) of the first named pattern found in the text, or ("", "")."""
    for name, pattern in patterns:
        match = pattern.search(text)
        if match:
            return name, _snippet(text, match)
    return "", ""


def _stated_state(text, fields, senescent=_SENESCENT):
    """Return (state, kind of control, evidence) from wording that states the senescence state, or None."""
    for name, value in fields:
        # A field named after senescence states it with yes or no, or with the state itself.
        if "senescen" in name.lower() and "senescen" not in value.lower():
            if _NO.match(value.strip()):
                return "not senescent", "non-senescent as stated", f"{name}: {value}"
            if _YES.match(value.strip()):
                return "senescent", "", f"{name}: {value}"
    special = _SPECIAL.search(text)
    if special:
        return "unclear", "", _snippet(text, special)
    negative = _NOT_SENESCENT.search(text)
    positive = senescent.search(text)
    kind, wording = _first(STATED_CONTROLS, text)
    # Wording for both states in one sample, as in "senescent versus proliferating", does not label it.
    if positive and (negative or kind):
        return "unclear", "", _snippet(text, positive)
    if negative:
        return "not senescent", "non-senescent as stated", _snippet(text, negative)
    if positive:
        return "senescent", "", _snippet(text, positive)
    if kind:
        return "not senescent", kind, wording
    return None


def _varying_parts(samples):
    """Split the text of every sample of a dataset into the part that differs between its samples and the rest.

    A word that every sample of the dataset carries, such as the name of the experiment, says nothing about one
    sample. Returns, per sample, (varying text, varying fields, constant text).
    """
    texts = [_join_negations(sample["title"]) for sample in samples]
    tokens = [re.findall(r"[A-Za-z0-9+:]+", text) for text in texts]
    fields = [[(name, _join_negations(value)) for name, value in sample["label_fields"]
               if len(value) <= LABEL_VALUE_MAX] for sample in samples]
    if len(samples) == 1:
        return [(f"{texts[0]} {' '.join(value for _, value in fields[0])}", fields[0], "")]
    lowered = [{token.lower() for token in sample_tokens} for sample_tokens in tokens]
    constant_tokens = set.intersection(*lowered)
    constant_fields = set.intersection(*(set(sample_fields) for sample_fields in fields))
    parts = []
    for sample_tokens, sample_fields in zip(tokens, fields):
        varying_fields = [field for field in sample_fields if field not in constant_fields]
        varying = " ".join(token for token in sample_tokens if token.lower() not in constant_tokens)
        constant = " ".join(token for token in sample_tokens if token.lower() in constant_tokens)
        parts.append((f"{varying} {' '.join(value for _, value in varying_fields)}", varying_fields,
                      f"{constant} {' '.join(value for _, value in constant_fields)}"))
    return parts


def _label_dataset(samples):
    """Set the senescence columns of the samples of one `direct` dataset."""
    parts = _varying_parts(samples)
    for sample, (varying, varying_fields, _) in zip(samples, parts):
        stated = _stated_state(varying, varying_fields)
        inducer, wording = _first(INDUCERS, varying)
        if sample.get("two_channel"):
            sample.update(senescence_state="unclear", senescence_rule=RULE_TWO_CHANNEL)
        elif stated:
            state, kind, evidence = stated
            sample.update(senescence_state=state, control_kind=kind, senescence_evidence=evidence,
                          senescence_rule=RULE_MIXED if state == "unclear" else RULE_STATED)
        elif inducer and not _first(CONTROLS, _REAGENT_CONTROL.sub(" ", varying))[0]:
            sample.update(senescence_state="inducer applied", senescence_evidence=wording,
                          senescence_rule=RULE_INDUCER)
        else:
            sample["senescence_state"] = "unclear"
    has_senescent = any(sample["senescence_state"] == "senescent" for sample in samples)
    for sample, (varying, _, constant) in zip(samples, parts):
        if sample["senescence_state"] == "unclear" and not sample.get("senescence_rule"):
            kind, wording = _first(CONTROLS, _REAGENT_CONTROL.sub(" ", varying))
            # Wording shared by all samples labels them only when nothing else did.
            shared = _stated_state(constant, [], _SENESCENT_SHARED)
            if has_senescent and kind:
                sample.update(senescence_state="not senescent", control_kind=kind, senescence_evidence=wording,
                              senescence_rule=RULE_CONTROL)
            elif shared and shared[0] != "unclear":
                state, kind, evidence = shared
                sample.update(senescence_state=state, control_kind=kind, senescence_evidence=evidence,
                              senescence_rule=RULE_SHARED)
        if sample["senescence_state"] in ("senescent", "inducer applied"):
            sample["senescence_inducer"] = _first(INDUCERS, f"{varying} {constant}")[0]


def _geo_assay(row):
    if row["type"] == "SRA":
        single_cell = "single cell" in row["library_source"]
        return f"{row['library_strategy']}{', single cell' if single_cell else ''}"
    return GEO_ARRAY_ASSAYS.get(row["type"], row["type"])


def _geo_datasets(samples_dir):
    """Yield (series, samples) for every GEO series; a sample is a dict with its title and repository fields."""
    with gzip.open(samples_dir / GEO_SAMPLES, "rt", newline="") as samples_file, \
            gzip.open(samples_dir / GEO_CHARACTERISTICS, "rt", newline="") as characteristics_file:
        # Both tables list the series and their samples in the same order.
        characteristics = itertools.groupby(csv.DictReader(characteristics_file, delimiter="\t"),
                                            key=lambda row: (row["series"], row["sample"]))
        current = next(characteristics, None)
        for series, rows in itertools.groupby(csv.DictReader(samples_file, delimiter="\t"),
                                              key=lambda row: row["series"]):
            samples = []
            for row in rows:
                described = []
                if current and current[0] == (series, row["sample"]):
                    described = [(item["key"], item["value"]) for item in current[1] if item["channel"] == "ch1"]
                    current = next(characteristics, None)
                described.append(("source name", row["source_name_ch1"]))
                samples.append({"sample": row["sample"], "title": row["title"], "assay": _geo_assay(row),
                                "two_channel": row["channel_count"] == "2",
                                "fields": [("organism", row["organism_ch1"])] + described,
                                "label_fields": described})
            yield series, samples
        if current:
            raise SystemExit(f"characteristics of {current[0]} do not follow the order of the sample table")


def _sheet_datasets(sheets_path):
    """Yield (dataset, samples) for every ArrayExpress experiment with a sample sheet."""
    datasets = defaultdict(dict)
    with gzip.open(sheets_path, "rt", newline="") as handle:
        for (dataset, sheet, row), cells in itertools.groupby(
                csv.DictReader(handle, delimiter="\t"), key=lambda cell: (cell["dataset"], cell["sheet"], cell["row"])):
            fields = [(cell["column"], cell["value"]) for cell in cells]
            described = [field for field in fields if _SHEET_SAMPLE_FIELD.match(field[0])]
            named = defaultdict(list)
            for name, value in fields:
                named[_field_name(name)].append(value)
            sample = title = named["source name"][0] if named["source name"] else f"{sheet} row {row}"
            assay = SEPARATOR.join(dict.fromkeys(named["library strategy"] + named["technology type"]))
            entry = datasets[dataset].setdefault(sample, {"sample": sample, "title": title, "assay": assay,
                                                          "fields": [], "label_fields": []})
            # A sample sheet has one row per data file, so a sample can have several rows.
            entry["fields"] += [field for field in fields if field not in entry["fields"]]
            entry["label_fields"] += [field for field in described if field not in entry["label_fields"]]
    for dataset, samples in sorted(datasets.items()):
        yield dataset, list(samples.values())


def build_sample_table(out_dir, inventory_path):
    """Write the sample table, the coverage of its columns and the counts of its senescence states."""
    with inventory_path.open(newline="") as handle:
        inventory = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    coverage = defaultdict(Counter)
    covered_datasets = defaultdict(set)
    states = defaultdict(lambda: [0, set()])
    seen = set()
    with gzip.open(out_dir / SAMPLES_NAME, "wt", newline="") as out_file:
        writer = csv.DictWriter(out_file, COLUMNS, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()

        def write(dataset, level, rows):
            known = inventory.get(dataset, {})
            group = (known.get("repository", "") if level == "dataset" or not dataset.startswith("GSE") else "GEO",
                     known.get("relevance", ""), level)
            covered_datasets[group].add(dataset)
            for row in rows:
                row.update(dataset=dataset, repository=group[0], relevance=group[1], level=level)
                writer.writerow({column: row.get(column, "") for column in COLUMNS})
                coverage[group]["rows"] += 1
                coverage[group].update(column for column in COVERAGE_COLUMNS[5:] if row.get(column))
                if row.get("senescence_state"):
                    entry = states[(group[0], row["senescence_state"], row.get("senescence_rule", ""))]
                    entry[0] += 1
                    entry[1].add(dataset)

        for dataset, samples in itertools.chain(_geo_datasets(out_dir), _sheet_datasets(out_dir / SHEETS_NAME)):
            seen.add(dataset)
            for sample in samples:
                sample.update(_derived(sample["fields"]))
            if inventory.get(dataset, {}).get("relevance") == DIRECT:
                _label_dataset(samples)
            write(dataset, "sample", samples)
        for dataset, known in sorted(inventory.items()):
            if dataset not in seen and int(known["files_downloaded"] or 0) > 0:
                write(dataset, "dataset", [{"title": known["title"], "organism": known["organism"],
                                            "assay": known["assay"]}])

    with (out_dir / COVERAGE_NAME).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COVERAGE_COLUMNS)
        for group, counts in sorted(coverage.items()):
            writer.writerow((*group, counts["rows"], len(covered_datasets[group]),
                             *(counts[column] for column in COVERAGE_COLUMNS[5:])))
    with (out_dir / STATES_NAME).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(STATE_COLUMNS)
        for group, (samples, datasets) in sorted(states.items()):
            writer.writerow((*group, samples, len(datasets)))
    print(f"{sum(counts['rows'] for counts in coverage.values())} rows of {len(seen)} datasets with samples and "
          f"{sum(len(names) for group, names in covered_datasets.items() if group[2] == 'dataset')} datasets without "
          f"a sample sheet written to {out_dir / SAMPLES_NAME}")


# 4. The label inputs.

LABEL_SERIES_FIELDS = ("title", "summary", "overall_design")
# Name in the output, and the column of the GEO sample table it is copied from.
LABEL_SAMPLE_FIELDS = (
    ("title", "title"),
    ("source_name", "source_name_ch1"),
    ("organism", "organism_ch1"),
    ("molecule", "molecule_ch1"),
    ("treatment_protocol", "treatment_protocol_ch1"),
    ("growth_protocol", "growth_protocol_ch1"),
    ("description", "description"),
    ("library_source", "library_source"),
    ("library_selection", "library_selection"),
)
CHANNEL = "ch1"
# The repository of an ArrayExpress experiment in the sample table.
ARRAYEXPRESS_REPOSITORY = "ArrayExpress"
IDF_PATTERN = "*.idf.txt"
# Series field, and the IDF row whose values it holds.
IDF_SERIES_ROWS = (("title", "Investigation Title"), ("summary", "Experiment Description"),
                   ("overall_design", "Experimental Design"))
# Sample field, and the sample sheet column it is copied from.
SHEET_FIELDS = (("organism", "Characteristics[organism]"), ("description", "Description"),
                ("library_source", "Comment[LIBRARY_SOURCE]"), ("library_selection", "Comment[LIBRARY_SELECTION]"))
# Sample field, and the IDF protocol type whose descriptions it holds. GEO has no field for the sample collection
# protocol, in which ArrayExpress submitters often describe the culture and the treatment.
PROTOCOL_FIELDS = (("treatment_protocol", "treatment protocol"), ("growth_protocol", "growth protocol"),
                   ("sample_collection_protocol", "sample collection protocol"))
ARRAYEXPRESS_SAMPLE_FIELDS = ("sample", "characteristics", "title", "source_name", "organism", "molecule",
                              "treatment_protocol", "growth_protocol", "sample_collection_protocol", "description",
                              "library_source", "library_selection")
# Sample sheet columns that describe the sample; they are given as characteristics, with their column names.
_DESCRIBING_COLUMN = re.compile(r"(?:characteristics|factor ?value|unit) ?\[", re.I)
SOURCE_NAME = "sourcename"
EXTRACT_NAME = "extractname"
MATERIAL_TYPE = "materialtype"
PROTOCOL_REF = "protocolref"


def _read_rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _column_key(column):
    """Return a sample sheet column name in one form: lower case, without spaces."""
    return column.replace(" ", "").lower()


def _join(values):
    return SEPARATOR.join(dict.fromkeys(value for value in values if value))


def read_idf(path):
    """Return the rows of an investigation description (IDF): its tag mapped to its values, empty values kept."""
    with path.open(newline="", encoding="utf-8") as handle:
        return {row[0].strip(): [value.strip() for value in row[1:]]
                for row in csv.reader(handle, delimiter="\t") if row and row[0].strip()}


def arrayexpress_series(idf):
    """Return the series fields of an ArrayExpress experiment from its IDF."""
    return {field: _join(idf.get(tag, [])) for field, tag in IDF_SERIES_ROWS}


def _protocols(idf):
    """Return, per protocol accession of the IDF, its type and its description."""
    rows = (idf.get(tag, []) for tag in ("Protocol Name", "Protocol Type", "Protocol Description"))
    return {name: (kind, description) for name, kind, description in itertools.zip_longest(*rows, fillvalue="")
            if name}


def _extract_material_types(row):
    """Return the material types that a sample sheet row gives after its extract name, that is of the extract."""
    keys = [_column_key(column) for column, _ in row]
    if EXTRACT_NAME not in keys:
        return []
    start = keys.index(EXTRACT_NAME)
    return [value for key, (_, value) in zip(keys[start:], row[start:]) if key == MATERIAL_TYPE]


def arrayexpress_sample(sample, rows, idf):
    """Return the label fields of one sample of an ArrayExpress experiment.

    `rows` are the sample sheet rows of the sample (one per data file), each a list of (column, value) in the order of
    the sheet; `idf` is the experiment's IDF as read_idf() returns it.
    """
    cells = [cell for row in rows for cell in row]
    fields = {"sample": sample, "title": sample, "source_name": sample,
              "characteristics": [list(cell) for cell in dict.fromkeys(cells) if _DESCRIBING_COLUMN.match(cell[0])],
              "molecule": _join(value for row in rows for value in _extract_material_types(row))}
    for field, column in SHEET_FIELDS:
        fields[field] = _join(value for name, value in cells if _column_key(name) == _column_key(column))
    protocols = _protocols(idf)
    references = [value for name, value in cells if _column_key(name) == PROTOCOL_REF]
    for field, kind in PROTOCOL_FIELDS:
        fields[field] = _join(protocols[reference][1] for reference in references if protocols[reference][0] == kind)
    return {field: fields[field] for field in ARRAYEXPRESS_SAMPLE_FIELDS}


def read_sheet_rows(path, datasets):
    """Return, per (ArrayExpress dataset, source name), the sample sheet rows of the sample from the table of sample
    sheets; each row is a list of (column, value) in the order of the sheet."""
    rows = defaultdict(list)
    for cell in _read_rows(path):
        if cell["dataset"] in datasets and cell["data_folder"] == ARRAYEXPRESS_FOLDER:
            rows[(cell["dataset"], cell["sheet"], cell["row"])].append((cell["column"], cell["value"]))
    samples = defaultdict(list)
    for (dataset, sheet, number), row in rows.items():
        names = [value for column, value in row if _column_key(column) == SOURCE_NAME]
        if len(names) != 1:
            raise SystemExit(f"{dataset}: row {number} of {sheet} has {len(names)} source names")
        samples[(dataset, names[0])].append(row)
    return samples


def find_idf(folder):
    idfs = sorted(folder.glob(IDF_PATTERN))
    if len(idfs) != 1:
        raise SystemExit(f"{folder}: {len(idfs)} IDF files")
    return idfs[0]


def prepare_label_inputs(dataset_labels, out_dir, arrayexpress_dir):
    """Write the label inputs of the datasets of the table of dataset labels, one JSON line per dataset."""
    datasets = [row["dataset"] for row in _read_rows(dataset_labels)]
    wanted = set(datasets)
    samples = {dataset: {} for dataset in datasets}
    arrayexpress = set()
    for row in _read_rows(out_dir / SAMPLES_NAME):
        if row["dataset"] in wanted and row["level"] == "sample" and row["assay"] in BULK_RNA_SEQ_ASSAYS:
            if row["sample"] in samples[row["dataset"]]:
                raise SystemExit(f"{row['dataset']}: sample {row['sample']} has two rows in the sample table")
            samples[row["dataset"]][row["sample"]] = {"sample": row["sample"], "characteristics": []}
            if row["repository"] == ARRAYEXPRESS_REPOSITORY:
                arrayexpress.add(row["dataset"])
    print(f"PROGRESS 1/5 sample table read: {sum(map(len, samples.values()))} samples of {len(datasets)} datasets",
          flush=True)

    series = {row["series"]: {field: row[field] for field in LABEL_SERIES_FIELDS}
              for row in _read_rows(out_dir / GEO_SERIES) if row["series"] in wanted}
    print(f"PROGRESS 2/5 series fields read: {len(series)} series", flush=True)

    seen = set()
    for row in _read_rows(out_dir / GEO_SAMPLES):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None:
            sample.update({name: row[column] for name, column in LABEL_SAMPLE_FIELDS})
            seen.add((row["series"], row["sample"]))
    print(f"PROGRESS 3/5 sample fields read: {len(seen)} samples", flush=True)

    for row in _read_rows(out_dir / GEO_CHARACTERISTICS):
        sample = samples.get(row["series"], {}).get(row["sample"])
        if sample is not None and row["channel"] == CHANNEL:
            sample["characteristics"].append([row["key"], row["value"]])
    print("PROGRESS 4/5 characteristics read", flush=True)

    sheet_rows = read_sheet_rows(out_dir / SHEETS_NAME, arrayexpress)
    for dataset in sorted(arrayexpress):
        idf = read_idf(find_idf(arrayexpress_dir / dataset))
        series[dataset] = arrayexpress_series(idf)
        for sample in samples[dataset]:
            if (dataset, sample) in sheet_rows:
                samples[dataset][sample] = arrayexpress_sample(sample, sheet_rows[(dataset, sample)], idf)
                seen.add((dataset, sample))
    print(f"PROGRESS 5/5 sample sheets and IDF read: {len(arrayexpress)} ArrayExpress experiments", flush=True)

    missing_series = wanted - set(series)
    missing_samples = {(dataset, sample) for dataset in datasets for sample in samples[dataset]} - seen
    if missing_series or missing_samples:
        raise SystemExit(f"without fields: {len(missing_series)} series, {len(missing_samples)} samples; "
                         f"the first: {sorted(missing_series)[:3]} {sorted(missing_samples)[:3]}")

    with gzip.open(out_dir / LABEL_INPUTS_NAME, "wt") as handle:
        for dataset in datasets:
            record = {"dataset": dataset, "series": series[dataset], "samples": list(samples[dataset].values())}
            handle.write(json.dumps(record, ensure_ascii=False) + "\n")
    print(f"{len(datasets)} datasets written to {out_dir / LABEL_INPUTS_NAME}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--downloads-dir", type=Path, required=True, help="the folder of download.sh")
    parser.add_argument("--inventory", type=Path, required=True,
                        help="table of datasets: accession, repository, relevance, title, organism, assay and "
                             "files_downloaded")
    parser.add_argument("--dataset-labels", type=Path, required=True,
                        help="table with a column `dataset`: the datasets whose label inputs are written")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    harvest_geo_records(args.downloads_dir / GEO_FOLDER, args.out_dir, args.workers)
    harvest_sample_sheets(args.downloads_dir / ARRAYEXPRESS_FOLDER, args.out_dir)
    build_sample_table(args.out_dir, args.inventory)
    prepare_label_inputs(args.dataset_labels, args.out_dir, args.downloads_dir / ARRAYEXPRESS_FOLDER)


if __name__ == "__main__":
    main()
