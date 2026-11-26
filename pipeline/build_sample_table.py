#!/usr/bin/env python3
"""Build one sample table from the copied repository fields, with columns and senescence labels derived by rules.

Every sample gets organism, assay, tissue, cell type, cell line, sex, age, genotype, treatment, time and passage,
taken from the repository fields by name. Samples of datasets with relevance `direct` also get a senescence state,
an inducer and a kind of control, derived from the wording of their own fields. A dataset without a sample sheet
gets one row that describes the dataset.
"""

import argparse
import csv
import gzip
import itertools
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

GEO_SAMPLES = "geo_samples.tsv.gz"
GEO_CHARACTERISTICS = "geo_sample_characteristics.tsv.gz"
SHEETS_PATTERN = "sample_sheets_*.tsv.gz"
SEPARATOR = "; "
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
SENNET_ORGANISMS = {"human": "Homo sapiens", "mouse": "Mus musculus"}
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

csv.field_size_limit(sys.maxsize)


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


def _sheet_datasets(samples_dir):
    """Yield (dataset, samples) for every dataset with a sample sheet outside GEO."""
    datasets = defaultdict(dict)
    for path in sorted(samples_dir.glob(SHEETS_PATTERN)):
        with gzip.open(path, "rt", newline="") as handle:
            for (dataset, folder, sheet, row), cells in itertools.groupby(
                    csv.DictReader(handle, delimiter="\t"),
                    key=lambda cell: (cell["dataset"], cell["data_folder"], cell["sheet"], cell["row"])):
                fields = [(cell["column"], cell["value"]) for cell in cells]
                if folder == "sennet":
                    # The metadata tables of SenNet describe the assay run, not the sample; the record does both.
                    if not sheet.endswith(".json"):
                        continue
                    exact = dict(reversed(fields))
                    sample = exact.get("sennet_id", sheet)
                    title = exact.get("title", "")
                    assay = exact.get("dataset_type", "")
                    kind = next((value.lower() for name, value in fields if name.endswith("source_type")), "")
                    fields.append(("organism", SENNET_ORGANISMS.get(kind, kind)))
                    # A SenNet record describes the project in free text; only its title describes the sample.
                    described = []
                else:
                    described = [field for field in fields if _SHEET_SAMPLE_FIELD.match(field[0])]
                    named = defaultdict(list)
                    for name, value in fields:
                        named[_field_name(name)].append(value)
                    sample = title = named["source name"][0] if named["source name"] else f"{sheet} row {row}"
                    assay = SEPARATOR.join(dict.fromkeys(named["library strategy"] + named["technology type"]))
                    if folder == "proteomics_metabolomics":
                        assay = "mass spectrometry proteomics"
                entry = datasets[dataset].setdefault(sample, {"sample": sample, "title": title, "assay": assay,
                                                              "fields": [], "label_fields": []})
                # A sample sheet has one row per data file, so a sample can have several rows.
                entry["fields"] += [field for field in fields if field not in entry["fields"]]
                entry["label_fields"] += [field for field in described if field not in entry["label_fields"]]
    for dataset, samples in sorted(datasets.items()):
        yield dataset, list(samples.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-dir", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--states", type=Path, required=True)
    args = parser.parse_args()

    with args.inventory.open(newline="") as handle:
        inventory = {row["accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    coverage = defaultdict(Counter)
    covered_datasets = defaultdict(set)
    states = defaultdict(lambda: [0, set()])
    seen = set()
    with gzip.open(args.out, "wt", newline="") as out_file:
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

        for dataset, samples in itertools.chain(_geo_datasets(args.samples_dir), _sheet_datasets(args.samples_dir)):
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

    with args.coverage.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COVERAGE_COLUMNS)
        for group, counts in sorted(coverage.items()):
            writer.writerow((*group, counts["rows"], len(covered_datasets[group]),
                             *(counts[column] for column in COVERAGE_COLUMNS[5:])))
    with args.states.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(STATE_COLUMNS)
        for group, (samples, datasets) in sorted(states.items()):
            writer.writerow((*group, samples, len(datasets)))
    print(f"{sum(counts['rows'] for counts in coverage.values())} rows of {len(seen)} datasets with samples and "
          f"{sum(len(names) for group, names in covered_datasets.items() if group[2] == 'dataset')} datasets without "
          f"a sample sheet written to {args.out}")


if __name__ == "__main__":
    main()
