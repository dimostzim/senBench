#!/usr/bin/env python3
"""Build one count table per dataset (genes by samples; raw read counts, or expected counts after rounding), check it
and record its source; then write the table of all count tables and the table of their columns.

A table is built from GEO's uniform counts, from our counts of the raw reads, or, through the stored layout of the
dataset, from the submitters' supplementary files: delimited text, workbooks and the members of tar archives.
"""

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
import zipfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import pandas as pd
import xlrd

GEO_UNIFORM = "GEO uniform"
COUNTED_FROM_RAW_READS = "counted from raw reads"
AUTHORS_RAW_COUNTS = "authors raw counts"
AUTHORS_EXPECTED_COUNTS = "authors expected counts, rounded"
AUTHORS_NORMALISED_ONLY = "authors normalised values only"
NO_COUNT_TABLE = "none"
MATCH_ACCESSION = "accession"
MATCH_EXACT = "exact"
MATCH_AGENT = "agent"
NO_COLUMN = "no column"
# Outcomes of a sample in the table of column reassignments; `reassigned` is also the `match` of such a column.
KEPT = "kept"
REASSIGNED = "reassigned"
SUBMITTERS_RAW_COUNTS = "submitters' raw counts"
LEFT_OUT = "left out"
ENSEMBL_GENE = "Ensembl gene"
NCBI_GENE_ID = "NCBI gene id"
GENE_SYMBOL = "gene symbol"
OTHER_IDENTIFIER = "other"
ENSEMBL_GENE_PATTERN = r"ENS[A-Z]*G\d{11}(\.\d+)?"
NCBI_GENE_ID_PATTERN = r"\d+"
# The start of an Ensembl or RefSeq transcript accession.
TRANSCRIPT_PATTERN = r"ENS[A-Z]*T\d{11}|[NX][MR]_\d+"
# Rows that HTSeq and STAR write next to the genes: reads that were not assigned to a gene.
UNASSIGNED_READS_PATTERN = (r"__.*|N_(unmapped|multimapping|noFeature|ambiguous)"
                            r"|no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique")
# Up to this share of the rows of a table may lack a gene identifier (they are left out) or repeat one (they are
# kept); a table with more such rows does not hold one count per gene.
MAX_FAULTY_ROW_SHARE = 0.05
UNASSIGNED_ROWS = "rows of reads that were not assigned to a gene"
UNIDENTIFIED_ROWS = "rows without a gene identifier"
CHECK_FAILED = "check failed:"
UNIFORM_TABLE_PRESENT = "present"
SUPPLEMENTARY_FOLDER = "suppl"
# The folder of GEO's uniform count tables in the folder of download.sh.
UNIFORM_FOLDER = "geo_uniform_counts"
GENE_COLUMN = "gene_identifier"
TABLES_NAME = "count_tables.tsv"
COLUMNS_NAME = "count_table_columns.tsv.gz"
UNITS_NAME = "units.tsv"
TABLE_COLUMNS = ("dataset", "count_table_state", "source_file", "source_sha256", "gene_identifier_kind", "genes",
                 "samples_of_dataset", "samples_in_table", "column_match", "note")
COLUMN_COLUMNS = ("dataset", "sample", "file", "member", "sheet", "column", "column_name", "match")
CHUNK_BYTES = 1 << 20
# The `content` of a stored layout. In a layout of values per sample, each source states its `values`.
RAW_READ_COUNTS = "raw read counts"
NORMALISED_VALUES_ONLY = "normalised values only"
VALUES_PER_SAMPLE = "values per sample"
# The kinds of content of a supplementary file or tar member.
TEXT, WORKBOOK, TAR, OTHER = "text", "workbook", "tar", "other"
# GEO's list of the members of a tar; not a data file.
MEMBER_LIST_NAME = "filelist.txt"
SEPARATORS = {"tab": "\t", "comma": ",", "semicolon": ";", "whitespace": r"\s+"}
_GZIP_START = b"\x1f\x8b"
_WORKBOOK_STARTS = (b"PK\x03\x04", b"\xd0\xcf\x11\xe0")
_TAR_MARK = slice(257, 262)
_SNIFF_BYTES = 8192
# The first non-empty lines of a text content, which give its delimiter.
_FIRST_ROWS = 4
_MAX_LINE_CHARS = 50_000_000

csv.field_size_limit(_MAX_LINE_CHARS)


class ContentProblem(Exception):
    """A file or tar member that cannot be read as its kind."""


@dataclass(frozen=True)
class Table:
    """The values of some columns of one file, member or sheet, as text."""
    # Column number (counted in the rows of values, from 1) to the name of the column; empty without a names row.
    column_names: dict
    values: pd.DataFrame


def list_data_files(folder):
    """Return the data files of a supplementary folder, sorted by name; a missing folder has none."""
    if not folder.is_dir():
        return []
    return sorted(path for path in folder.iterdir() if path.is_file() and path.name != MEMBER_LIST_NAME)


def open_content(opener):
    """Return the kind of a file or tar member and a handle at the start of its content, uncompressed."""
    handle = opener()
    try:
        if handle.peek(len(_GZIP_START)).startswith(_GZIP_START):
            handle = io.BufferedReader(gzip.GzipFile(fileobj=handle))
        start = handle.peek(_SNIFF_BYTES)[:_SNIFF_BYTES]
    except (gzip.BadGzipFile, EOFError) as problem:
        raise ContentProblem(f"gzip: {problem}") from None
    if start.startswith(_WORKBOOK_STARTS):
        return WORKBOOK, handle
    if start[_TAR_MARK] == b"ustar":
        return TAR, handle
    if start and b"\x00" not in start:
        return TEXT, handle
    return OTHER, handle


def open_tar(handle):
    """Open the content of a tar for reading single members."""
    try:
        return tarfile.open(fileobj=handle, mode="r:")
    except tarfile.TarError as problem:
        raise ContentProblem(f"tar: {problem}") from None


def _read_lines(handle):
    """Yield (line number, line) for the lines of a text content that are not empty."""
    text = io.TextIOWrapper(handle, encoding="utf-8-sig", errors="replace")
    number = 0
    try:
        while line := text.readline(_MAX_LINE_CHARS):
            number += 1
            if len(line) == _MAX_LINE_CHARS:
                raise ContentProblem(f"line {number} has more than {_MAX_LINE_CHARS} characters")
            if line.strip():
                yield number, line.rstrip("\r\n")
    except (gzip.BadGzipFile, EOFError) as problem:
        raise ContentProblem(f"gzip: {problem}") from None


def _find_delimiter(lines):
    for name in ("tab", "comma", "semicolon"):
        if any(SEPARATORS[name] in line for line in lines):
            return name
    return "whitespace"


def _split_line(line, delimiter):
    if delimiter == "whitespace":
        return line.split()
    return next(csv.reader([line], delimiter=SEPARATORS[delimiter]))


def _read_text_rows(handle):
    """Return the delimiter of a text content and its first non-empty lines as (line number, fields)."""
    lines = [numbered for numbered, _ in zip(_read_lines(handle), range(_FIRST_ROWS))]
    delimiter = _find_delimiter([line for _, line in lines])
    return delimiter, [(number, _split_line(line, delimiter)) for number, line in lines]


def _read_sheet(handle, sheet):
    """Return one sheet of a workbook as text."""
    try:
        frame = pd.read_excel(io.BytesIO(handle.read()), sheet_name=sheet, header=None, dtype=str, na_filter=False)
    except (gzip.BadGzipFile, EOFError, zipfile.BadZipFile, xlrd.XLRDError, ValueError, KeyError) as problem:
        raise ContentProblem(f"workbook: {problem}") from None
    # Some workbooks give a number or a missing value for a cell although text was asked for.
    return frame.map(lambda cell: "" if pd.isna(cell) else str(cell))


def _name_columns(names, width, columns):
    """Return the names of columns of the rows of values; a names row may lack the name of the first column."""
    offset = width - len(names)
    if offset not in (0, 1):
        raise ContentProblem(f"the row of column names has {len(names)} fields, the first row of values {width}")
    return {column: names[column - 1 - offset].strip() if column > offset else "" for column in columns}


def _read_text_table(opener, column_names_row, first_values_row, columns):
    rows = {}
    delimiter = _read_text_rows(open_content(opener)[1])[0]
    for number, line in _read_lines(open_content(opener)[1]):
        if number in (column_names_row, first_values_row):
            rows[number] = _split_line(line, delimiter)
        if number >= first_values_row:
            break
    if first_values_row not in rows or (column_names_row and column_names_row not in rows):
        raise ContentProblem(f"no fields in row {column_names_row} or row {first_values_row}")
    names = _name_columns(rows[column_names_row], len(rows[first_values_row]), columns) if column_names_row else {}
    try:
        values = pd.read_csv(open_content(opener)[1], sep=SEPARATORS[delimiter], header=None,
                             skiprows=first_values_row - 1, usecols=[column - 1 for column in columns], dtype=str,
                             na_filter=False, skip_blank_lines=False, encoding="utf-8-sig",
                             encoding_errors="replace")
    except (gzip.BadGzipFile, EOFError, ValueError) as problem:
        raise ContentProblem(f"text: {problem}") from None
    return names, values.rename(columns=lambda position: position + 1)


def _read_sheet_table(opener, sheet, column_names_row, first_values_row, columns):
    frame = _read_sheet(open_content(opener)[1], sheet)
    if first_values_row > len(frame) or max(columns) > frame.shape[1]:
        raise ContentProblem(f"the sheet has {len(frame)} rows and {frame.shape[1]} columns")
    names = ({column: frame.iat[column_names_row - 1, column - 1].strip() for column in columns}
             if column_names_row else {})
    values = frame.iloc[first_values_row - 1:, [column - 1 for column in columns]]
    return names, values.set_axis(columns, axis=1)


def read_table(opener, sheet, column_names_row, first_values_row, columns):
    """Return the names and the values of the given columns of a text content or of a sheet of a workbook.

    Rows and columns are numbered from 1; `column_names_row` is 0 where there is no row of column names.
    Rows in which every given column is empty are left out.
    """
    kind = open_content(opener)[0]
    if kind == TEXT and not sheet:
        names, values = _read_text_table(opener, column_names_row, first_values_row, columns)
    elif kind == WORKBOOK and sheet:
        names, values = _read_sheet_table(opener, sheet, column_names_row, first_values_row, columns)
    else:
        raise ContentProblem(f"the content is of kind {kind} and the sheet is {sheet!r}")
    return Table(names, values[(values != "").any(axis=1)].reset_index(drop=True))


class _CheckProblem(Exception):
    """A table that is not one count or one expected count per gene and sample."""


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _check_expected_counts(genes, numbers, where):
    """Raise unless values that are not all whole numbers are expected counts per gene.

    Expected counts keep whole numbers for the genes whose reads are assigned without ambiguity; values that were
    divided by a factor or transformed do not.
    """
    transcripts = genes.str.match(TRANSCRIPT_PATTERN).sum()
    if transcripts > MAX_FAULTY_ROW_SHARE * len(genes):
        raise _CheckProblem(f"{where}: values are not whole numbers, and {transcripts} of {len(genes)} rows are "
                            f"named by a transcript accession")
    most_common = numbers[numbers > 0].mode().iloc[0].dropna()
    scaled = most_common[most_common % 1 != 0]
    if not scaled.empty:
        raise _CheckProblem(f"{where}: values are not whole numbers, and the most common value above zero of "
                            f"{scaled.index[0]} is {scaled.iat[0]}, not a whole number")


def _check_counts(frame, where):
    """Return the rows of genes of a frame with whole numbers as counts, the numbers of rows left out by kind, and
    the number of values that were not whole numbers and were rounded half to even.

    Raises where a value is neither a count nor an expected count.
    """
    genes = frame[GENE_COLUMN].str.strip()
    unassigned = genes.str.fullmatch(UNASSIGNED_READS_PATTERN)
    unidentified = genes == ""
    if unidentified.sum() > MAX_FAULTY_ROW_SHARE * len(genes):
        raise _CheckProblem(f"{where}: {unidentified.sum()} of {len(genes)} rows have no gene identifier")
    left_out = collections.Counter({UNASSIGNED_ROWS: int(unassigned.sum()),
                                    UNIDENTIFIED_ROWS: int(unidentified.sum())})
    of_genes = ~(unassigned | unidentified)
    genes = genes[of_genes].reset_index(drop=True)
    text = frame[of_genes].drop(columns=GENE_COLUMN).reset_index(drop=True)
    if text.empty:
        raise _CheckProblem(f"{where}: no rows of values")
    numbers = text.apply(pd.to_numeric, errors="coerce")
    for problem, wrong in (("are not numbers", numbers.isna()), ("are negative", numbers < 0)):
        if wrong.any(axis=None):
            sample = wrong.any().idxmax()
            raise _CheckProblem(f"{where}: {wrong.sum().sum()} values {problem}, the first "
                                f"{text[sample][wrong[sample]].iat[0]!r} of {sample}")
    # pd.to_numeric can miss a number by its last bit, which decides the rounding of a value next to a half.
    numbers = text.astype("float64")
    not_whole = int((numbers % 1 != 0).sum().sum())
    if not_whole:
        _check_expected_counts(genes, numbers, where)
    counts = numbers.round().astype("int64")
    counts.insert(0, GENE_COLUMN, genes)
    return counts, left_out, not_whole


def _join_frames(frames):
    """Join checked frames by gene identifier; return the table and the number of genes that a frame lacks."""
    first = frames[0][GENE_COLUMN]
    if all(first.equals(frame[GENE_COLUMN]) for frame in frames):
        return pd.concat([frames[0]] + [frame.drop(columns=GENE_COLUMN) for frame in frames[1:]], axis=1), 0
    for frame in frames:
        if frame[GENE_COLUMN].duplicated().any():
            raise _CheckProblem("the sources hold different genes and one of them repeats a gene identifier")
    joined = pd.concat([frame.set_index(GENE_COLUMN) for frame in frames], axis=1, sort=False).astype("Int64")
    lacking = int(joined.isna().any(axis=1).sum())
    return joined.rename_axis(GENE_COLUMN).reset_index(), lacking


def _classify_identifiers(genes, symbols):
    """Return the kind of most gene identifiers, and how many are of that kind."""
    counts = {ENSEMBL_GENE: genes.str.fullmatch(ENSEMBL_GENE_PATTERN).sum(),
              NCBI_GENE_ID: genes.str.fullmatch(NCBI_GENE_ID_PATTERN).sum(), GENE_SYMBOL: genes.isin(symbols).sum()}
    kind = max(counts, key=counts.get)
    if counts[kind] * 2 > len(genes):
        return kind, f"{counts[kind]} of {len(genes)} gene identifiers are of the kind {kind}"
    return OTHER_IDENTIFIER, f"no kind covers more than half of the gene identifiers ({counts})"


def match_column(sample, title, member, column_name):
    """Return how a column is known to belong to a sample."""
    if column_name == sample or re.match(rf"{sample}(?!\d)", member):
        return MATCH_ACCESSION
    if column_name and column_name == title.strip():
        return MATCH_EXACT
    return MATCH_AGENT


def _read_source(folder, source, titles):
    """Return the checked frames of one source of a layout, the rows left out, the number of values rounded, and
    one row per sample column."""
    path = folder / source["file"]
    opener = partial(path.open, "rb")
    kind, handle = open_content(opener)
    tar = open_tar(handle) if kind == TAR else None
    by_member = {}
    for item in source["sample_columns"]:
        by_member.setdefault(item["member"], []).append(item)
    frames = []
    left_out = collections.Counter()
    rounded = 0
    matches = []
    for member, items in by_member.items():
        where = member or source["file"]
        table = read_table(partial(tar.extractfile, member) if member else opener, source["sheet"],
                           source["column_names_row"], source["first_values_row"],
                           [source["gene_identifier_column"]] + [item["column"] for item in items])
        names = {source["gene_identifier_column"]: GENE_COLUMN, **{item["column"]: item["sample"] for item in items}}
        frame, member_left_out, member_rounded = _check_counts(table.values.rename(columns=names), where)
        frames.append(frame)
        left_out += member_left_out
        rounded += member_rounded
        for item in items:
            name = table.column_names.get(item["column"], "")
            matches.append({"sample": item["sample"], "file": source["file"], "member": member,
                            "sheet": source["sheet"], "column": item["column"], "column_name": name,
                            "match": match_column(item["sample"], titles[item["sample"]], member, name)})
    return frames, left_out, rounded, matches


def _read_uniform_table(path, titles, source_columns):
    """Return the checked frame of a GEO uniform table, the rows left out, one row per column of a sample of the
    dataset, and the number of columns that are not samples of the dataset.

    `source_columns` maps a sample of the dataset to the sample accession of the column it takes; a sample that it
    does not name has no column.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str, na_filter=False)
    samples = [sample for sample in titles if source_columns.get(sample) in frame.columns]
    matches = [{"sample": sample, "file": path.name, "member": "", "sheet": "",
                "column": frame.columns.get_loc(source_columns[sample]) + 1, "column_name": source_columns[sample],
                "match": MATCH_ACCESSION if source_columns[sample] == sample else REASSIGNED}
               for sample in samples]
    others = sum(column not in titles for column in frame.columns[1:])
    frame = frame[[frame.columns[0]] + [source_columns[sample] for sample in samples]]
    frame.columns = [GENE_COLUMN] + samples
    frame, left_out, rounded = _check_counts(frame, path.name)
    if rounded:
        raise _CheckProblem(f"{path.name}: {rounded} values are not whole numbers")
    return frame, left_out, matches, others


def _summarise_matches(matches):
    counts = pd.Series([match["match"] for match in matches]).value_counts()
    return "; ".join(f"{kind}: {count}" for kind, count in counts.items())


def _write_table(frames, left_out, matches, titles, symbols, target):
    """Check the joined table, write it, and return the table fields that describe it."""
    table, lacking = _join_frames(frames)
    genes = table[GENE_COLUMN]
    repeated = int(genes.duplicated(keep=False).sum())
    if repeated > MAX_FAULTY_ROW_SHARE * len(genes):
        raise _CheckProblem(f"{repeated} of {len(genes)} rows repeat a gene identifier")
    kind, kind_note = _classify_identifiers(genes, symbols)
    table = table[[GENE_COLUMN] + [sample for sample in titles if sample in table.columns]]
    partial_target = target.with_name(target.name + ".partial")
    table.to_csv(partial_target, sep="\t", index=False, compression="gzip")
    partial_target.replace(target)
    notes = [kind_note] + [f"{count} {kind} were left out, summed over the files and members read"
                           for kind, count in left_out.items() if count]
    if repeated:
        notes.append(f"{repeated} rows repeat a gene identifier")
    if lacking:
        notes.append(f"{lacking} genes have no count in some samples, because the sources hold different genes")
    if len(table.columns) - 1 < len(titles):
        notes.append(f"{len(titles) - len(table.columns) + 1} samples have no column")
    return {"gene_identifier_kind": kind, "genes": len(genes), "samples_in_table": len(table.columns) - 1,
            "column_match": _summarise_matches(matches), "note": "; ".join(notes)}


def _read_layout(record):
    """Return the agent's layout from its stored record; a layout of values per sample whose sources all hold raw read
    counts is taken as a layout of raw read counts."""
    layout = record["answer"]
    if layout["content"] == VALUES_PER_SAMPLE and all(source["values"] == RAW_READ_COUNTS
                                                       for source in layout["sources"]):
        return {**layout, "content": RAW_READ_COUNTS}
    return layout


def build_dataset(dataset, titles, uniform_table, source_columns, suppl_folder, layout_record, symbols, out_dir,
                  uniform_state=GEO_UNIFORM):
    """Build the count table of one dataset; return its row of the table of count tables and its column rows.

    `titles` maps every sample of the dataset to its title; `uniform_table` is a table in the layout of GEO uniform
    tables, None where there is none or where the table is built from the submitters' raw counts instead, and
    `uniform_state` is its `count_table_state`. `source_columns` maps a sample to the accession of the column of that
    table that it takes (a sample not named has no column). `layout_record` is the stored layout of the dataset.
    """
    row = dict.fromkeys(TABLE_COLUMNS, "")
    row.update({"dataset": dataset, "samples_of_dataset": len(titles), "samples_in_table": 0})
    target = out_dir / f"{dataset}.tsv.gz"
    if uniform_table is not None:
        row.update({"count_table_state": uniform_state, "source_file": uniform_table.name,
                    "source_sha256": _sha256(uniform_table)})
        if not source_columns:
            row["note"] = "no table: every sample is left out by the column reassignments"
            return row, []
        frame, left_out, matches, others = _read_uniform_table(uniform_table, titles, source_columns)
        row.update(_write_table([frame], left_out, matches, titles, symbols, target))
        if others:
            row["note"] += f"; {others} columns of the source are not samples of the dataset"
        reassigned = sum(match["match"] == REASSIGNED for match in matches)
        if reassigned or len(source_columns) < len(titles):
            row["note"] += (f"; by the column reassignments, {reassigned} samples take the column of another "
                            f"accession and {len(titles) - len(source_columns)} samples are left out")
        return row, matches
    if not list_data_files(suppl_folder):
        row.update({"count_table_state": NO_COUNT_TABLE, "note": "no supplementary file"})
        return row, []
    layout = _read_layout(layout_record)
    agent_note = f"agent: {layout['note']}"
    if layout["content"] != RAW_READ_COUNTS:
        state = AUTHORS_NORMALISED_ONLY if layout["content"] == NORMALISED_VALUES_ONLY else NO_COUNT_TABLE
        row.update({"count_table_state": state, "note": agent_note})
        return row, []
    if not layout["sources"]:
        row.update({"count_table_state": AUTHORS_RAW_COUNTS,
                    "note": f"no table: no column was matched to a sample; {agent_note}"})
        return row, []
    files = list(dict.fromkeys(source["file"] for source in layout["sources"]))
    row["source_file"] = ";".join(files)
    try:
        frames, left_out, rounded, matches = [], collections.Counter(), 0, []
        for source in layout["sources"]:
            source_frames, source_left_out, source_rounded, source_matches = _read_source(suppl_folder, source,
                                                                                          titles)
            frames += source_frames
            left_out += source_left_out
            rounded += source_rounded
            matches += source_matches
        described = _write_table(frames, left_out, matches, titles, symbols, target)
    except (_CheckProblem, ContentProblem) as problem:
        row.update({"count_table_state": NO_COUNT_TABLE, "note": f"{CHECK_FAILED} {problem}; {agent_note}"})
        return row, []
    row.update({"count_table_state": AUTHORS_EXPECTED_COUNTS if rounded else AUTHORS_RAW_COUNTS,
                "source_sha256": ";".join(_sha256(suppl_folder / name) for name in files)})
    row.update(described)
    if rounded:
        values = sum(len(frame) * (frame.shape[1] - 1) for frame in frames)
        row["note"] += (f"; {rounded} of {values} values ({rounded / values:.1%}) were not whole numbers and were "
                        f"rounded half to even")
    row["note"] += f"; {agent_note}"
    return row, matches


def read_reassignments(path):
    """Return the datasets whose table is built from the submitters' raw counts, and per other dataset of the table
    of column reassignments the accession of the GEO uniform column that each sample with a column takes."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    unknown = set(table.outcome) - {KEPT, REASSIGNED, SUBMITTERS_RAW_COUNTS, LEFT_OUT}
    if unknown:
        raise ValueError(f"{path}: unknown outcomes {sorted(unknown)}")
    from_submitters = set(table.dataset[table.outcome == SUBMITTERS_RAW_COUNTS])
    uniform = table[~table.dataset.isin(from_submitters)]
    taken = uniform[uniform.outcome.isin([KEPT, REASSIGNED])]
    repeated = taken[taken.duplicated(["dataset", "count_source_column"], keep=False)]
    if not repeated.empty:
        raise ValueError(f"{path}: a GEO column is given to several samples: {repeated[['dataset', 'sample']].values}")
    source_columns = {dataset: {} for dataset in uniform.dataset}
    for row in taken.itertuples():
        source_columns[row.dataset][row.sample] = row.count_source_column
    return from_submitters, source_columns


def _build_unit(datasets_titles, repositories, uniform_files, downloads_dir, layouts, symbols, out_dir, from_submitters,
                reassigned_columns, counted_dir, dataset):
    """Build one dataset; return (dataset, status, note) for the units table and its record (None if it failed)."""
    suppl_folder = downloads_dir / repositories[dataset] / dataset / SUPPLEMENTARY_FOLDER
    uses_uniform = dataset in uniform_files and dataset not in from_submitters
    uniform_table = downloads_dir / UNIFORM_FOLDER / dataset / uniform_files[dataset] if uses_uniform else None
    uniform_state = GEO_UNIFORM
    if (counted_dir / f"{dataset}.tsv.gz").exists():
        uniform_table, uniform_state = counted_dir / f"{dataset}.tsv.gz", COUNTED_FROM_RAW_READS
    if uniform_table is None and list_data_files(suppl_folder) and dataset not in layouts:
        return dataset, "failed", "no layout from the agent", None
    titles = datasets_titles[dataset]
    source_columns = reassigned_columns.get(dataset, {sample: sample for sample in titles})
    row, matches = build_dataset(dataset, titles, uniform_table, source_columns, suppl_folder, layouts.get(dataset),
                                 symbols, out_dir, uniform_state)
    matched = {match["sample"]: match for match in matches}
    columns = [{"dataset": dataset, **matched.get(sample, {"sample": sample, "match": NO_COLUMN})} for sample in titles]
    note = f"{row['count_table_state']}; {row['samples_in_table']} of {row['samples_of_dataset']} samples"
    return dataset, "ok", note, {"table": row, "columns": columns}


def write_summary_tables(records, out_dir):
    """Write the table of count tables and the table of their columns from the records of all datasets, in the order
    of the datasets."""
    records = sorted(records, key=lambda record: record["table"]["dataset"])
    tables = pd.DataFrame([record["table"] for record in records], columns=TABLE_COLUMNS)
    tables.to_csv(out_dir / TABLES_NAME, sep="\t", index=False)
    columns = pd.DataFrame([column for record in records for column in record["columns"]], columns=COLUMN_COLUMNS)
    columns["column"] = columns["column"].astype("Int64")
    columns.to_csv(out_dir / COLUMNS_NAME, sep="\t", index=False, compression="gzip")
    return tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", type=Path, required=True,
                        help="table with the columns `dataset` and `repository`, the folder of --downloads-dir that "
                             "holds the dataset's folder")
    parser.add_argument("--uniform-state", type=Path, required=True,
                        help="table with the columns `dataset`, `geo_uniform_count_table` and `url`")
    parser.add_argument("--label-inputs", type=Path, required=True, help="label_inputs.jsonl.gz of the sample table")
    parser.add_argument("--downloads-dir", type=Path, required=True, help="the folder of download.sh")
    parser.add_argument("--annotation", type=Path, required=True, help="GEO's human gene annotation table")
    parser.add_argument("--layouts", type=Path, required=True, help="JSON: the stored layout of each dataset")
    parser.add_argument("--records-not-rebuilt", type=Path, required=True,
                        help="JSON: the records of the datasets whose tables are not built, for the summary tables")
    parser.add_argument("--column-reassignments", type=Path, required=True,
                        help="its datasets take their GEO uniform columns as it assigns them, or are built from the "
                             "submitters' raw counts with their stored layout")
    parser.add_argument("--counted-dir", type=Path, required=True,
                        help="our count tables of raw reads, `<dataset>.tsv.gz`: a dataset with one is built from it, "
                             "with the state `counted from raw reads`")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()

    with args.datasets.open(newline="") as handle:
        repositories = {row["dataset"]: row["repository"] for row in csv.DictReader(handle, delimiter="\t")}
    datasets = list(repositories)
    with args.uniform_state.open(newline="") as handle:
        uniform_files = {row["dataset"]: row["url"].rsplit("/", 1)[1]
                         for row in csv.DictReader(handle, delimiter="\t")
                         if row["geo_uniform_count_table"] == UNIFORM_TABLE_PRESENT}
    with gzip.open(args.label_inputs, "rt") as handle:
        titles = {record["dataset"]: {sample["sample"]: sample["title"] for sample in record["samples"]}
                  for record in map(json.loads, handle)}
    symbols = set(pd.read_csv(args.annotation, sep="\t", usecols=["Symbol"], dtype=str, na_filter=False)["Symbol"])
    from_submitters, reassigned_columns = read_reassignments(args.column_reassignments)
    layouts = json.loads(args.layouts.read_text())
    records = list(json.loads(args.records_not_rebuilt.read_text()).values())
    args.out_dir.mkdir(parents=True, exist_ok=True)

    build = partial(_build_unit, titles, repositories, uniform_files, args.downloads_dir, layouts, symbols,
                    args.out_dir, from_submitters, reassigned_columns, args.counted_dir)
    failed = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool, \
            (args.out_dir / UNITS_NAME).open("w", newline="") as units_handle:
        for done, (dataset, status, note, record) in enumerate(pool.map(build, datasets), start=1):
            csv.writer(units_handle, delimiter="\t", lineterminator="\n").writerow((dataset, status, note))
            print(f"PROGRESS {done}/{len(datasets)} {dataset} {status}: {note}", flush=True)
            if status == "failed":
                failed.append(dataset)
            else:
                records.append(record)
    tables = write_summary_tables(records, args.out_dir)
    print(tables[tables["dataset"].isin(datasets)]["count_table_state"].value_counts().to_string())
    if failed:
        raise SystemExit(f"failed units ({len(failed)}): {' '.join(failed)}")


if __name__ == "__main__":
    main()
