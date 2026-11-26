#!/usr/bin/env python3
"""List the direct human bulk RNA-seq datasets of the sample table, with the states the rules of the table found."""

import argparse
from pathlib import Path

import pandas as pd

HUMAN = "Homo sapiens"
# Assay values of the sample table that are bulk RNA-seq: GEO's library strategy, and the library strategy of an
# ArrayExpress sample sheet joined to its technology type.
BULK_RNA_SEQ_ASSAYS = ("RNA-Seq", "RNA-Seq; sequencing assay")
# Library source of a bulk RNA-seq sample in an ArrayExpress sample sheet; single-cell and genomic libraries share the
# assay values above.
BULK_LIBRARY_SOURCE = "TRANSCRIPTOMIC"
# Sample sheet columns in one form, lower case and without spaces.
_SOURCE_NAME = "sourcename"
_LIBRARY_SOURCE = "comment[library_source]"
_SHEET_ROW = ["dataset", "data_folder", "sheet", "row"]


def _read_library_sources(sheets_path):
    """Return the library source of every sample whose sample sheet states one, by dataset and source name."""
    cells = pd.read_csv(sheets_path, sep="\t", dtype=str, keep_default_na=False)
    column = cells.column.str.replace(" ", "").str.lower()
    names = cells[column == _SOURCE_NAME][_SHEET_ROW + ["value"]].rename(columns={"value": "sample"})
    sources = cells[column == _LIBRARY_SOURCE][_SHEET_ROW + ["value"]].rename(columns={"value": "library_source"})
    rows = names.merge(sources, on=_SHEET_ROW, validate="one_to_one")
    # A sample has one sheet row per data file; a sample with two library sources stops the merge in list_datasets.
    return rows[["dataset", "sample", "library_source"]].drop_duplicates()


def list_datasets(samples_path, sheets_path):
    """Return one row per dataset whose most frequent organism is human.

    A sample whose sample sheet states a library source other than `BULK_LIBRARY_SOURCE` is not bulk RNA-seq.
    """
    samples = pd.read_csv(samples_path, sep="\t", dtype=str, keep_default_na=False)
    rows = samples[(samples.relevance == "direct") & (samples.level == "sample")
                   & samples.assay.isin(BULK_RNA_SEQ_ASSAYS)]
    rows = rows.merge(_read_library_sources(sheets_path), on=["dataset", "sample"], how="left",
                      validate="many_to_one")
    rows = rows[rows.library_source.isna() | (rows.library_source == BULK_LIBRARY_SOURCE)]
    grouped = rows.groupby("dataset", sort=True)
    table = pd.DataFrame({
        "sample_rows": grouped.size(),
        "most_frequent_organism": grouped.organism.agg(lambda values: values.value_counts().index[0]),
        "rule_senescent": grouped.senescence_state.agg(lambda values: (values == "senescent").sum()),
        "rule_not_senescent": grouped.senescence_state.agg(lambda values: (values == "not senescent").sum()),
    })
    table = table[table.most_frequent_organism == HUMAN].drop(columns="most_frequent_organism")
    table["both_states_by_rules"] = ((table.rule_senescent > 0) & (table.rule_not_senescent > 0)).map(
        {True: "yes", False: "no"})
    return table.reset_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--sample-sheets", type=Path, required=True, help="the table of sample sheets")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    table = list_datasets(args.samples, args.sample_sheets)
    table.to_csv(args.out, sep="\t", index=False)
    print(f"{len(table)} datasets, {table.sample_rows.sum()} sample rows, "
          f"{(table.both_states_by_rules == 'yes').sum()} with both states by the rules")


if __name__ == "__main__":
    main()
