"""Check that the dataset list takes a sample whose sample sheet states a library source only if it is transcriptomic,
and that it keeps samples without a sample sheet."""

import pandas as pd
import pytest

from list_bulk_human_datasets import list_datasets

# GSE1 has no sample sheet; E-MTAB-1 is bulk, E-MTAB-2 single cell, E-MTAB-3 holds one bulk and one genomic library.
SAMPLES = pd.DataFrame({
    "dataset": ["GSE1", "GSE1", "E-MTAB-1", "E-MTAB-1", "E-MTAB-2", "E-MTAB-2", "E-MTAB-3", "E-MTAB-3"],
    "sample": ["GSM1", "GSM2", "A", "B", "A", "B", "A", "B"],
    "relevance": "direct",
    "level": "sample",
    "assay": ["RNA-Seq", "RNA-Seq"] + ["RNA-Seq; sequencing assay"] * 6,
    "organism": "Homo sapiens",
    "senescence_state": ["senescent", "not senescent"] * 4,
})
LIBRARY_SOURCES = [("E-MTAB-1", "A", "TRANSCRIPTOMIC"), ("E-MTAB-1", "B", "TRANSCRIPTOMIC"),
                   ("E-MTAB-2", "A", "TRANSCRIPTOMIC SINGLE CELL"), ("E-MTAB-2", "B", "TRANSCRIPTOMIC SINGLE CELL"),
                   ("E-MTAB-3", "A", "TRANSCRIPTOMIC"), ("E-MTAB-3", "B", "GENOMIC")]


def _write_tables(tmp_path, library_sources):
    """Write the sample table and a table of sample sheets with one sheet row per (dataset, sample, library source)."""
    cells = []
    for row, (dataset, sample, source) in enumerate(library_sources):
        sheet_row = {"dataset": dataset, "data_folder": "biostudies_arrayexpress", "sheet": f"{dataset}.sdrf.txt",
                     "row": str(row)}
        # Submitters write the column names in either case.
        cells += [{**sheet_row, "column": "Source name" if dataset == "E-MTAB-3" else "Source Name", "value": sample},
                  {**sheet_row, "column": "Comment[LIBRARY_SOURCE]", "value": source}]
    samples_path, sheets_path = tmp_path / "samples.tsv", tmp_path / "sheets.tsv"
    SAMPLES.to_csv(samples_path, sep="\t", index=False)
    pd.DataFrame(cells).to_csv(sheets_path, sep="\t", index=False)
    return samples_path, sheets_path


def test_list_keeps_bulk_libraries_and_samples_without_a_sheet(tmp_path):
    # A second data file of a sample repeats its library source.
    table = list_datasets(*_write_tables(tmp_path, LIBRARY_SOURCES + [("E-MTAB-1", "A", "TRANSCRIPTOMIC")]))
    assert dict(zip(table.dataset, table.sample_rows)) == {"E-MTAB-1": 2, "E-MTAB-3": 1, "GSE1": 2}


def test_list_refuses_a_sample_with_two_library_sources(tmp_path):
    with pytest.raises(pd.errors.MergeError):
        list_datasets(*_write_tables(tmp_path, LIBRARY_SOURCES + [("E-MTAB-1", "A", "GENOMIC")]))
