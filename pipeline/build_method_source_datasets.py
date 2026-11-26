#!/usr/bin/env python3
"""Collect the lists of all methods in one file: the datasets that the authors of a trained model or a gene set used."""

import argparse
from pathlib import Path

import pandas as pd

from accessions import ACCESSION

COLUMNS = ["method", "accession", "role", "source"]
DERIVED_FROM = "derived from"
APPLIED_ONLY = "applied only"
# A laboratory group with a listed accession in one of these roles is in the list; `applied only` is reported only.
COUNTED_ROLES = ("training", "validation or model selection", "test", "unknown role", DERIVED_FROM)


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def list_gene_set_sources(gene_sets):
    """Return one row per gene set and accession named in its `derived_from` text, under the name of the gene set."""
    rows = gene_sets.assign(accession=gene_sets.derived_from.map(lambda text: sorted(set(ACCESSION.findall(text)))))
    rows = rows.explode("accession").dropna(subset="accession")
    return pd.DataFrame({
        "method": rows.gene_set, "accession": rows.accession, "role": DERIVED_FROM,
        "source": rows.source + "; `derived_from` of the gene set table: " + rows.derived_from,
    }).drop_duplicates(["method", "accession"])


def validate_roles(lists):
    unknown = sorted(set(lists.role) - set(COUNTED_ROLES) - {APPLIED_ONLY})
    if unknown:
        raise ValueError(f"roles outside the allowed list: {unknown}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool-datasets", type=Path, required=True, help="tool_training_datasets.tsv")
    parser.add_argument("--gene-sets", type=Path, required=True, help="gene_sets.tsv, with the column `derived_from`")
    parser.add_argument("--other-lists", type=Path, nargs="+", required=True,
                        help="tables with the columns method, accession, role, source (the markeR case study)")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists")

    tools = _read(args.tool_datasets).rename(columns={"tool": "method"})
    lists = pd.concat([tools[COLUMNS], list_gene_set_sources(_read(args.gene_sets)),
                       *(_read(path)[COLUMNS] for path in args.other_lists)], ignore_index=True)
    validate_roles(lists)
    lists.to_csv(args.out, sep="\t", index=False)
    print(f"{len(lists)} rows, {lists.method.nunique()} methods, {lists.accession.nunique()} accessions")


if __name__ == "__main__":
    main()
