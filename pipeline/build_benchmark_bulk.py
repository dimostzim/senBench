#!/usr/bin/env python3
"""Build the bulk benchmark table: labelled samples of cultured human cells with raw counts, one count table per dataset."""

import argparse
import re
from pathlib import Path

import pandas as pd

from build_count_tables import (AUTHORS_EXPECTED_COUNTS, AUTHORS_RAW_COUNTS, COUNTED_FROM_RAW_READS, GEO_UNIFORM,
                                NO_COLUMN)
from build_sample_table import BULK_RNA_SEQ_ASSAYS

HUMAN = "Homo sapiens"
STATES = ("senescent", "not senescent")
RAW_COUNT_STATES = (GEO_UNIFORM, AUTHORS_RAW_COUNTS, AUTHORS_EXPECTED_COUNTS, COUNTED_FROM_RAW_READS)
LABEL_COLUMNS = ("state", "state_basis", "inducer", "control_kind", "cell_type", "cell_line", "time_after_induction",
                 "other_treatment")
ENSEMBL_GENE = re.compile(r"ENSG\d+")
NCBI_GENE = re.compile(r"\d+")
SAMPLES_NAME = "benchmark_bulk_samples.tsv"
DATASETS_NAME = "benchmark_bulk_datasets.tsv"
LEFT_OUT_NAME = "benchmark_bulk_left_out.tsv"
COUNTS_NAME = "benchmark_bulk_counts.tsv.gz"
COUNTS_FOLDER = "counts"
GENOMIC_DNA_REASON = "reads on nearly every gene, as with genomic DNA in the libraries"
GEO_SAMPLE = re.compile(r"GSM\d+")


def sample_key(dataset, sample):
    """Return the key that names one sample across datasets.

    A GEO sample accession names the same sample in every series that holds it, such as a SuperSeries and its
    SubSeries. Any other sample name, such as the source name of an ArrayExpress sample sheet, names a sample of its own
    dataset only, so its key is the dataset and the sample.
    """
    return sample if GEO_SAMPLE.fullmatch(sample) else f"{dataset}/{sample}"


def sample_keys(datasets, samples):
    """Return sample_key() of each pair of a dataset and a sample."""
    return [sample_key(dataset, sample) for dataset, sample in zip(datasets, samples, strict=True)]


def _read(path, **options):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **options)


def select_samples(labels, dataset_labels, organisms, count_tables, count_columns):
    """Return the labelled human samples that have a column of raw counts, and the samples left out with the reason."""
    table = labels[["dataset", "sample", *LABEL_COLUMNS]].merge(organisms, on=["dataset", "sample"], how="left",
                                                                 validate="1:1")
    if table.organism.isna().any():
        missing = table.loc[table.organism.isna(), ["dataset", "sample"]].values.tolist()
        raise ValueError(f"{len(missing)} labelled samples have no bulk RNA-seq row in the sample table: {missing[:5]}")
    table = table.merge(dataset_labels[["dataset", "cultured_cells"]], on="dataset", validate="m:1")
    table = table.merge(count_tables[["dataset", "count_table_state", "gene_identifier_kind"]], on="dataset",
                        how="left", validate="m:1")
    table = table.merge(count_columns[["dataset", "sample", "match"]].rename(columns={"match": "column_match"}),
                        on=["dataset", "sample"], how="left", validate="1:1")
    reasons = pd.Series("", index=table.index)
    for reason, applies in (
        ("the state is not `senescent` or `not senescent`", ~table.state.isin(STATES)),
        ("the organism is not human", table.organism != HUMAN),
        ("`cultured_cells` of the dataset is not `yes`", table.cultured_cells != "yes"),
        ("the dataset has no count table yet", table.count_table_state.isna()),
        ("the dataset has no raw counts or rounded expected counts", ~table.count_table_state.isin(RAW_COUNT_STATES)),
        ("the count table has no column for the sample", table.column_match.isna() | (table.column_match == NO_COLUMN)),
    ):
        reasons = reasons.mask((reasons == "") & applies.fillna(False), reason)
    return table[reasons == ""].copy(), table.assign(reason=reasons)[reasons != ""]


def assign_samples_to_datasets(eligible):
    """Keep every sample once: the dataset with the most eligible samples takes its samples first.

    A sample of a SuperSeries is also a sample of a SubSeries. A sample whose state is `senescent` in one dataset
    and `not senescent` in another is left out. A dataset stays only if it keeps samples of both states.
    """
    keys = pd.Series(sample_keys(eligible.dataset, eligible["sample"]), index=eligible.index)
    states_per_sample = eligible.groupby(keys).state.nunique()
    contradicted = keys.map(states_per_sample) > 1
    left_out = [eligible[contradicted].assign(reason="the state differs between the datasets that hold the sample")]
    eligible, keys = eligible[~contradicted], keys[~contradicted]
    sizes = eligible.groupby("dataset").size()
    order = sorted(sizes.index, key=lambda dataset: (-sizes[dataset], dataset))
    taken = set()
    kept = []
    for dataset in order:
        rows = eligible[(eligible.dataset == dataset) & ~keys.isin(taken)]
        if set(rows.state) == set(STATES):
            kept.append(rows)
            taken.update(keys[rows.index])
    kept = pd.concat(kept) if kept else eligible.iloc[:0]
    rest = eligible[~keys.isin(taken)]
    left_out.append(rest.assign(reason="no dataset keeps the sample together with samples of both states"))
    return kept, pd.concat(left_out)


def name_studies(all_labels, kept):
    """Give every kept dataset a study: datasets that share a sample, directly or through other datasets, are one study."""
    parent = {dataset: dataset for dataset in all_labels.dataset.unique()}

    def find(dataset):
        while parent[dataset] != dataset:
            parent[dataset] = parent[parent[dataset]]
            dataset = parent[dataset]
        return dataset

    keys = pd.Series(sample_keys(all_labels.dataset, all_labels["sample"]), index=all_labels.index)
    for datasets in all_labels.groupby(keys).dataset.agg(list):
        for other in datasets[1:]:
            parent[find(other)] = find(datasets[0])
    members = {}
    for dataset in parent:
        members.setdefault(find(dataset), []).append(dataset)
    sizes = all_labels.groupby("dataset").size()
    # A study is named after its dataset with the most samples, which is the SuperSeries where there is one.
    names = {root: max(datasets, key=lambda dataset: (sizes[dataset], dataset)) for root, datasets in members.items()}
    return kept.dataset.map(lambda dataset: names[find(dataset)])


def map_to_ncbi_gene_ids(counts, annotation):
    """Return the counts with NCBI gene ids as rows; counts of rows that map to the same gene are added."""
    by_ensembl = annotation[annotation.EnsemblGeneID != ""].drop_duplicates("EnsemblGeneID").set_index(
        "EnsemblGeneID").GeneID
    by_symbol = annotation.drop_duplicates("Symbol").set_index("Symbol").GeneID
    known = set(annotation.GeneID)

    def gene_id(identifier):
        ensembl = ENSEMBL_GENE.match(identifier)
        if ensembl:
            return by_ensembl.get(ensembl.group())
        if NCBI_GENE.fullmatch(identifier):
            return identifier if identifier in known else None
        return by_symbol.get(identifier)

    mapped = counts.index.map(gene_id)
    found = pd.notna(mapped)
    result = counts[found].groupby(mapped[found]).sum()
    result.index.name = "GeneID"
    return result, int((~found).sum())


def build_dataset_counts(dataset, samples, counts_dir, annotation):
    """Return the count table of the kept samples of one dataset with NCBI gene ids, and what was left out of it."""
    table = _read(counts_dir / f"{dataset}.tsv.gz")
    table = table.set_index(table.columns[0])[list(samples)]
    complete = (table != "").all(axis=1)
    counts, not_mapped = map_to_ncbi_gene_ids(table[complete].astype("int64"), annotation)
    return counts, {"dataset": dataset, "genes": len(counts), "rows_with_an_empty_cell_left_out": int((~complete).sum()),
                    "rows_without_an_ncbi_gene_id_left_out": not_mapped}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-labels", type=Path, required=True)
    parser.add_argument("--dataset-labels", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True, help="the sample table, for the organism")
    parser.add_argument("--count-tables", type=Path, required=True)
    parser.add_argument("--count-columns", type=Path, required=True)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True, help="GEO's gene annotation table")
    parser.add_argument("--min-genes", type=int, required=True,
                        help="a dataset whose count table maps to fewer NCBI gene ids is left out")
    parser.add_argument("--max-median-genes-with-reads", type=int, required=True,
                        help="a dataset whose median sample has reads on more genes is left out")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists")

    labels = _read(args.sample_labels)
    organisms = _read(args.samples, usecols=["dataset", "level", "sample", "organism", "assay"])
    organisms = organisms[(organisms.level == "sample") & organisms.assay.isin(BULK_RNA_SEQ_ASSAYS)][
        ["dataset", "sample", "organism"]]
    eligible, left_out = select_samples(labels, _read(args.dataset_labels), organisms, _read(args.count_tables),
                                        _read(args.count_columns))
    kept, not_assigned = assign_samples_to_datasets(eligible)
    kept = kept.assign(study=name_studies(labels, kept)).sort_values(["study", "dataset", "sample"])
    annotation = _read(args.annotation, usecols=["GeneID", "Symbol", "EnsemblGeneID"])

    (args.out_dir / COUNTS_FOLDER).mkdir(parents=True)
    all_counts, dataset_rows = [], []
    too_few_genes, reads_on_nearly_every_gene = [], []
    for dataset, rows in kept.groupby("dataset"):
        counts, row = build_dataset_counts(dataset, rows["sample"], args.counts_dir, annotation)
        if len(counts) < max(args.min_genes, 1):
            too_few_genes.append(dataset)
            continue
        if (counts > 0).sum().median() > args.max_median_genes_with_reads:
            reads_on_nearly_every_gene.append(dataset)
            continue
        counts.to_csv(args.out_dir / COUNTS_FOLDER / f"{dataset}.tsv.gz", sep="\t")
        all_counts.append(counts.set_axis(sample_keys([dataset] * counts.shape[1], counts.columns), axis=1))
        dataset_rows.append({**row, "study": rows.study.iloc[0], "count_table_state": rows.count_table_state.iloc[0],
                             "senescent": int((rows.state == "senescent").sum()),
                             "not_senescent": int((rows.state == "not senescent").sum())})
        print(f"PROGRESS {len(dataset_rows)}/{kept.dataset.nunique()} {dataset}: {len(rows)} samples", flush=True)

    few_genes = kept[kept.dataset.isin(too_few_genes)].assign(
        reason=f"the count table maps to fewer than {max(args.min_genes, 1)} NCBI gene ids")
    genomic_dna = kept[kept.dataset.isin(reads_on_nearly_every_gene)].assign(reason=GENOMIC_DNA_REASON)
    kept = kept[~kept.dataset.isin(too_few_genes + reads_on_nearly_every_gene)]
    columns = ["sample", "dataset", "study", *LABEL_COLUMNS, "count_table_state", "gene_identifier_kind", "column_match"]
    kept[columns].to_csv(args.out_dir / SAMPLES_NAME, sep="\t", index=False)
    pd.DataFrame(dataset_rows).to_csv(args.out_dir / DATASETS_NAME, sep="\t", index=False)
    pd.concat([left_out, not_assigned, few_genes, genomic_dna])[["dataset", "sample", "state", "reason"]].to_csv(
        args.out_dir / LEFT_OUT_NAME, sep="\t", index=False)
    pd.concat(all_counts, axis=1).astype("Int64").to_csv(args.out_dir / COUNTS_NAME, sep="\t")
    print(f"{len(kept)} samples ({(kept.state == 'senescent').sum()} senescent, "
          f"{(kept.state == 'not senescent').sum()} not senescent) of {kept.dataset.nunique()} datasets in "
          f"{kept.study.nunique()} studies")


if __name__ == "__main__":
    main()
