#!/usr/bin/env python3
"""Build a dataset for Agentomics-ML from the final training data and `test` of a round of the bulk benchmark.

The dataset folder holds `train/`, `validation/` and `test/`; Agentomics hides `test/` from the agent and scores the
final model of every run on it. A second folder holds the key from each anonymous id to its sample; it is never given
to the agent.

Every split holds `input/data.csv` (column `id` and one column per gene: the read count as in the count table, not
normalized; `NA` where the sample's count table has no row for the gene) and `labels.csv` (`id`, `label`). The columns
are every gene with a value in at least one sample of the training data (train and validation together; never `test`).
Validation holds a drawn share of the laboratory groups of `development`; train holds the rest of the training data.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from accessions import sample_keys

SEED = 20261124
TRAIN = "train"
VALIDATION = "validation"
TEST = "test"
DEVELOPMENT = "development"
LABELS = {"senescent": "senescent", "not senescent": "not_senescent"}
INPUT_FILE = Path("input") / "data.csv"
LABELS_FILE = "labels.csv"
ID_COLUMN = "id"
LABEL_COLUMN = "label"
KEY_COLUMNS = [ID_COLUMN, "split", "study", "dataset", "sample", "laboratory_group", "part", LABEL_COLUMN]
DESCRIPTION = """# Senescence in bulk RNA-seq of cultured human cells

## Task
Binary classification: for each sample, predict whether the cultured cells are senescent.

## Background
Cellular senescence is a lasting stop of cell division that cells enter after stress or damage, such as many rounds
of division (replicative senescence), DNA-damaging drugs, irradiation or an activated oncogene. Senescent cells stay
alive, and they change their shape and metabolism and secrete inflammatory factors. Cells that have stopped dividing
for a time and can divide again (quiescent cells) are not senescent.

## Samples
- Each sample is the RNA-seq of one culture of human cells from a public study in GEO or ArrayExpress.
- The samples come from hundreds of studies and laboratories. The cells include fibroblasts, cancer cell lines,
  endothelial cells, mesenchymal stromal cells and epithelial cells, and senescence was induced in many ways.
- A study usually holds both senescent samples and their controls.

## Target
- `label` in `labels.csv` is `senescent` or `not_senescent`.
- `senescent`: the cells were made senescent, and the study shows or states that they are senescent at the time
  they were sampled.
- `not_senescent`: a control of the senescent cells of the same study: proliferating, quiescent, early passage or
  young, untreated, vehicle or empty vector.
- The labels were read from the public records and the papers. No person checked them, and some may be wrong.

## Data format
- `input/data.csv` has a column `id`, which matches `id` in `labels.csv`, and one column per gene.
- A gene column is named by the gene's symbol. Where a symbol is shared by several genes, the name is
  `<symbol>_<NCBI Gene id>`; where a gene has no symbol, it is `GeneID_<NCBI Gene id>`.
- The values are the read counts of each gene, whole numbers that are not normalized, so they depend on how deeply
  each sample was sequenced.
- The counts come from different sources: GEO's own counts, the authors' raw counts, the authors' rounded expected
  counts, or counts made from the raw reads.
- A cell is `NA` where the count table of the sample's study has no row for the gene; some tables lack many genes.
  `NA` is not a zero.

## Splits
- Train and validation share no laboratory.
- The validation samples, and the samples on which the final model will be tested, come from laboratories that have
  no sample in train.
"""


def list_samples(training, test):
    """Return one row per labelled sample of the training data and of `test`, with its key and label."""
    training = training[training.state_for_training.isin(LABELS)].assign(label=lambda t: t.state_for_training.map(LABELS))
    test = test[test.state.isin(LABELS)].assign(label=lambda t: t.state.map(LABELS), part=TEST, split=TEST)
    samples = pd.concat([training.assign(split=TRAIN), test], ignore_index=True)
    samples["key"] = sample_keys(samples.dataset, samples["sample"])
    if samples.key.duplicated().any():
        raise ValueError(f"a sample is listed twice: {samples.key[samples.key.duplicated()].tolist()[:5]}")
    return samples[["key", "study", "dataset", "sample", "laboratory_group", "part", "split", LABEL_COLUMN]]


def draw_validation(samples, share, seed):
    """Return the samples with the drawn laboratory groups of `development` moved from train to validation."""
    groups = np.sort(samples.loc[samples.part == DEVELOPMENT, "laboratory_group"].unique())
    drawn = np.random.default_rng(seed).permutation(groups)[:round(share * len(groups))]
    in_validation = (samples.split == TRAIN) & samples.laboratory_group.isin(drawn)
    return samples.assign(split=samples.split.mask(in_validation, VALIDATION))


def validate_splits(samples):
    """Raise when a laboratory group or a study has samples in two splits."""
    for column in ("laboratory_group", "study"):
        splits = samples.groupby(column).split.nunique()
        if (splits > 1).any():
            raise ValueError(f"a {column} is in two splits: {splits[splits > 1].index.tolist()[:5]}")


def choose_genes(counts, keys):
    """Return the genes with a value in at least one of the samples named by `keys`."""
    return counts.index[counts[keys].notna().any(axis=1)]


def validate_counts(counts):
    """Raise when a count is negative or not a whole number."""
    values = counts.to_numpy()
    values = values[~np.isnan(values)]
    if (values < 0).any() or (values % 1 != 0).any():
        raise ValueError("a count is negative or not a whole number")


def name_genes(gene_ids, annotation):
    """Return a column name per gene: its symbol; symbol and gene id where the symbol is not unique; `GeneID_<id>`
    where it has none."""
    symbols = annotation.reindex(gene_ids).fillna("")
    repeated = symbols.duplicated(keep=False)
    return [f"GeneID_{gene}" if not symbol else f"{symbol}_{gene}" if clash else symbol
            for gene, symbol, clash in zip(gene_ids, symbols, repeated)]


def assign_ids(samples, seed):
    """Return the samples in a shuffled order, with anonymous ids that say nothing of split, study or label."""
    order = np.random.default_rng(seed).permutation(len(samples))
    shuffled = samples.iloc[order].reset_index(drop=True)
    return shuffled.assign(**{ID_COLUMN: [f"sample_{number:04d}" for number in range(1, len(shuffled) + 1)]})


def write_split(folder, samples, values):
    """Write `input/data.csv` and `labels.csv` of one split."""
    (folder / INPUT_FILE).parent.mkdir(parents=True)
    data = values[samples.key].T.set_axis(samples[ID_COLUMN].values)
    data.rename_axis(ID_COLUMN).to_csv(folder / INPUT_FILE, float_format="%.0f", na_rep="NA")
    samples[[ID_COLUMN, LABEL_COLUMN]].to_csv(folder / LABELS_FILE, index=False)


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True, help="benchmark_bulk_counts.tsv.gz of the round")
    parser.add_argument("--training", type=Path, required=True, help="split/training_data_samples.tsv")
    parser.add_argument("--test", type=Path, required=True, help="split/test_samples.tsv")
    parser.add_argument("--annotation", type=Path, required=True, help="GEO's gene annotation (GeneID, Symbol)")
    parser.add_argument("--validation-share", type=float, required=True,
                        help="share of the laboratory groups of `development` that go to validation")
    parser.add_argument("--dataset-dir", type=Path, required=True, help="must not exist")
    parser.add_argument("--key-dir", type=Path, required=True, help="must not exist")
    arguments = parser.parse_args()
    for folder in (arguments.dataset_dir, arguments.key_dir):
        if folder.exists():
            raise SystemExit(f"error: {folder} exists")

    samples = draw_validation(list_samples(_read(arguments.training), _read(arguments.test)),
                              arguments.validation_share, SEED)
    validate_splits(samples)
    counts = pd.read_csv(arguments.counts, sep="\t", index_col=0)
    absent = sorted(set(samples.key) - set(counts.columns))
    if absent:
        raise ValueError(f"{len(absent)} samples have no column in the count table: {absent[:5]}")
    genes = choose_genes(counts, samples.key[samples.split != TEST])
    values = counts.loc[genes, samples.key]
    validate_counts(values)
    annotation = pd.read_csv(arguments.annotation, sep="\t", usecols=["GeneID", "Symbol"], index_col="GeneID").Symbol
    values.index = name_genes(genes, annotation)
    samples = assign_ids(samples, SEED)

    arguments.dataset_dir.mkdir(parents=True)
    for split in (TRAIN, VALIDATION, TEST):
        write_split(arguments.dataset_dir / split, samples[samples.split == split], values)
    metadata = {"task_type": "classification", "positive_class": LABELS["senescent"],
                "negative_class": LABELS["not senescent"]}
    (arguments.dataset_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (arguments.dataset_dir / "dataset_description.md").write_text(DESCRIPTION)
    arguments.key_dir.mkdir(parents=True)
    samples[KEY_COLUMNS].sort_values(ID_COLUMN).to_csv(arguments.key_dir / "id_key.tsv", sep="\t", index=False)

    print(f"genes: {len(genes)}; empty cells: {int(values.isna().sum().sum())} of {values.size}")
    print(samples.groupby("split").agg(samples=(ID_COLUMN, "size"), studies=("study", "nunique"),
                                       laboratory_groups=("laboratory_group", "nunique"),
                                       senescent=(LABEL_COLUMN, lambda labels: int((labels == "senescent").sum())))
          .to_string())
    for folder in (arguments.dataset_dir, arguments.key_dir):
        for path in sorted(p for p in folder.rglob("*") if p.is_file()):
            print(f"{_sha256(path)}  {path}")


if __name__ == "__main__":
    main()
