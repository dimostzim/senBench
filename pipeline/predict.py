#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"
TEST_SAMPLES = BUILD / "all980_v10" / "split" / "test_samples.tsv"
COUNTS = BUILD / "all980_v10" / "counts"
AGENTOMICS_DATASET = BUILD / "agentomics" / "senescence_bulk_all980_v10"
AGENTOMICS_KEY = BUILD / "agentomics" / "senescence_bulk_all980_v10_key" / "id_key.tsv"
ANNOTATION = ROOT / "downloads" / "geo_uniform_counts" / "annotation" / "Human.GRCh38.p13.annot.tsv.gz"
AGENTOMICS = "Agentomics"
TOOL_SCORES = {"hUSI": "hUSI_unscaled", "SenCID": "RecSID_Score", "SenCat": "rna_score", "SENCAN": "SENCAN_score",
               "markeR": "HernandezSegura:ssGSEA", "SenFinder": "Senescence_probability"}
SCORED_PER_DATASET = "hUSI"
WORKERS = len(os.sched_getaffinity(0))


def read_rows(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def predict_agentomics(work_dir):
    output = work_dir / "inference.csv"
    subprocess.run([str(BUILD / "envs" / AGENTOMICS / "bin" / "python"), str(ROOT / "agentomics_model" / "inference.py"),
                    "--input", str(AGENTOMICS_DATASET / "test" / "input"), "--output", str(output)], check=True)
    key = {row["id"]: row for row in read_rows(AGENTOMICS_KEY)}
    rows = read_rows(output, ",")
    if {row["id"] for row in rows} != {id_ for id_, row in key.items() if row["split"] == "test"}:
        raise SystemExit(f"{output}: the ids differ from those of `test` in the dataset's key")
    return {(key[row["id"]]["dataset"], key[row["id"]]["sample"]): row["probability_1"] for row in rows}


def write_one_sample_table(counts, sample, target):
    with gzip.open(counts, "rt", newline="") as source, gzip.open(target, "wt", newline="") as handle:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        column = header.index(sample)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((header[0], sample))
        writer.writerows((row[0], row[column]) for row in reader)


def list_units(tool, inputs_dir):
    by_dataset = {}
    for row in read_rows(TEST_SAMPLES):
        by_dataset.setdefault(row["dataset"], []).append(row["sample"])
    if tool == SCORED_PER_DATASET:
        return [(dataset, COUNTS / f"{dataset}.tsv.gz", samples) for dataset, samples in by_dataset.items()]
    units = []
    for dataset, samples in by_dataset.items():
        for position, sample in enumerate(samples, start=1):
            table = inputs_dir / f"{dataset}__{position}.tsv.gz"
            write_one_sample_table(COUNTS / f"{dataset}.tsv.gz", sample, table)
            units.append((dataset, table, [sample]))
    return units


def predict_tool(tool, work_dir):
    (work_dir / "inputs").mkdir()
    (work_dir / "logs").mkdir()
    units = list_units(tool, work_dir / "inputs")

    def run_unit(unit):
        dataset, table, samples = unit
        name = table.name.removesuffix(".tsv.gz")
        out = work_dir / "runs" / name
        out.parent.mkdir(exist_ok=True)
        with (work_dir / "logs" / f"{name}.log").open("w") as log:
            result = subprocess.run(["bash", str(ROOT / "tools" / tool / "run_tool.sh"), str(table), str(ANNOTATION),
                                     str(out)], stdout=log, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            return name, f"failed: exit code {result.returncode}", {}
        scores = {(dataset, row["sample"]): row["score"] for row in read_rows(out / "scores.tsv")
                  if row["score_name"] == TOOL_SCORES[tool] and row["sample"] in samples}
        if len(scores) != len(samples):
            return name, f"failed: {len(scores)} of {len(samples)} samples have a {TOOL_SCORES[tool]}", {}
        return name, "ok", scores

    scores = {}
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for done, (name, status, unit_scores) in enumerate(pool.map(run_unit, units), start=1):
            print(f"PROGRESS {done}/{len(units)} {name} {status}", flush=True)
            scores.update(unit_scores)
    return scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("method", choices=[AGENTOMICS, *TOOL_SCORES])
    method = parser.parse_args().method
    work_dir = BUILD / "predict" / method
    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True)
    scores = predict_agentomics(work_dir) if method == AGENTOMICS else predict_tool(method, work_dir)
    predictions = ROOT / "predictions" / f"{method}.tsv"
    with predictions.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("dataset", "sample", "score"))
        writer.writerows((dataset, sample, repr(float(score))) for (dataset, sample), score in sorted(scores.items()))
    print(f"{len(scores)} samples scored; {predictions}")


if __name__ == "__main__":
    main()
