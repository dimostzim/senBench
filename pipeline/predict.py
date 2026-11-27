#!/usr/bin/env python3
"""Score the labelled samples of `test` of a build with Agentomics or with one of the six tools, and write the method's
predictions.tsv (columns dataset, sample and score, one row per scored sample), replacing the shipped file.

Agentomics: agentomics_model/inference.py scores `test/input` of the build's Agentomics dataset, and the dataset's key
maps each anonymous id to its dataset and sample. Its environment is created from agentomics_model/environment.yml
when it does not exist yet.

A tool: its wrapper, set up by tools/<tool>/setup.sh, is called once per sample on a count table that holds only that
sample, cut from the round's count table of its dataset, so that every score is computed from one sample alone. hUSI's
wrapper refuses a single sample, so it is called once per dataset on the round's count table; the score it reports
(`hUSI_unscaled`) depends on the sample alone. A sample on which the wrapper stops gets no row.

The work files go to <build>/predict/<method>/, which must not exist. Needs only Python 3.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ROUND = "all980_v10"
PREDICTIONS_NAME = "predictions.tsv"
AGENTOMICS = "agentomics"
# Every method of the results table, by its folder name, with its name in the results.
METHODS = {AGENTOMICS: "Agentomics", "husi": "hUSI", "sencid": "SenCID", "sencat": "SenCat", "sencan": "SENCAN",
           "marker": "markeR", "senfinder": "SenFinder"}
# The score of each tool that the benchmark reports, as the score_name of the wrapper's scores.tsv.
TOOL_SCORES = {"husi": "hUSI_unscaled", "sencid": "RecSID_Score", "sencat": "rna_score", "sencan": "SENCAN_score",
               "marker": "HernandezSegura:ssGSEA", "senfinder": "Senescence_probability"}
PER_DATASET = "husi"
ANNOTATION = Path("geo_uniform_counts") / "annotation" / "Human.GRCh38.p13.annot.tsv.gz"
# The numeric libraries of a wrapper use one thread, as in the benchmark's runs, which ran many wrappers at a time.
ONE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def prediction_file(method):
    """Return the predictions.tsv of a method: in agentomics_model/ or in the tool's folder."""
    return ROOT / ("agentomics_model" if method == AGENTOMICS else f"tools/{method}") / PREDICTIONS_NAME


def read_rows(path, opener=open):
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_predictions(method, scores):
    """Write one row per (dataset, sample) of `scores`, in their order, each score as Python's shortest decimal."""
    with prediction_file(method).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("dataset", "sample", "score"))
        writer.writerows((dataset, sample, repr(float(score))) for (dataset, sample), score in sorted(scores.items()))
    print(f"{len(scores)} samples scored; {prediction_file(method)}")


def create_environment(prefix, conda):
    """Create the environment of the Agentomics model once; a folder left by a creation that failed is created again."""
    complete = prefix / ".senbench_complete"
    if complete.exists():
        return
    shutil.rmtree(prefix, ignore_errors=True)
    environment = {**os.environ, "CONDA_ALWAYS_YES": "true", "MAMBA_ALWAYS_YES": "true"}
    subprocess.run([conda, "env", "create", "--prefix", str(prefix), "--file",
                    str(ROOT / "agentomics_model" / "environment.yml")], check=True, env=environment)
    complete.touch()


def predict_agentomics(build_dir, conda, work_dir):
    """Return the probability of senescent of every `test` sample, by dataset and sample."""
    prefix = build_dir / "environments" / AGENTOMICS
    create_environment(prefix, conda)
    dataset = build_dir / "agentomics" / f"senescence_bulk_{ROUND}"
    output = work_dir / "inference.csv"
    subprocess.run([str(prefix / "bin" / "python"), str(ROOT / "agentomics_model" / "inference.py"), "--input",
                    str(dataset / "test" / "input"), "--output", str(output)], check=True)
    key = {row["id"]: row for row in read_rows(dataset.with_name(dataset.name + "_key") / "id_key.tsv")}
    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if {row["id"] for row in rows} != {id_ for id_, row in key.items() if row["split"] == "test"}:
        raise SystemExit(f"{output}: the ids differ from those of `test` in the dataset's key")
    return {(key[row["id"]]["dataset"], key[row["id"]]["sample"]): row["probability_1"] for row in rows}


def tool_command(tool, tools_dir):
    """Return the call of a tool's wrapper, without its input, and the variables it needs."""
    wrapper = ["bash", str(ROOT / "tools" / tool / "run_tool.sh")]
    envs = tools_dir / "envs"
    commands = {
        "husi": wrapper + ["--env-prefix", str(envs / "husi"), "--tool-repository", str(tools_dir / "HUSI")],
        # SenCID's wrapper reads its paths from its config.sh, so it runs from the copy that setup.sh writes.
        "sencid": ["bash", str(tools_dir / "wrappers" / "sencid" / "run_tool.sh")],
        "sencat": wrapper,
        "sencan": wrapper + ["--env-prefix", str(envs / "sencan"), "--tool-repository", str(tools_dir / "HUSI")],
        "marker": wrapper + ["--env-prefix", str(envs / "marker")],
        "senfinder": wrapper + ["--set-absent-genes-to-zero", "--env-prefix", str(envs / "senfinder"),
                                "--tool-repository", str(tools_dir / "SenFinder")],
    }
    variables = {**ONE_THREAD, **({"SENESCENCE_TOOLS_DIR": str(tools_dir)} if tool == "sencat" else {})}
    return commands[tool], variables


def write_one_sample_table(counts, sample, target):
    """Write the column of one sample of a count table, with its gene column, as a table of its own."""
    with gzip.open(counts, "rt", newline="") as source, gzip.open(target, "wt", newline="") as handle:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader)
        column = header.index(sample)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((header[0], sample))
        writer.writerows((row[0], row[column]) for row in reader)


def list_units(tool, test_samples, counts_dir, inputs_dir):
    """Return the units of a tool: (name, count table, the samples it scores as (dataset, sample))."""
    by_dataset = {}
    for row in test_samples:
        by_dataset.setdefault(row["dataset"], []).append(row["sample"])
    if tool == PER_DATASET:
        return [(dataset, counts_dir / f"{dataset}.tsv.gz", [(dataset, sample) for sample in samples])
                for dataset, samples in by_dataset.items()]
    units = []
    for dataset, samples in by_dataset.items():
        for position, sample in enumerate(samples, start=1):
            name = f"{dataset}__{position}"
            table = inputs_dir / f"{name}.tsv.gz"
            write_one_sample_table(counts_dir / f"{dataset}.tsv.gz", sample, table)
            units.append((name, table, [(dataset, sample)]))
    return units


def predict_tool(tool, build_dir, downloads_dir, tools_dir, workers, work_dir):
    """Return the reported score of every `test` sample on which the tool's wrapper did not stop."""
    command, variables = tool_command(tool, tools_dir)
    annotation = downloads_dir / ANNOTATION
    round_dir = build_dir / ROUND
    test_samples = read_rows(round_dir / "split" / "test_samples.tsv")
    (work_dir / "inputs").mkdir(parents=True)
    units = list_units(tool, test_samples, round_dir / "counts", work_dir / "inputs")
    environment = {**os.environ, **variables}

    def run_unit(unit):
        name, table, samples = unit
        out = work_dir / "runs" / name
        out.parent.mkdir(parents=True, exist_ok=True)
        with (work_dir / "logs" / f"{name}.log").open("w") as log:
            result = subprocess.run(command + ["--counts", str(table), "--annotation", str(annotation), "--out-dir",
                                               str(out)], stdout=log, stderr=subprocess.STDOUT, env=environment)
        if result.returncode != 0:
            return name, "failed", f"exit code {result.returncode}", {}
        names = {sample: (dataset, sample) for dataset, sample in samples}
        scores = {names[row["sample"]]: row["score"] for row in read_rows(out / "scores.tsv")
                  if row["score_name"] == TOOL_SCORES[tool] and row["sample"] in names}
        if len(scores) != len(samples):
            return name, "failed", f"{len(scores)} of {len(samples)} samples have a {TOOL_SCORES[tool]}", {}
        return name, "ok", "", scores

    (work_dir / "logs").mkdir()
    scores = {}
    with ThreadPoolExecutor(max_workers=workers) as pool, (work_dir / "units.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for done, (name, status, note, unit_scores) in enumerate(pool.map(run_unit, units), start=1):
            writer.writerow((name, status, note))
            print(f"PROGRESS {done}/{len(units)} {name} {status}{': ' + note if note else ''}", flush=True)
            scores.update(unit_scores)
    return scores


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("method", choices=list(METHODS))
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build", help="the folder of build.sh")
    parser.add_argument("--downloads-dir", type=Path, default=ROOT / "downloads",
                        help="the folder of download.sh, for GEO's annotation table")
    parser.add_argument("--tools-dir", type=Path, default=ROOT / "build" / "tools",
                        help="the install folder of the tools' setup.sh")
    parser.add_argument("--conda", default="conda", help="conda, mamba or micromamba, to create the environment of "
                                                         "the Agentomics model")
    parser.add_argument("--workers", type=int, default=8, help="wrapper calls at a time")
    args = parser.parse_args()
    work_dir = args.build_dir / "predict" / args.method
    if work_dir.exists():
        raise SystemExit(f"{work_dir} exists; delete it to score again")
    work_dir.mkdir(parents=True)
    if args.method == AGENTOMICS:
        scores = predict_agentomics(args.build_dir.resolve(), args.conda, work_dir)
    else:
        scores = predict_tool(args.method, args.build_dir.resolve(), args.downloads_dir.resolve(),
                              args.tools_dir.resolve(), args.workers, work_dir)
    write_predictions(args.method, scores)


if __name__ == "__main__":
    main()
