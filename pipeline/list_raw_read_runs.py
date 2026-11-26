#!/usr/bin/env python3
"""List the raw read files of GEO samples from ENA, as a manifest for download_manifest.py.

A sample's runs are the runs of the SRA experiments that its GEO record names in `relations`; ENA's record of each run
must name the same experiment, with the GEO sample as the experiment's alias (or the sample followed by `_r` and a
number, as SRA names the experiments of newer GEO submissions). Where ENA has not made FASTQ files of a run yet, the
run's SRA file at NCBI is listed instead.
"""

import argparse
import csv
import gzip
import io
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
NCBI_LOCATOR = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
ENA_FIELDS = ("run_accession", "experiment_accession", "experiment_alias", "library_layout", "library_source",
              "library_strategy", "scientific_name", "read_count", "base_count", "fastq_ftp", "fastq_bytes",
              "fastq_md5")
CLIENT = "dataset-catalog"
# EBI and NCBI ask for at most 3 requests a second.
REQUEST_PAUSE_SECONDS = 0.4
EXPERIMENT_PATTERN = re.compile(r"sra\?term=([SED]RX\d+)")
TRANSCRIPTOMIC = "TRANSCRIPTOMIC"
PAIRED = "PAIRED"
MANIFEST_COLUMNS = ("accession", "sample", "experiment", "run", "library_layout", "library_strategy", "read_count",
                    "base_count", "file", "url", "bytes", "checksum_type", "checksum", "access", "note")
LEFT_OUT_COLUMNS = ("accession", "sample", "experiment", "run", "reason")

csv.field_size_limit(sys.maxsize)


def read_experiments(geo_samples, samples):
    """Return the SRA experiments that the GEO record of each sample names."""
    wanted = {sample for _, sample in samples}
    experiments = {}
    with gzip.open(geo_samples, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] in wanted:
                experiments[row["sample"]] = sorted(set(EXPERIMENT_PATTERN.findall(row["relations"])))
    missing = wanted - set(experiments)
    if missing:
        raise SystemExit(f"{geo_samples}: no record of {len(missing)} samples, for example {sorted(missing)[0]}")
    return experiments


def _get(url):
    result = subprocess.run(["curl", "-fsS", "-A", CLIENT, "--retry", "5", "--retry-delay", "10", url],
                            capture_output=True, text=True, check=True)
    time.sleep(REQUEST_PAUSE_SECONDS)
    return result.stdout


def fetch_runs(experiment):
    """Return ENA's rows of the runs of one experiment."""
    query = f"accession={experiment}&result=read_run&format=tsv&fields={','.join(ENA_FIELDS)}"
    return list(csv.DictReader(io.StringIO(_get(f"{ENA_FILEREPORT}?{query}")), delimiter="\t"))


def locate_sra_file(run):
    """Return the open address, bytes and md5 of the SRA file of a run at NCBI (not its file without qualities)."""
    files = [file for file in json.loads(_get(f"{NCBI_LOCATOR}?acc={run}"))["result"][0].get("files", [])
             if file["type"] == "sra" and file["name"] == run]
    links = [location["link"] for file in files for location in file["locations"]
             if "link" in location and not location.get("payRequired")]
    if not links:
        raise ValueError("ENA has no FASTQ file, and NCBI has no open SRA file")
    return links[0], str(files[0]["size"]), files[0]["md5"]


def select_files(run):
    """Return the (address, bytes, md5) of the read files to count: both files of a paired run, without a file of
    unpaired reads, or the one file of a single run; or the run's SRA file at NCBI."""
    if not run["fastq_ftp"]:
        return [locate_sra_file(run["run_accession"])]
    fields = (run[field].split(";") for field in ("fastq_ftp", "fastq_bytes", "fastq_md5"))
    files = [(f"https://{url}", size, md5) for url, size, md5 in zip(*fields)]
    if run["library_layout"] == PAIRED:
        files = [file for file in files if re.search(r"_[12]\.fastq\.gz$", file[0])]
        expected = 2
    else:
        expected = 1
    if len(files) != expected:
        raise ValueError(f"{len(files)} read files for a {run['library_layout']} run: {run['fastq_ftp']!r}")
    return files


def list_sample_runs(dataset, sample, experiments):
    """Return the manifest rows and the rows left out of one sample."""
    rows, left_out = [], []
    if not experiments:
        return rows, [(dataset, sample, "", "", "the GEO record names no SRA experiment")]
    for experiment in experiments:
        runs = fetch_runs(experiment)
        if not runs:
            left_out.append((dataset, sample, experiment, "", "ENA has no run of the experiment"))
        for run in runs:
            problem = None
            if run["experiment_accession"] != experiment or not re.fullmatch(rf"{sample}(_r\d+)?",
                                                                               run["experiment_alias"]):
                problem = (f"ENA's run names experiment {run['experiment_accession']} with alias "
                           f"{run['experiment_alias']}")
            elif run["library_source"] != TRANSCRIPTOMIC:
                problem = f"library source {run['library_source']}"
            # ENA leaves the organism of some runs empty.
            elif run["scientific_name"] not in ("Homo sapiens", ""):
                problem = f"organism {run['scientific_name']}"
            if problem is None:
                try:
                    files = select_files(run)
                except ValueError as error:
                    problem = str(error)
            if problem:
                left_out.append((dataset, sample, experiment, run["run_accession"], problem))
                continue
            for url, size, md5 in files:
                rows.append((dataset, sample, experiment, run["run_accession"], run["library_layout"],
                             run["library_strategy"], run["read_count"], run["base_count"],
                             f"{run['run_accession']}/{url.rsplit('/', 1)[1]}", url, size, "md5", md5, "open",
                             "" if run["fastq_ftp"] else "SRA file of NCBI; ENA has no FASTQ file"))
    return rows, left_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=Path, required=True, help="table with the columns `dataset` and `sample`")
    parser.add_argument("--geo-samples", type=Path, required=True, help="geo_samples.tsv.gz with `relations`")
    parser.add_argument("--out", type=Path, required=True, help="manifest of the read files")
    parser.add_argument("--left-out", type=Path, required=True, help="table of the samples and runs left out")
    args = parser.parse_args()
    for path in (args.out, args.left_out):
        if path.exists():
            raise SystemExit(f"{path} exists")

    with args.samples.open(newline="") as handle:
        samples = [(row["dataset"], row["sample"]) for row in csv.DictReader(handle, delimiter="\t")]
    experiments = read_experiments(args.geo_samples, samples)
    rows, left_out = [], []
    for done, (dataset, sample) in enumerate(samples, start=1):
        sample_rows, sample_left_out = list_sample_runs(dataset, sample, experiments[sample])
        rows += sample_rows
        left_out += sample_left_out
        print(f"PROGRESS {done}/{len(samples)} {dataset} {sample}: {len({row[3] for row in sample_rows})} runs, "
              f"{len(sample_left_out)} left out", flush=True)
    for path, columns, table in ((args.out, MANIFEST_COLUMNS, rows), (args.left_out, LEFT_OUT_COLUMNS, left_out)):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(table)
    gib = sum(int(row[10]) for row in rows) / 2**30
    print(f"{len(rows)} files of {len({row[3] for row in rows})} runs, {gib:.2f} GiB; {len(left_out)} left out")


if __name__ == "__main__":
    main()
