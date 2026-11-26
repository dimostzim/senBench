#!/usr/bin/env python3
"""Count the raw reads of GEO samples per NCBI gene id the way GEO counts its uniform tables.

Each run is aligned with HISAT2 (an SRA file is first written as FASTQ files) and counted with featureCounts; a run
at or below the minimum alignment rate is left out, and the counts of a sample's runs are summed. One count table per
dataset and featureCounts option set is written: NCBI gene ids by samples, rounded half up.

HISAT2 chooses among equally good alignments of a read differently with another number of threads, so the counts of a
sample depend on it; the column `threads` of the manifest gives the number for each sample.
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from download_manifest import dataset_folder
from list_raw_read_runs import PAIRED

GENE_COLUMN = "GeneID"
SAMPLES_FOLDER = "samples"
TABLES_FOLDER = "tables"
UNITS_NAME = "units.tsv"
RUNS_NAME = "runs.tsv"
RUN_COLUMNS = ("dataset", "sample", "run", "library_layout", "reads", "alignment_rate", "kept")
ALIGNMENT_RATE_PATTERN = re.compile(r"Overall alignment rate: ([\d.]+)%")
TOTAL_READS_PATTERN = re.compile(r"Total (?:reads|pairs): (\d+)")
FASTQ_SUFFIX = ".fastq.gz"


class _CountProblem(Exception):
    """A run that could not be aligned or counted, or counts that fail a check."""


def _run(command, log):
    """Run a command with its standard error appended to a log; raise with the log's end if it fails."""
    with log.open("a") as handle:
        result = subprocess.run(command, stderr=handle)
    if result.returncode:
        raise _CountProblem(f"{command[0]} exit {result.returncode}: {log.read_text()[-500:]!r}")


def convert_sra_file(sra_file, layout, threads, work, log):
    """Write the reads of an SRA file as FASTQ files in `work`, split as ENA splits them; return the files to align."""
    _run(["fasterq-dump", "--split-3", "--threads", str(threads), "--temp", str(work), "--outdir", str(work),
          str(sra_file)], log)
    run = sra_file.name
    return [work / f"{run}_1.fastq", work / f"{run}_2.fastq"] if layout == PAIRED else [work / f"{run}.fastq"]


def align_run(hisat2_program, files, layout, index, threads, bam, log):
    """Align one run with HISAT2 into an unsorted BAM; return the number of reads (pairs) and the alignment rate."""
    summary = bam.with_suffix(".summary.txt")
    reads = ["-1", str(files[0]), "-2", str(files[1])] if layout == PAIRED else ["-U", str(files[0])]
    with log.open("a") as handle:
        hisat2 = subprocess.Popen([str(hisat2_program), "-p", str(threads), "--new-summary", "--summary-file",
                                   str(summary), "-x", str(index), *reads], stdout=subprocess.PIPE, stderr=handle)
        # An uncompressed BAM costs no CPU time; it is removed after counting.
        samtools = subprocess.run(["samtools", "view", "-u", "-o", str(bam), "-"], stdin=hisat2.stdout,
                                  stderr=handle)
        hisat2.stdout.close()
        hisat2.wait()
    if hisat2.returncode or samtools.returncode:
        raise _CountProblem(f"hisat2 exit {hisat2.returncode}, samtools exit {samtools.returncode}: "
                            f"{log.read_text()[-500:]!r}")
    text = summary.read_text()
    return int(TOTAL_READS_PATTERN.search(text).group(1)), float(ALIGNMENT_RATE_PATTERN.search(text).group(1))


def count_run(bam, layout, features, options, threads, out, log, gene_lengths):
    """Count one aligned run with featureCounts; return the counts per gene and featureCounts' summary."""
    _run(["featureCounts", "-a", str(features), "-F", "SAF", "-T", str(threads), *options[layout.lower()],
          "-o", str(out), str(bam)], log)
    counts = pd.read_csv(out, sep="\t", comment="#", dtype={"Geneid": str}).set_index("Geneid")
    # The lengths that featureCounts gives its genes show that it counted the features GEO counted.
    if set(counts.index) != set(gene_lengths.index) or (counts.Length != gene_lengths[counts.index]).any():
        raise _CountProblem(f"{out}: the genes or their lengths differ from GEO's annotation table")
    summary = pd.read_csv(out.with_name(out.name + ".summary"), sep="\t", index_col=0).iloc[:, 0]
    return counts.iloc[:, -1], {status: int(value) for status, value in summary.items() if value}


def alignment_folder(work_dir, dataset, sample):
    """Return the folder of a sample's temporary alignments. HISAT2's wrapper splits its arguments at spaces, so the
    folder name keeps only letters, digits, `.`, `_` and `-` of the sample name."""
    return work_dir / dataset / dataset_folder(sample)


def validate_alignment_folders(work_dir, samples):
    """Refuse two samples (dataset, sample) whose temporary alignments would share a folder."""
    folders = [alignment_folder(work_dir, dataset, sample) for dataset, sample in samples]
    shared = [sample for sample, folder in zip(samples, folders) if folders.count(folder) > 1]
    if shared:
        raise ValueError(f"samples whose alignment folders would have the same name: {shared}")


def count_sample(dataset, sample, runs, reads_dir, hisat2_program, index, features, option_sets, gene_lengths,
                 threads, min_alignment_rate, out_dir, work_dir, keep_alignments):
    """Align and count the runs of one sample and write its counts; return (status, note, run rows)."""
    work = alignment_folder(work_dir, dataset, sample)
    work.mkdir(parents=True, exist_ok=True)
    log = work / "log.txt"
    run_rows, sums, record = [], {}, {"runs": {}}
    for run, rows in runs.groupby("run", sort=False):
        layout = rows.library_layout.iat[0]
        files = [reads_dir / dataset / file for file in rows.file]
        converted = not files[0].name.endswith(FASTQ_SUFFIX)
        if converted:
            files = convert_sra_file(files[0], layout, threads, work, log)
        bam = work / f"{run}.bam"
        reads, rate = align_run(hisat2_program, files, layout, index, threads, bam, log)
        if converted:
            for fastq in work.glob(f"{run}*.fastq"):
                fastq.unlink()
        kept = rate > min_alignment_rate
        run_rows.append((dataset, sample, run, layout, reads, rate, "yes" if kept else "no"))
        record["runs"][run] = {"library_layout": layout, "reads": reads, "alignment_rate": rate, "kept": kept,
                               "featurecounts": {}}
        if kept:
            for name, options in option_sets.items():
                counts, summary = count_run(bam, layout, features, options, threads, work / f"{run}.{name}.txt",
                                            log, gene_lengths)
                sums[name] = counts if name not in sums else sums[name] + counts
                record["runs"][run]["featurecounts"][name] = summary
        if not keep_alignments:
            bam.unlink()
    target = out_dir / SAMPLES_FOLDER / dataset
    target.mkdir(parents=True, exist_ok=True)
    (target / f"{sample}.json").write_text(json.dumps(record, indent=1))
    kept_runs = sum(row[-1] == "yes" for row in run_rows)
    if not kept_runs:
        return "skipped", f"no run above {min_alignment_rate}% aligned ({len(run_rows)} runs)", run_rows
    table = pd.DataFrame(sums).rename_axis(GENE_COLUMN)
    partial = target / f"{sample}.tsv.gz.partial"
    table.to_csv(partial, sep="\t", compression="gzip")
    partial.replace(target / f"{sample}.tsv.gz")
    if not keep_alignments:
        shutil.rmtree(work)
    first = next(iter(option_sets))
    return "ok", f"{kept_runs} of {len(run_rows)} runs; {table[first].sum():.0f} counted ({first})", run_rows


def write_tables(manifest, statuses, option_sets, out_dir):
    """Write per dataset and option set the table of the samples counted, in the order of the manifest."""
    for dataset, rows in manifest.groupby("accession", sort=False):
        samples = [sample for sample in dict.fromkeys(rows["sample"]) if statuses[(dataset, sample)] == "ok"]
        if not samples:
            continue
        counted = {sample: pd.read_csv(out_dir / SAMPLES_FOLDER / dataset / f"{sample}.tsv.gz", sep="\t",
                                       dtype={GENE_COLUMN: str}).set_index(GENE_COLUMN) for sample in samples}
        for name in option_sets:
            table = pd.DataFrame({sample: counts[name] for sample, counts in counted.items()})
            folder = out_dir / TABLES_FOLDER / name
            folder.mkdir(parents=True, exist_ok=True)
            np.floor(table + 0.5).astype("int64").to_csv(folder / f"{dataset}.tsv.gz", sep="\t", compression="gzip")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True, help="manifest of list_raw_read_runs.py")
    parser.add_argument("--reads-dir", type=Path, required=True, help="output folder of download_manifest.py")
    parser.add_argument("--hisat2", type=Path, required=True, help="the hisat2 program")
    parser.add_argument("--index", type=Path, required=True, help="prefix of the HISAT2 index")
    parser.add_argument("--features", type=Path, required=True, help="SAF file of make_geo_gene_features.py")
    parser.add_argument("--annotation", type=Path, required=True, help="GEO's human gene annotation table")
    parser.add_argument("--option-sets", type=Path, required=True,
                        help="JSON: per name the featureCounts options for paired and for single reads")
    parser.add_argument("--use", nargs="+", required=True, help="names of the option sets to count with")
    parser.add_argument("--datasets", nargs="+", help="count only the samples of these datasets of --runs")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--units-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True, help="folder of the temporary alignments")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--min-alignment-rate", type=float, default=50.0)
    parser.add_argument("--keep-alignments", action="store_true", help="keep the BAM files in --work-dir")
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} exists")

    manifest = pd.read_csv(args.runs, sep="\t", dtype=str, keep_default_na=False)
    if (manifest.groupby(["accession", "sample"]).threads.nunique() != 1).any():
        raise SystemExit(f"{args.runs}: the runs of a sample have different numbers of threads")
    if args.datasets:
        unknown = set(args.datasets) - set(manifest.accession)
        if unknown:
            raise SystemExit(f"{args.runs} has no runs of {sorted(unknown)}")
        manifest = manifest[manifest.accession.isin(args.datasets)]
    option_sets = {name: json.loads(args.option_sets.read_text())[name] for name in args.use}
    gene_lengths = pd.read_csv(args.annotation, sep="\t", usecols=["GeneID", "Length"],
                               dtype={"GeneID": str}).set_index("GeneID").Length
    samples = list(manifest.groupby(["accession", "sample"], sort=False))
    validate_alignment_folders(args.work_dir, [key for key, _ in samples])
    args.out_dir.mkdir(parents=True)

    def run_unit(item):
        (dataset, sample), runs = item
        try:
            status, note, run_rows = count_sample(dataset, sample, runs, args.reads_dir, args.hisat2, args.index,
                                                  args.features, option_sets, gene_lengths, int(runs.threads.iat[0]),
                                                  args.min_alignment_rate, args.out_dir, args.work_dir,
                                                  args.keep_alignments)
        except _CountProblem as problem:
            status, note, run_rows = "failed", str(problem), []
        return dataset, sample, status, note, run_rows

    statuses = {}
    units = args.units_dir / UNITS_NAME
    with ThreadPoolExecutor(max_workers=args.workers) as pool, units.open("w", newline="") as units_handle, \
            (args.out_dir / RUNS_NAME).open("w", newline="") as runs_handle:
        units_writer = csv.writer(units_handle, delimiter="\t", lineterminator="\n")
        runs_writer = csv.writer(runs_handle, delimiter="\t", lineterminator="\n")
        runs_writer.writerow(RUN_COLUMNS)
        for done, (dataset, sample, status, note, run_rows) in enumerate(pool.map(run_unit, samples), start=1):
            statuses[(dataset, sample)] = status
            units_writer.writerow((f"{dataset}/{sample}", status, note))
            units_handle.flush()
            runs_writer.writerows(run_rows)
            runs_handle.flush()
            print(f"PROGRESS {done}/{len(samples)} {dataset} {sample} {status}: {note}", flush=True)
    write_tables(manifest, statuses, option_sets, args.out_dir)
    failed = [unit for unit, status in statuses.items() if status == "failed"]
    if failed:
        sys.exit(f"failed units ({len(failed)}): {' '.join('/'.join(unit) for unit in failed)}")


if __name__ == "__main__":
    main()
