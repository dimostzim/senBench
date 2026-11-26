#!/usr/bin/env python3
"""Check on made-up inputs which gene features make_geo_gene_features.py writes, which read files and runs
list_raw_read_runs.py takes, and how count_raw_reads.py sums and rounds the counts of a sample."""

import gzip
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import count_raw_reads as c  # noqa: E402
import list_raw_read_runs as l  # noqa: E402
import make_geo_gene_features as f  # noqa: E402

GTF_ROWS = [
    ("NC_1", "gene", 100, 500, "+", 'gene_id "A"; db_xref "GeneID:1";'),
    ("NC_1", "exon", 100, 200, "+", 'gene_id "A"; transcript_id "t1";'),
    ("NC_1", "exon", 400, 500, "+", 'gene_id "A"; transcript_id "t1";'),
    ("NC_1", "gene", 600, 700, "-", 'gene_id "B"; db_xref "GeneID:2";'),
    ("NC_1", "CDS", 600, 697, "-", 'gene_id "B"; transcript_id "t2";'),
    ("NC_1", "gene", 800, 900, "+", 'gene_id "C"; db_xref "GeneID:3";'),
    ("NC_2", "gene", 100, 300, "+", 'gene_id "A"; db_xref "GeneID:1";'),
    ("NC_2", "exon", 100, 300, "+", 'gene_id "A"; transcript_id "t3";'),
    ("NW_9", "gene", 100, 300, "+", 'gene_id "D"; db_xref "GeneID:4";'),
    ("NW_9", "exon", 100, 300, "+", 'gene_id "D"; transcript_id "t4";'),
]


def read_features():
    """Return the SAF rows of a made-up GTF: NC_1 and NC_2 are sequences of the analysis set, NW_9 is not."""
    with tempfile.TemporaryDirectory() as folder:
        gtf = Path(folder) / "genes.gtf.gz"
        with gzip.open(gtf, "wt") as handle:
            for name, kind, start, end, strand, attributes in GTF_ROWS:
                handle.write(f"{name}\tRefSeq\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
        return f.read_features(gtf, {"NC_1": "chr1", "NC_2": "chrY"})


def select_files(layout, names):
    run = {"run_accession": "SRR1", "library_layout": layout, "fastq_ftp": ";".join(names),
           "fastq_bytes": ";".join("1" for _ in names), "fastq_md5": ";".join("m" for _ in names)}
    try:
        return [url for url, _, _ in l.select_files(run)]
    except ValueError as problem:
        return str(problem).split(":")[0]


def list_runs(alias, source="TRANSCRIPTOMIC", organism="Homo sapiens"):
    """Return the runs taken and the reasons of the runs left out, for a sample whose experiment has one run."""
    run = {"run_accession": "SRR1", "experiment_accession": "SRX1", "experiment_alias": alias,
           "library_layout": "SINGLE", "library_source": source, "library_strategy": "RNA-Seq",
           "scientific_name": organism, "read_count": "10", "base_count": "500", "fastq_ftp": "host/SRR1.fastq.gz",
           "fastq_bytes": "1", "fastq_md5": "m"}
    l.fetch_runs = lambda experiment: [run]
    rows, left_out = l.list_sample_runs("GSE1", "GSM12", ["SRX1"])
    return [row[3] for row in rows], [row[-1] for row in left_out]


def validate_alignment_folders(samples):
    """Return the error of count_raw_reads.py for these samples, or `accepted`."""
    try:
        c.validate_alignment_folders(Path("w"), samples)
    except ValueError as problem:
        return str(problem).split(":")[0]
    return "accepted"


def write_tables(values):
    """Return the table that count_raw_reads.py writes for one sample with these summed counts of three genes."""
    with tempfile.TemporaryDirectory() as folder:
        out = Path(folder)
        (out / c.SAMPLES_FOLDER / "GSE1").mkdir(parents=True)
        pd.DataFrame({c.GENE_COLUMN: ["1", "2", "3"], "set": values}).to_csv(
            out / c.SAMPLES_FOLDER / "GSE1" / "GSM1.tsv.gz", sep="\t", index=False)
        manifest = pd.DataFrame({"accession": ["GSE1", "GSE1"], "sample": ["GSM1", "GSM2"]})
        c.write_tables(manifest, {("GSE1", "GSM1"): "ok", ("GSE1", "GSM2"): "skipped"}, {"set": {}}, out)
        table = pd.read_csv(out / c.TABLES_FOLDER / "set" / "GSE1.tsv.gz", sep="\t", dtype={c.GENE_COLUMN: str})
    return table.to_dict("list")


CASES = {
    "exons of a gene on every sequence, the span of a gene with a coding sequence only, no gene without either": (
        read_features(),
        [("1", "chr1", 100, 200, "+"), ("1", "chr1", 400, 500, "+"), ("2", "chr1", 600, 700, "-"),
         ("1", "chrY", 100, 300, "+")]),
    "both files of a paired run, without the file of unpaired reads": (
        select_files("PAIRED", ["h/SRR1.fastq.gz", "h/SRR1_1.fastq.gz", "h/SRR1_2.fastq.gz"]),
        ["https://h/SRR1_1.fastq.gz", "https://h/SRR1_2.fastq.gz"]),
    "the one file of a single run": (select_files("SINGLE", ["h/SRR1.fastq.gz"]), ["https://h/SRR1.fastq.gz"]),
    "a paired run with one file fails": (select_files("PAIRED", ["h/SRR1.fastq.gz"]),
                                         "0 read files for a PAIRED run"),
    "an experiment named by the sample is taken": (list_runs("GSM12"), (["SRR1"], [])),
    "an experiment named by the sample and _r1 is taken": (list_runs("GSM12_r1"), (["SRR1"], [])),
    "an experiment named by another sample is left out": (
        list_runs("GSM123"), ([], ["ENA's run names experiment SRX1 with alias GSM123"])),
    "a run of another library source is left out": (list_runs("GSM12", source="GENOMIC"),
                                                    ([], ["library source GENOMIC"])),
    "a run without an organism at ENA is taken, one of another organism is left out": (
        (list_runs("GSM12", organism=""), list_runs("GSM12", organism="Mus musculus")),
        ((["SRR1"], []), ([], ["organism Mus musculus"]))),
    "summed counts are rounded half up, and a skipped sample has no column": (
        write_tables([0.5, 2.5, 2.33]), {c.GENE_COLUMN: ["1", "2", "3"], "GSM1": [1, 3, 2]}),
    "the alignment folder of a sample name with a space has none; a GEO sample keeps its name": (
        (c.alignment_folder(Path("w"), "E-MTAB-1", "Sample 1"), c.alignment_folder(Path("w"), "GSE1", "GSM1")),
        (Path("w/E-MTAB-1/Sample_1"), Path("w/GSE1/GSM1"))),
    "two samples of a dataset whose alignment folders would have the same name are refused, others are accepted": (
        (validate_alignment_folders([("E-MTAB-1", "Sample 1"), ("E-MTAB-1", "Sample_1")]),
         validate_alignment_folders([("E-MTAB-1", "Sample 1"), ("E-MTAB-2", "Sample 1"), ("E-MTAB-1", "Sample 2")])),
        ("samples whose alignment folders would have the same name", "accepted")),
}

failed = 0
for name, (got, expected) in CASES.items():
    ok = got == expected
    failed += not ok
    print(("ok   " if ok else "FAIL ") + name + ("" if ok else f"\n      got      {got}\n      expected {expected}"))
sys.exit(1 if failed else 0)
