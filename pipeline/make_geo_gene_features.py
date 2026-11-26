#!/usr/bin/env python3
"""Write the gene features that GEO's uniform counts are counted on, as a featureCounts SAF file.

The features are the exons of each NCBI gene id in the GTF file of the annotation release, on the sequences of the
analysis set that the reads are aligned to; a gene without an exon but with a coding sequence is its whole span, and
a gene with neither is left out. Sequences are renamed from RefSeq accessions to the analysis set's UCSC-style names
with the release's assembly report.
"""

import argparse
import csv
import gzip
import re
from pathlib import Path

# Assembly units whose sequences the analysis set without alternate loci holds.
ANALYSIS_SET_UNITS = ("Primary Assembly", "non-nuclear")
SAF_COLUMNS = ("GeneID", "Chr", "Start", "End", "Strand")
GENE_ID_PATTERN = re.compile(r'db_xref "GeneID:(\d+)"')
GTF_GENE_ID_PATTERN = re.compile(r'gene_id "([^"]+)"')


def read_sequence_names(assembly_report, analysis_set_index):
    """Return the UCSC-style name of each RefSeq sequence of the analysis set's assembly units."""
    in_index = {line.split("\t", 1)[0] for line in analysis_set_index.read_text().splitlines()}
    names = {}
    with assembly_report.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            refseq, unit, ucsc = fields[6], fields[7], fields[9]
            if unit not in ANALYSIS_SET_UNITS:
                continue
            if ucsc not in in_index:
                raise SystemExit(f"{refseq} ({ucsc}) of the {unit} unit is not a sequence of {analysis_set_index}")
            names[refseq] = ucsc
    return names


def read_features(gtf, names):
    """Return the SAF rows: the exons of every gene on the named sequences, or the span of a gene without exons but
    with a coding sequence."""
    genes = {}
    exons = {}
    coding = set()
    with gzip.open(gtf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] not in names or fields[2] not in ("gene", "exon", "CDS"):
                continue
            # A GTF gene_id is unique on one sequence; exon rows need not carry the NCBI gene id.
            key = (fields[0], GTF_GENE_ID_PATTERN.search(fields[8]).group(1))
            feature = (names[fields[0]], int(fields[3]), int(fields[4]), fields[6])
            if fields[2] == "gene":
                if key in genes:
                    raise SystemExit(f"{gtf}: gene {key} occurs twice")
                genes[key] = (GENE_ID_PATTERN.search(fields[8]).group(1), feature)
            elif fields[2] == "exon":
                exons.setdefault(key, []).append(feature)
            else:
                coding.add(key)
    unknown = (set(exons) | coding) - set(genes)
    if unknown:
        raise SystemExit(f"{gtf}: {len(unknown)} genes of exon or CDS rows have no gene row, for example "
                         f"{sorted(unknown)[0]}")
    return [(gene_id, *feature) for key, (gene_id, span) in genes.items() if key in exons or key in coding
            for feature in exons.get(key, [span])]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True, help="GTF file of the annotation release")
    parser.add_argument("--assembly-report", type=Path, required=True, help="assembly report of the release")
    parser.add_argument("--analysis-set-index", type=Path, required=True, help=".fai of the analysis set")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists")

    rows = read_features(args.gtf, read_sequence_names(args.assembly_report, args.analysis_set_index))
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SAF_COLUMNS)
        writer.writerows(rows)
    print(f"{len(rows)} features of {len({row[0] for row in rows})} genes")


if __name__ == "__main__":
    main()
