import csv

import pytest

from downloads import (ARRAYEXPRESS, GEO, MANIFESTS_FOLDER, RAW_READS, RECORD_NAME, REFERENCE, list_missing,
                       select_kinds, summarise_downloads, write_manifests)

COLUMNS = ("kind", "accession", "file", "bytes", "checksum_type", "checksum", "url", "access", "our_bytes", "our_md5")
ROWS = [dict(zip(COLUMNS, values)) for values in (
    (GEO, "GSE1", "soft/GSE1_family.soft.gz", "", "", "", "https://x/1", "open", "10", "a1"),
    (GEO, "GSE1", "suppl/GSE1_counts.txt.gz", "", "", "", "https://x/2", "open", "20", "a2"),
    (ARRAYEXPRESS, "E-MTAB-1", "E-MTAB-1.sdrf.txt", "5", "", "", "https://x/3", "open", "5", ""),
    (REFERENCE, "analysis_set", "index.tar.gz", "30", "md5", "a4", "https://x/4", "open", "30", "a4"),
    (RAW_READS, "GSE2", "SRR1/SRR1.fastq.gz", "40", "md5", "a5", "https://x/5", "open", "40", "a5"),
    (RAW_READS, "GSE3", "SRR2/SRR2.fastq.gz", "50", "md5", "a6", "https://x/6", "open", "50", "a6"),
)]


def place(data_dir, row, md5=None):
    """Put a file where download_manifest.py puts it and record it as that script does."""
    path = data_dir / row["kind"] / row["accession"] / row["file"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x")
    record = data_dir / row["kind"] / RECORD_NAME
    new = not record.exists()
    with record.open("a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if new:
            writer.writerow(("accession", "file", "bytes", "checksum_type", "checksum", "url", "status"))
        writer.writerow((row["accession"], row["file"], "1", "md5", md5 or row["our_md5"], row["url"], "ok"))


@pytest.mark.parametrize("without_raw_reads, kinds", [
    (False, ["geo", "geo_uniform_counts", "biostudies_arrayexpress", "reference", "raw_reads"]),
    (True, ["geo", "geo_uniform_counts", "biostudies_arrayexpress"]),
])
def test_without_raw_reads_leaves_out_the_reads_and_their_reference(without_raw_reads, kinds):
    assert select_kinds(without_raw_reads) == kinds


def test_manifests_hold_the_rows_of_their_kind_and_a_kind_left_out_loses_its_manifest(tmp_path):
    write_manifests(ROWS, tmp_path, select_kinds(False))
    by_kind = write_manifests(ROWS, tmp_path, select_kinds(True))
    assert sorted(path.name for path in (tmp_path / MANIFESTS_FOLDER).iterdir()) == [
        "biostudies_arrayexpress.tsv", "geo.tsv", "geo_uniform_counts.tsv"]
    with (tmp_path / MANIFESTS_FOLDER / "geo.tsv").open(newline="") as handle:
        assert [row["file"] for row in csv.DictReader(handle, delimiter="\t")] == [row["file"] for row in ROWS[:2]]
    assert len(by_kind["geo_uniform_counts"]) == 0


def test_summary_tells_missing_files_and_files_that_differ_from_ours(tmp_path):
    write_manifests(ROWS, tmp_path, select_kinds(True))
    place(tmp_path, ROWS[0])
    place(tmp_path, ROWS[1], md5="other")
    place(tmp_path, ROWS[4])
    summary = summarise_downloads(ROWS, tmp_path)
    assert [(row["file"], row["status"], row["same_as_ours"]) for row in summary] == [
        ("soft/GSE1_family.soft.gz", "ok", "yes"), ("suppl/GSE1_counts.txt.gz", "ok", "no"),
        ("E-MTAB-1.sdrf.txt", "missing", "unknown")]


def test_a_build_needs_the_reference_and_the_reads_only_of_the_datasets_it_counts(tmp_path):
    for row in ROWS[:3]:
        place(tmp_path, row)
    assert list_missing(ROWS, tmp_path, set()) == []
    assert [row["file"] for row in list_missing(ROWS, tmp_path, {"GSE3"})] == ["index.tar.gz", "SRR2/SRR2.fastq.gz"]
