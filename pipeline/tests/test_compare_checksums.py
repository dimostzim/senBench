import gzip
import hashlib

from compare_checksums import DIFFERS, EQUAL, MISSING, NOT_BUILT, UNCOMPRESSED, compare_outputs


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def test_outputs_are_compared_by_file_or_by_uncompressed_content(tmp_path):
    (tmp_path / "round").mkdir()
    (tmp_path / "round" / "samples.tsv").write_bytes(b"a\tb\n")
    (tmp_path / "round" / "changed.tsv").write_bytes(b"new\n")
    # gzip stores the time of writing, so two writings of one table differ as files and not in their content.
    with gzip.GzipFile(tmp_path / "round" / "counts.tsv.gz", "wb", mtime=12345) as handle:
        handle.write(b"GeneID\tGSM1\n1\t5\n")
    checksums = [
        {"path": "round/samples.tsv", "hashed": "file", "sha256": sha256(b"a\tb\n")},
        {"path": "round/changed.tsv", "hashed": "file", "sha256": sha256(b"old\n")},
        {"path": "round/counts.tsv.gz", "hashed": UNCOMPRESSED, "sha256": sha256(b"GeneID\tGSM1\n1\t5\n")},
        {"path": "round/absent.tsv", "hashed": "file", "sha256": sha256(b"")},
        {"path": "reference/features.saf", "hashed": "file", "sha256": sha256(b"")},
    ]
    rows = compare_outputs(checksums, tmp_path, ["reference"])
    assert [row["outcome"] for row in rows] == [EQUAL, DIFFERS, EQUAL, MISSING, NOT_BUILT]
    assert rows[1]["sha256_of_this_build"] == sha256(b"new\n")
