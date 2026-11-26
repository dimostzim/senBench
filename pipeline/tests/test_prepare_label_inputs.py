"""Check the label fields of an ArrayExpress sample, made from a made-up sample sheet and IDF."""

import pytest

from prepare_label_inputs import arrayexpress_sample, arrayexpress_series, read_idf

IDF = "\n".join([
    "Investigation Title\tRNA-seq of senescent fibroblasts",
    "Experimental Design\tstimulus or stress design\ttime series design",
    "Experiment Description\tFibroblasts were irradiated and became senescent.",
    "Protocol Name\tP-1\tP-2\tP-3\tP-4",
    "Protocol Type\tgrowth protocol\ttreatment protocol\tsample collection protocol\tnucleic acid extraction protocol",
    "Protocol Description\tDMEM with 10% serum\t10 Gy\tCells were harvested 10 days after irradiation\tRNeasy",
])
# One sample with two rows (one per read file), as an SDRF of paired reads has them.
ROW = [("Source Name", "IR_1"), ("Characteristics[organism]", "Homo sapiens"),
       ("Characteristics[sampling time point]", "10"), ("Unit[time unit]", "day"), ("Material Type", "cell"),
       ("Protocol REF", "P-1"), ("Protocol REF", "P-2"), ("Protocol REF", "P-3"), ("Protocol REF", "P-4"),
       ("Extract Name", "IR_1"), ("Material Type", "total RNA"), ("Comment[LIBRARY_SOURCE]", "TRANSCRIPTOMIC"),
       ("Factor Value[irradiate]", "10 Gy"), ("Factor Value[sampling time point]", "10"), ("Unit[time unit]", "day")]


@pytest.fixture
def idf(tmp_path):
    path = tmp_path / "E-MTAB-1.idf.txt"
    path.write_text(IDF + "\n")
    return read_idf(path)


def test_series_fields(idf):
    assert arrayexpress_series(idf) == {"title": "RNA-seq of senescent fibroblasts",
                                        "summary": "Fibroblasts were irradiated and became senescent.",
                                        "overall_design": "stimulus or stress design; time series design"}


def test_sample_fields(idf):
    second_row = ROW + [("Comment[FASTQ_URI]", "IR_1_2.fastq.gz")]
    sample = arrayexpress_sample("IR_1", [ROW, second_row], idf)
    assert sample["characteristics"] == [
        ["Characteristics[organism]", "Homo sapiens"], ["Characteristics[sampling time point]", "10"],
        ["Unit[time unit]", "day"], ["Factor Value[irradiate]", "10 Gy"], ["Factor Value[sampling time point]", "10"]]
    assert (sample["title"], sample["source_name"], sample["organism"]) == ("IR_1", "IR_1", "Homo sapiens")
    assert sample["molecule"] == "total RNA"
    assert (sample["growth_protocol"], sample["treatment_protocol"], sample["sample_collection_protocol"]) == (
        "DMEM with 10% serum", "10 Gy", "Cells were harvested 10 days after irradiation")
    assert (sample["library_source"], sample["library_selection"], sample["description"]) == ("TRANSCRIPTOMIC", "", "")


def test_protocol_missing_from_the_idf(idf):
    with pytest.raises(KeyError):
        arrayexpress_sample("IR_1", [ROW + [("Protocol REF", "P-9")]], idf)
