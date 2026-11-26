"""Repository accession patterns shared by the inventory scripts, and the key that names a sample across datasets."""

import re

ACCESSION = re.compile(
    r"\b("
    r"GSE\d+|E-[A-Z]{4}-\d+|PRJ(?:NA|EB|DB|CA)\d+|[SED]RP\d+|"
    r"PXD\d+|MSV\d+|JPST\d+|IPX\d+|HRA\d+|CRA\d+|OMIX\d+|"
    r"S-BSST\d+|MTBLS\d+|ST\d{6}|EGAS\d+|phs\d{6}"
    r")\b"
)
GEO_SAMPLE = re.compile(r"GSM\d+")


def sample_key(dataset, sample):
    """Return the key that names one sample across datasets.

    A GEO sample accession names the same sample in every series that holds it, such as a SuperSeries and its
    SubSeries. Any other sample name, such as the source name of an ArrayExpress sample sheet, names a sample of its own
    dataset only, so its key is the dataset and the sample.
    """
    return sample if GEO_SAMPLE.fullmatch(sample) else f"{dataset}/{sample}"


def sample_keys(datasets, samples):
    """Return sample_key() of each pair of a dataset and a sample."""
    return [sample_key(dataset, sample) for dataset, sample in zip(datasets, samples, strict=True)]


def source_keys(text):
    """Return the accessions named in a record name, or its leading token when it names none."""
    found = sorted(set(ACCESSION.findall(text)))
    return found or [re.split(r"[/:]", text.strip(), maxsplit=1)[0]]
