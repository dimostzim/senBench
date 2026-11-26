#!/usr/bin/env python3
"""The layout in which an agent states where the supplementary files of a dataset hold values per gene and sample of
any kind (raw read counts, expected counts, normalised values or log values), and how each column is identified."""

from count_layout import RAW_READ_COUNTS, LayoutProblem, check_layout

VALUES_PER_SAMPLE = "values per sample"
NOTHING_PER_SAMPLE = "nothing per sample"
EXPECTED_COUNTS = "expected counts"
NORMALISED_VALUES = "normalised values"
LOG_VALUES = "log values"
BASIS_NAMES = "names"
BASIS_ORDER_ONLY = "order only"


def output_values_layout():
    """Return the JSON Schema of one layout."""
    def number(minimum):
        return {"type": "integer", "minimum": minimum}

    def fixed_object(properties):
        return {"type": "object", "properties": properties, "required": list(properties),
                "additionalProperties": False}

    sample_column = fixed_object({"sample": {"type": "string"}, "member": {"type": "string"}, "column": number(1),
                                  "basis": {"type": "string", "enum": [BASIS_NAMES, BASIS_ORDER_ONLY]}})
    values = [RAW_READ_COUNTS, EXPECTED_COUNTS, NORMALISED_VALUES, LOG_VALUES]
    source = fixed_object({"file": {"type": "string"}, "sheet": {"type": "string"}, "column_names_row": number(0),
                           "first_values_row": number(1), "gene_identifier_column": number(1),
                           "values": {"type": "string", "enum": values},
                           "sample_columns": {"type": "array", "items": sample_column, "minItems": 1}})
    return fixed_object({"content": {"type": "string", "enum": [VALUES_PER_SAMPLE, NOTHING_PER_SAMPLE]},
                         "sources": {"type": "array", "items": source}, "note": {"type": "string"}})


def check_values_layout(layout, files, samples):
    """Raise unless the layout names files, members, sheets and samples that exist, each sample and column once."""
    if bool(layout["sources"]) != (layout["content"] == VALUES_PER_SAMPLE):
        raise LayoutProblem(f"content is {layout['content']!r} with {len(layout['sources'])} sources")
    # The checks of files, members, sheets, samples and columns are those of a layout of raw read counts.
    check_layout({**layout, "content": RAW_READ_COUNTS}, files, samples)
