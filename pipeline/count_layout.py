#!/usr/bin/env python3
"""The layout in which an agent states where the supplementary files of a dataset hold raw read counts."""

RAW_READ_COUNTS = "raw read counts"
NORMALISED_VALUES_ONLY = "normalised values only"
NOTHING_USABLE = "nothing usable"
_READ_KINDS = ("text", "workbook")
_TAR = "tar"


class LayoutProblem(Exception):
    """A layout that does not fit the description of the files or the samples of the dataset."""


def output_layout():
    """Return the JSON Schema of one layout."""
    def number(minimum):
        return {"type": "integer", "minimum": minimum}

    def fixed_object(properties):
        return {"type": "object", "properties": properties, "required": list(properties),
                "additionalProperties": False}

    sample_column = fixed_object({"sample": {"type": "string"}, "member": {"type": "string"}, "column": number(1)})
    source = fixed_object({"file": {"type": "string"}, "sheet": {"type": "string"}, "column_names_row": number(0),
                           "first_values_row": number(1), "gene_identifier_column": number(1),
                           "sample_columns": {"type": "array", "items": sample_column, "minItems": 1}})
    contents = [RAW_READ_COUNTS, NORMALISED_VALUES_ONLY, NOTHING_USABLE]
    return fixed_object({"content": {"type": "string", "enum": contents},
                         "sources": {"type": "array", "items": source}, "note": {"type": "string"}})


def _find_sheets(files, source, member):
    """Return the sheet names of the file or member that a sample column names; none for a text."""
    name = source["file"]
    if name not in files:
        raise LayoutProblem(f"no file {name!r}")
    described = files[name]
    if (described["kind"] == _TAR) != bool(member):
        raise LayoutProblem(f"{name!r} is of kind {described['kind']} and the member is {member!r}")
    if member:
        members = {item["member"]: item for item in described["members"]}
        if member not in members:
            raise LayoutProblem(f"{name!r} has no member {member!r}")
        described = members[member]
    if described["kind"] not in _READ_KINDS:
        raise LayoutProblem(f"{member or name!r} is of kind {described['kind']}")
    return [sheet["sheet"] for sheet in described.get("sheets", [])]


def check_layout(layout, files, samples):
    """Raise unless the layout names files, members, sheets and samples that exist, each sample and column once."""
    if layout["sources"] and layout["content"] != RAW_READ_COUNTS:
        raise LayoutProblem(f"content is {layout['content']!r} with {len(layout['sources'])} sources")
    files = {described["file"]: described for described in files}
    seen_samples = set()
    seen_columns = set()
    for source in layout["sources"]:
        if source["first_values_row"] <= source["column_names_row"]:
            raise LayoutProblem(f"{source['file']!r}: the first row of values is not after the row of column names")
        for item in source["sample_columns"]:
            sheets = _find_sheets(files, source, item["member"])
            if (sheets or source["sheet"]) and source["sheet"] not in sheets:
                raise LayoutProblem(f"{item['member'] or source['file']!r} has no sheet {source['sheet']!r}")
            column = (source["file"], item["member"], source["sheet"], item["column"])
            if item["column"] == source["gene_identifier_column"] or column in seen_columns:
                raise LayoutProblem(f"column {item['column']} of {item['member'] or source['file']!r} is given twice")
            if item["sample"] not in samples or item["sample"] in seen_samples:
                raise LayoutProblem(f"{item['sample']!r} is not a sample of the dataset or is given twice")
            seen_columns.add(column)
            seen_samples.add(item["sample"])
