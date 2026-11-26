#!/usr/bin/env python3
"""Read the supplementary files of a GEO dataset: delimited text, workbooks and the members of tar archives."""

import csv
import gzip
import io
import tarfile
import zipfile
from dataclasses import dataclass

import pandas as pd
import xlrd

TEXT, WORKBOOK, TAR, OTHER = "text", "workbook", "tar", "other"
# GEO's list of the members of a tar; not a data file.
MEMBER_LIST_NAME = "filelist.txt"
SEPARATORS = {"tab": "\t", "comma": ",", "semicolon": ";", "whitespace": r"\s+"}
_GZIP_START = b"\x1f\x8b"
_WORKBOOK_STARTS = (b"PK\x03\x04", b"\xd0\xcf\x11\xe0")
_TAR_MARK = slice(257, 262)
_SNIFF_BYTES = 8192
# The first non-empty lines of a text content, which give its delimiter and are shown in its description.
_FIRST_ROWS = 4
_MAX_LINE_CHARS = 50_000_000
_MAX_FIELD_CHARS = 300

csv.field_size_limit(_MAX_LINE_CHARS)


class ContentProblem(Exception):
    """A file or tar member that cannot be read as its kind."""


@dataclass(frozen=True)
class Table:
    """The values of some columns of one file, member or sheet, as text."""
    # Column number (counted in the rows of values, from 1) to the name of the column; empty without a names row.
    column_names: dict
    values: pd.DataFrame


def list_data_files(folder):
    """Return the data files of a supplementary folder, sorted by name; a missing folder has none."""
    if not folder.is_dir():
        return []
    return sorted(path for path in folder.iterdir() if path.is_file() and path.name != MEMBER_LIST_NAME)


def open_content(opener):
    """Return the kind of a file or tar member and a handle at the start of its content, uncompressed."""
    handle = opener()
    try:
        if handle.peek(len(_GZIP_START)).startswith(_GZIP_START):
            handle = io.BufferedReader(gzip.GzipFile(fileobj=handle))
        start = handle.peek(_SNIFF_BYTES)[:_SNIFF_BYTES]
    except (gzip.BadGzipFile, EOFError) as problem:
        raise ContentProblem(f"gzip: {problem}") from None
    if start.startswith(_WORKBOOK_STARTS):
        return WORKBOOK, handle
    if start[_TAR_MARK] == b"ustar":
        return TAR, handle
    if start and b"\x00" not in start:
        return TEXT, handle
    return OTHER, handle


def open_tar(handle):
    """Open the content of a tar for reading single members."""
    try:
        return tarfile.open(fileobj=handle, mode="r:")
    except tarfile.TarError as problem:
        raise ContentProblem(f"tar: {problem}") from None


def _read_lines(handle):
    """Yield (line number, line) for the lines of a text content that are not empty."""
    text = io.TextIOWrapper(handle, encoding="utf-8-sig", errors="replace")
    number = 0
    try:
        while line := text.readline(_MAX_LINE_CHARS):
            number += 1
            if len(line) == _MAX_LINE_CHARS:
                raise ContentProblem(f"line {number} has more than {_MAX_LINE_CHARS} characters")
            if line.strip():
                yield number, line.rstrip("\r\n")
    except (gzip.BadGzipFile, EOFError) as problem:
        raise ContentProblem(f"gzip: {problem}") from None


def _find_delimiter(lines):
    for name in ("tab", "comma", "semicolon"):
        if any(SEPARATORS[name] in line for line in lines):
            return name
    return "whitespace"


def _split_line(line, delimiter):
    if delimiter == "whitespace":
        return line.split()
    return next(csv.reader([line], delimiter=SEPARATORS[delimiter]))


def _read_text_rows(handle):
    """Return the delimiter of a text content and its first non-empty lines as (line number, fields)."""
    lines = [numbered for numbered, _ in zip(_read_lines(handle), range(_FIRST_ROWS))]
    delimiter = _find_delimiter([line for _, line in lines])
    return delimiter, [(number, _split_line(line, delimiter)) for number, line in lines]


def _read_sheets(handle, sheet, row_count):
    """Return the first rows of one sheet, or of every sheet by name when `sheet` is None, as text."""
    try:
        sheets = pd.read_excel(io.BytesIO(handle.read()), sheet_name=sheet, header=None, nrows=row_count, dtype=str,
                               na_filter=False)
    except (gzip.BadGzipFile, EOFError, zipfile.BadZipFile, xlrd.XLRDError, ValueError, KeyError) as problem:
        raise ContentProblem(f"workbook: {problem}") from None
    if sheet is None:
        return {name: _cells_as_text(frame) for name, frame in sheets.items()}
    return _cells_as_text(sheets)


def _cells_as_text(frame):
    # Some workbooks give a number or a missing value for a cell although text was asked for.
    return frame.map(lambda cell: "" if pd.isna(cell) else str(cell))


def _shorten(field):
    half = _MAX_FIELD_CHARS // 2
    return field if len(field) <= _MAX_FIELD_CHARS else f"{field[:half]} ... {field[-half:]}"


def _describe_rows(rows, max_fields):
    """Return rows as the description gives them: fields by column number, empty fields left out."""
    return [{"row": number, "field_count": len(fields),
             "fields": {str(column): _shorten(field) for column, field in enumerate(fields[:max_fields], start=1)
                        if field.strip()}}
            for number, fields in rows]


def _describe_content(kind, handle, max_fields):
    """Return the kind of a content and its first rows; a tar is not opened here."""
    if kind == TEXT:
        delimiter, rows = _read_text_rows(handle)
        return {"kind": kind, "delimiter": delimiter, "rows": _describe_rows(rows, max_fields)}
    if kind == WORKBOOK:
        sheets = _read_sheets(handle, None, _FIRST_ROWS)
        return {"kind": kind, "sheets": [
            {"sheet": name, "rows": _describe_rows([(number, list(row)) for number, row in
                                                    enumerate(frame.itertuples(index=False), start=1)], max_fields)}
            for name, frame in sheets.items()]}
    return {"kind": OTHER}


def _describe_members(tar, max_fields):
    """Return the members of a tar, each with its size, its kind and its first rows; a tar inside is not opened."""
    members = []
    try:
        for member in tar:
            if not member.isfile():
                continue
            described = {"member": member.name, "bytes": member.size}
            try:
                described.update(_describe_content(*open_content(lambda: tar.extractfile(member)), max_fields))
            except ContentProblem as problem:
                described.update({"kind": OTHER, "not_read": str(problem)})
            members.append(described)
    except tarfile.TarError as problem:
        raise ContentProblem(f"tar: {problem}") from None
    return members


def describe_file(opener, max_fields):
    """Return the kind of a file and its first rows, at most `max_fields` fields of a row; for a tar its members."""
    try:
        kind, handle = open_content(opener)
        if kind == TAR:
            return {"kind": kind, "members": _describe_members(open_tar(handle), max_fields)}
        return _describe_content(kind, handle, max_fields)
    except ContentProblem as problem:
        return {"kind": OTHER, "not_read": str(problem)}


def _name_columns(names, width, columns):
    """Return the names of columns of the rows of values; a names row may lack the name of the first column."""
    offset = width - len(names)
    if offset not in (0, 1):
        raise ContentProblem(f"the row of column names has {len(names)} fields, the first row of values {width}")
    return {column: names[column - 1 - offset].strip() if column > offset else "" for column in columns}


def _read_text_table(opener, column_names_row, first_values_row, columns):
    rows = {}
    delimiter = _read_text_rows(open_content(opener)[1])[0]
    for number, line in _read_lines(open_content(opener)[1]):
        if number in (column_names_row, first_values_row):
            rows[number] = _split_line(line, delimiter)
        if number >= first_values_row:
            break
    if first_values_row not in rows or (column_names_row and column_names_row not in rows):
        raise ContentProblem(f"no fields in row {column_names_row} or row {first_values_row}")
    names = _name_columns(rows[column_names_row], len(rows[first_values_row]), columns) if column_names_row else {}
    try:
        values = pd.read_csv(open_content(opener)[1], sep=SEPARATORS[delimiter], header=None,
                             skiprows=first_values_row - 1, usecols=[column - 1 for column in columns], dtype=str,
                             na_filter=False, skip_blank_lines=False, encoding="utf-8-sig",
                             encoding_errors="replace")
    except (gzip.BadGzipFile, EOFError, ValueError) as problem:
        raise ContentProblem(f"text: {problem}") from None
    return names, values.rename(columns=lambda position: position + 1)


def _read_sheet_table(opener, sheet, column_names_row, first_values_row, columns):
    frame = _read_sheets(open_content(opener)[1], sheet, None)
    if first_values_row > len(frame) or max(columns) > frame.shape[1]:
        raise ContentProblem(f"the sheet has {len(frame)} rows and {frame.shape[1]} columns")
    names = ({column: frame.iat[column_names_row - 1, column - 1].strip() for column in columns}
             if column_names_row else {})
    values = frame.iloc[first_values_row - 1:, [column - 1 for column in columns]]
    return names, values.set_axis(columns, axis=1)


def read_table(opener, sheet, column_names_row, first_values_row, columns):
    """Return the names and the values of the given columns of a text content or of a sheet of a workbook.

    Rows and columns are numbered from 1; `column_names_row` is 0 where there is no row of column names.
    Rows in which every given column is empty are left out.
    """
    kind = open_content(opener)[0]
    if kind == TEXT and not sheet:
        names, values = _read_text_table(opener, column_names_row, first_values_row, columns)
    elif kind == WORKBOOK and sheet:
        names, values = _read_sheet_table(opener, sheet, column_names_row, first_values_row, columns)
    else:
        raise ContentProblem(f"the content is of kind {kind} and the sheet is {sheet!r}")
    return Table(names, values[(values != "").any(axis=1)].reset_index(drop=True))
