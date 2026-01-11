"""
FASTA file parsing and input validation utilities.

This module provides helper functions to validate FASTA file paths
and parse FASTA-formatted files into sequence objects.
"""

# importing class sequence from stats.py
# importing Path from pathlib to validate the input path and raise errors
from bio_seq_v1.stats import sequence
from pathlib import Path


def validate_input_path(filename):
    """
    Validate an input file path.

    This function checks whether the provided filename exists, is a file
    (and not a directory), and is not empty.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the input file to be validated.

    Returns
    -------
    pathlib.Path
        A Path object representing the validated input file.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    IsADirectoryError
        If the path exists but is a directory.
    ValueError
        If the file exists but is empty.
    """
    path = Path(filename)  # setting the filename's path to a variable
    if not path.exists():  # if path does not exist
        raise FileNotFoundError(f"File not found: {path}")  # raise the filenotfound error
    if not path.is_file():  # if path is not a file (maybe the path provided is that of a directory)
        raise IsADirectoryError(f"Not a file: {path}")  # raise isadirectoryerror
    if path.stat().st_size == 0:  # if the file size is zero (empty file)
        raise ValueError(f"File is empty.")  # raise Valueerror
    return path


def fasta_parser(filename):
    """
    Parse a FASTA-formatted file and return sequence objects.

    This function reads a FASTA file, validates its format, and converts
    each record into a `sequence` object containing an identifier and
    its corresponding nucleotide or amino acid sequence.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the FASTA file to be parsed.

    Returns
    -------
    list
        A list of `sequence` objects parsed from the FASTA file.

    Raises
    ------
    ValueError
        If sequence data appears before a header line.
        If invalid characters are found in a sequence.
    FileNotFoundError
        If the input file does not exist.
    IsADirectoryError
        If the input path is a directory.
    """
    all_seq = []  # creating an empty list to store sequence objects (id + sequence)
    current_id = None  # setting current id to none
    current_seq = []  # creating an empty list for only sequences
    valid = set("ACGTUNRYSWKMBDHV-.")  # creating a set for valid bases
    path = validate_input_path(filename)  # validate the input file path

    with path.open('r') as fasta:  # opening the file and reading it
        for linenum, line in enumerate(fasta, start=1):  # running the loop on lines and numerating each line
            line = line.strip()  # removing whitespace if any

            if not line:  # if the line is empty (i.e no whitespaces), skip it
                continue

            if line.startswith(">"):  # header line indicating a new sequence
                if current_id is not None:
                    all_seq.append(sequence(current_id, "".join(current_seq)))
                current_id = line[1:]  # store the sequence ID (excluding '>')
                current_seq = []  # reset sequence accumulator
            else:
                if current_id is None:
                    raise ValueError(f"Sequence data found before header on line {linenum}")
                seq_line = line.upper()
                if all(base in valid for base in seq_line):
                    current_seq.append(seq_line)
                else:
                    raise ValueError(
                        f"Invalid character in sequence '{current_id}':'{seq_line}' "
                        f"on line: {linenum}"
                    )

        if current_id is not None:  # save the last sequence record
            all_seq.append(sequence(current_id, "".join(current_seq)))

    return all_seq  # return the list of sequence objects
