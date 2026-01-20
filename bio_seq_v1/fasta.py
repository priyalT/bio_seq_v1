"""
FASTA file parsing and input validation utilities.

This module provides helper functions to validate FASTA file paths
and parse FASTA-formatted files into sequence objects.
"""

from bio_seq_v1.stats import sequence
from pathlib import Path
from typing import Optional


class FASTAParser:
    def __init__(self, path: Optional[str] = None, strict: bool = False, strict_file: bool = False, strict_seq: bool = False):
        self.path = Path(path) if path else None
        self.sequences = []
        self.errors = []
        self.warnings = []
        self.strict = strict
        self.strict_file = strict_file
        self.strict_seq = strict_seq

        if self.strict_file:
            self._strict_file_validate()

    @classmethod
    def strict_mode(cls, path: str):
        return cls(path, strict=True)

    def _strict_file_validate(self):
        if not self.path:
            raise ValueError("Strict mode requires a file path")
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} does not exist")
        if not self.path.is_file():
            raise IsADirectoryError(f"{self.path} is a directory")
        if self.path.stat().st_size == 0:
            raise ValueError(f"File is empty.")
        
    def _validate_sequence(self, line, linenum):
        valid_chars = {"A","C","G","T","N","R","Y","S","W","K","M","B","D","H","V","-","."}
        invalid_chars = set(line) - valid_chars
        if invalid_chars:
            msg = f"Invalid character(s) in sequence at line {linenum}: {''.join(invalid_chars)}"
            if self.strict_seq:
                raise ValueError(msg)
            else:
                self.errors.append(msg)
        
              
    def _parse_lines(self, lines):
            seq = []
            header = None
            for linenum, line in enumerate(lines, start =1):
                if not line:
                    continue
                line = line.strip()
                if not line():
                    self.errors.append(f"Empty or whitespace-only sequence at line {linenum}")
                    continue
                if line.startswith(">"):
                    if header is None and seq:
                        msg = f"Sequence found before header at line {linenum}: {seq}"
                        if self.strict:
                            raise ValueError(msg)
                        else:
                            self.errors.append(msg)
                    if header:
                        self.sequences.append(sequence(header, "".join(seq)))
                    header = line[1:]
                    seq = []
                    if header and not seq:
                        msg = f"Header '{header}' has no sequence"
                        if self.strict:
                            raise ValueError(msg)
                        else:
                            self.errors.append(msg)
                else:
                    if not header:
                        msg = f"Sequence line before any header at line {linenum}"
                        if self.strict:
                            raise ValueError(msg)
                        else:
                            self.errors.append(msg)
                    self._validate_sequence(line, linenum)
                    seq.append(line.upper())
                
            if seq:
                try:
                    self.sequences.append(sequence(header, "".join(seq)))
                except ValueError as e:
                    self.errors.append(str(e))
            elif seq and not header:
                msg = "Sequence without header at end of file"
            if not self.sequences:
                msg = "No sequences found (empty or whitespace-only input)"
                if self.strict:
                    raise ValueError(msg)
                else:
                    self.errors.append(msg)

    def parse_file(self):
        if not self.path:
            raise ValueError("No file path provided")
        with self.path.open("r") as f:
            self._parse_lines(f)

    def parse_string(self, fasta_str: str):
        self._parse_lines(fasta_str.splitlines())   

    def get_report(self):
        lines =[]
        if not self.errors and not self.warnings:
            lines.append("Parser successful.")
        elif self.errors:
            lines.append("Parser failed with errors.")
        else:
            lines.append("Parser successful with warnings.")

        lines.append(f"Sequences parsed: {len(self.sequences)}")

        if self.errors:
            lines.append("\nErrors: ")
            for err in self.errors:
                lines.append(f"-{err}")
        if self.warnings:
            lines.append("\nWarnings: ")
            for war in self.warnings:
                lines.append(f"-{war}")
        return "\n".join(lines)
        
