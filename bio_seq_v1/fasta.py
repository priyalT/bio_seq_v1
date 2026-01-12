"""
FASTA file parsing and input validation utilities.

This module provides helper functions to validate FASTA file paths
and parse FASTA-formatted files into sequence objects.
"""

from bio_seq_v1.stats import sequence
from pathlib import Path

class FASTAParser():
    def __init__(self, path: str, strict: bool = False):
        self.path = Path(path)
        self.sequences = []
        self.errors = []
        self.warnings = []
        self.strict = strict
        if self.strict:
            self._strict_validate()
        if self.path.exists():
            self._parse()

    @classmethod
    def strict(cls, path: str):
        return cls(path, strict=True)

    def _strict_validate(self):
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} does not exist")
        if not self.path.is_file():
            raise IsADirectoryError(f"{self.path} is a directory")
        if self.path.stat().st_size == 0:
            raise ValueError(f"File is empty.")
        
    def _validate_sequence(self, line, linenum):
        valid_chars = {"A","C","G","T","U","N","R","Y","S","W","K","M","B","D","H","V","-","."}
        invalid_chars = set(line) - valid_chars
        if invalid_chars:
            msg = f"Invalid characters at line {linenum}: {''.join(invalid_chars)}"
            if self.strict:
                raise ValueError(msg)
            else:
                self.errors.append(msg)
        
              
    def _parse(self):
        with self.path.open('r') as fasta:
            seq = []
            header = None
            for linenum, line in enumerate(fasta, start =1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is None and seq:
                        msg = f"Sequence found before header at line {linenum}: {seq}"
                        if self.strict:
                            raise ValueError(msg)
                        else:
                            self.errors.append(msg)
                    if header:
                        self.sequences.append(("".join(seq)))
                    header = line[1:]
                    seq = []
                else:
                    if not header:
                        msg = f"Sequence line before any header at line {linenum}"
                        if self.strict:
                            raise ValueError(msg)
                        else:
                            self.errors.append(msg)
                    self._validate_sequence(line, linenum)
                    seq.append(line)
        if seq:
            self.sequences.append("".join(seq))