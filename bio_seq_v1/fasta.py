from bio_seq_v1.stats import sequence
from pathlib import Path
from typing import Optional, Iterable


class FASTAParser:
    def __init__(
        self,
        path: Optional[str] = None,
        strict: bool = False,
        strict_file: bool = False,
        strict_seq: bool = False,
    ):
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
        return cls(path=path, strict=True, strict_file=True, strict_seq=True)

    def _strict_file_validate(self):
        if not self.path:
            raise ValueError("Strict mode requires a file path")
        if not self.path.exists():
            raise FileNotFoundError(f"{self.path} does not exist")
        if not self.path.is_file():
            raise IsADirectoryError(f"{self.path} is a directory")
        if self.path.stat().st_size == 0:
            raise ValueError("File is empty")

    def _validate_sequence(self, line: str, linenum: int):
        valid = set("ACGTNRYKMSWBDHV-.")
        invalid = set(line) - valid
        if invalid:
            msg = f"Invalid character(s) at line {linenum}: {''.join(sorted(invalid))}"
            raise ValueError(msg)
        

    def _parse_lines(self, lines: Iterable[str]):
        seq = []
        header = None

        for linenum, raw in enumerate(lines, start=1):
            line = raw.strip()

            if not line:
                if header is not None:
                    msg = f"Empty or whitespace-only sequence line at line {linenum}"
                    if self.strict or self.strict_seq:
                        raise ValueError(msg)
                    self.errors.append(msg)
                else:
                    self.warnings.append(
                        f"Empty or whitespace-only line ignored at line {linenum}"
                    )
                continue

            if line.startswith(">"):
                if header is not None:
                    if not seq:
                        msg = f"Header '{header}' has no sequence"
                        if self.strict or self.strict_seq:
                            raise ValueError(msg)
                        self.errors.append(msg)
                    else:
                        self.sequences.append(sequence(header, "".join(seq)))

                header = line[1:].strip()
                seq = []

                if not header:
                    msg = f"Empty FASTA header at line {linenum}"
                    if self.strict:
                        raise ValueError(msg)
                    self.errors.append(msg)

            else:
                if header is None:
                    msg = f"Sequence line before any header at line {linenum}"
                    if self.strict:
                        raise ValueError(msg)
                    self.errors.append(msg)
                    continue
                else:
                    try:
                        self._validate_sequence(line, linenum)
                    except ValueError as e:
                        if self.strict or self.strict_seq:
                            raise
                        self.errors.append(str(e))
                        continue

                seq.append(line.upper())

        if header is not None:
            if not seq:
                msg = f"Header '{header}' has no sequence at end of file"
                if self.strict or self.strict_seq:
                    raise ValueError(msg)
                self.errors.append(msg)
            else:
                self.sequences.append(sequence(header, "".join(seq)))

        if not self.sequences:
            msg = "No sequences found (empty or invalid input)"
            if self.strict:
                raise ValueError(msg)
            self.errors.append(msg)

    def parse_file(self):
        if not self.path:
            raise ValueError("No file path provided")
        with self.path.open("r") as f:
            self._parse_lines(f)

    def parse_string(self, fasta_str: str):
        self._parse_lines(fasta_str.splitlines())

    def get_report(self):
        lines = []

        if not self.errors and not self.warnings:
            lines.append("Parser successful.")
        elif self.errors:
            lines.append("Parser failed with errors.")
        else:
            lines.append("Parser successful with warnings.")

        lines.append(f"Sequences parsed: {len(self.sequences)}")

        if self.errors:
            lines.append("\nErrors:")
            for err in self.errors:
                lines.append(f"- {err}")

        if self.warnings:
            lines.append("\nWarnings:")
            for war in self.warnings:
                lines.append(f"- {war}")

        return "\n".join(lines)