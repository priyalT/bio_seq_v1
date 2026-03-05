import warnings


class sequence:
    """
    Represents a biological sequence (DNA/RNA) with utility methods
    for analysis such as length, base count, GC content, and reverse complement.

    Attributes:
        id (str): Identifier for the sequence.
        sequence (str): Uppercase string of sequence bases.
    """

    revcomp_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        "-": "-",
        ".": ".",
    }

    valid = "ACGTUNRYSWKMBDHV-."

    def __init__(self, id, sequence):
        """
        Initialize a Sequence object.

        Args:
            id (str): Identifier for the sequence.
            sequence (str): Sequence string containing valid bases.

        Raises:
            ValueError: If the sequence is empty or contains invalid characters.
        """

        if not sequence:
            raise ValueError(f"Sequence for ID '{id}' is empty")
        invalid_chars = set(sequence.upper()) - set(self.valid)
        if invalid_chars:
            raise ValueError(
                f"Sequence '{id}' contains invalid characters: {invalid_chars}"
            )
        self.id = id
        self.sequence = sequence.upper()

    def to_dict(self):
        return {
            "id": self.id,
            "length": self.sequence_length(),
            "gc content": self.gc_content(),
            "sequence": self.sequence,
            "reverse complement": self.rev_complement(),
        }

    def __getitem__(self, index):
        return self.sequence[index]

    def sequence_length(self):
        """
        Return the length of the sequence.

        Returns:
            int: Number of bases in the sequence.
        """

        return len(self.sequence)

    def base_count(self):
        """
        Count the occurrences of each valid base in the sequence.

        Returns:
            dict: Dictionary with bases as keys and their counts as values.

        Notes:
            If an invalid character is present, a warning is issued.
        """
        counts = {b: 0 for b in self.valid}
        for b in self.sequence:
            if b in counts:
                counts[b] += 1
            else:
                warnings.warn(
                    f"Invalid character in id: '{self.id}' and sequence: '{self.sequence}"
                )
        return counts

    def gc_content(self):
        """
        Calculate the GC content percentage of the sequence.

        Returns:
            float: GC content as a percentage.

        Raises:
            ValueError: If the sequence contains no valid bases for calculation.
        """
        if not self.sequence:
            return 0.0
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        total = self.sequence_length()
        return ((g + c) / total) * 100

    def rev_complement(self):
        """
        Compute the reverse complement of the sequence.

        Returns:
            str: Reverse complement string.

        Raises:
            ValueError: If the sequence contains invalid bases not in revcomp_dict.
        """
        reverse = self.sequence[::-1]
        try:
            complement = "".join(self.revcomp_dict[b] for b in reverse)
        except KeyError as e:
            raise ValueError(f"Invalid base for reverse complement: {e.args[0]}")
        return complement
