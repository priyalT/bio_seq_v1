class Match():
    def __init__(self, seq_id, position, matched_seq, strand_attributes):
        self.seq_id = seq_id
        self.position = position
        self.matched_seq = matched_seq
        self.strand_attributes = strand_attributes

    def to_dict(self):
        return {
            'seq_id' : self.seq_id,
            'position' : self.position,
            'matched_seq' : self.matched_seq,
            'strand_attributes' : self.strand_attributes

        }
    


class MotifFinder():
    IUPAC = {
        'A': {'A'},
        'C': {'C'},
        'G': {'G'},
        'T': {'T'},
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'}
    }
    def __init__(self, dna: list[str], k: int):
        if not dna:
            raise ValueError("DNA list cannot be empty")
        if k <= 0:
            raise ValueError("Motif length must be positive")
        self.dna = dna
        self.k = k
        if not all(isinstance(seq, str) for seq in dna):
            raise TypeError("All DNA sequences must be strings")
        self._validate_sequences()

    def _convert_iupac_to_regex (self, motif: str):
        motif = motif.upper()
        regex_parts = []
        for char in motif:
            if char not in self.IUPAC:
                raise ValueError(f"Invalid IUPAC code: {char}")
            bases = self.IUPAC[char]
            if len(bases) == 1:
                regex_parts.append(next(iter(bases)))
            else:
                regex_parts.append(f"[{''.join(sorted(bases))}]")
        return "".join(regex_parts)
    
    def _validate_sequences(self):
        valid = set(self.IUPAC.keys())
        for seq in self.dna:
            if not set(seq).issubset(valid):
                raise ValueError(f"Invalid character in sequence: {seq}")
            
    def char_match(self, motif_char: str, base: str) -> bool:
        return base in self.IUPAC[motif_char]

    def mismatches(self, motif: str, kmer: str) -> int:
        if len(motif) != len(kmer):
            raise ValueError("Motif and k-mer must be of the same length")
        
        mismatches = 0
        for m,b in zip(motif, kmer):
            if not self.char_match(m,b):
                mismatches += 1
        return mismatches
    
    def kmer_generation(self, seq :str) -> list[str]:
        return [

        ]
    


        

    
