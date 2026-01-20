from bio_seq_v1.stats import sequence
from typing import List
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
    def __init__(self, k: int):
        if k <= 0:
            raise ValueError("Motif length must be positive")
        self.k = k
    
    def _validate_seq_string(self, seq: str):
        valid_bases = set().union(*self.IUPAC.values())
        if not set(seq).issubset(valid_bases):
            raise ValueError("Invalid characters in sequence")

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

    def kmer_generation(self, seq: str):
        k = self.k
        kmers = []
        for base in range(len(seq) - k + 1):
            kmers.append((seq[base:base+k], base))
        return kmers

    def search_single(self, seq_obj: sequence, motif: str, max_mismatches: int):
        if len(motif) != self.k:
            raise ValueError("Motif length must match k")
        motif = motif.upper()
        if not set(motif).issubset(self.IUPAC):
            raise ValueError(f"Invalid IUPAC character in motif: {motif}")
        seq = seq_obj.sequence.upper()
        motif = motif.upper()
        self._validate_seq_string(seq)

        matches = []
        for kmer, pos in self.kmer_generation(seq):
            if self.mismatches(motif, kmer) <= max_mismatches:
                matches.append(
                    Match(seq_obj.id, pos, kmer, None)
                )
        return matches
    
    def search_fasta(self, fasta_sequences: List[sequence], motif: str, mismatches: int):
        all_matches = []
        for seq_obj in fasta_sequences:
            all_matches.extend(
                self.search_single(seq_obj, motif, mismatches)
            )
        return all_matches
    
    def search_both_strands(self, seq_obj: sequence, motif: str, mismatches: int):
        matches = []
        for match in self.search_single(seq_obj, motif, mismatches):
            match.strand_attributes = "+"
            matches.append(match)


        rev_seq = seq_obj.rev_complement()
        rev_obj = sequence(seq_obj.id, rev_seq)
        seq_len = len(seq_obj.sequence)
        for match in self.search_single(rev_obj, motif, mismatches):
            original_pos = seq_len - match.position - self.k
            matches.append(
                Match(seq_id = seq_obj.id,
                      position= original_pos,
                      matched_seq= match.matched_seq,
                      strand_attributes="-")
            )
        return matches