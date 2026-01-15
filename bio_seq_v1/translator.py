
from bio_seq_v1.stats import sequence
genetic_code = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine (Start)
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic Acid
    "GAT": "D", "GAC": "D",
    # Glutamic Acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*"
}

start_codons = {"ATG"}
stop_codons = {"TAA", "TAG", "TGA"}

class Translator():
    def __init__(self, genetic_code=genetic_code, start_codons = start_codons, stop_codons = stop_codons):
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        self.stop_codons = stop_codons

    def translate(self, seq, frame):
        protein = []
        for i in range(frame, len(seq), 3):
            codon = seq[i:i+3]
            if len(codon) < 3:
                break
            codon = codon.upper()
            if codon not in self.genetic_code:
                protein.append("X")
            else:
                protein.append(self.genetic_code[codon])
        return "".join(protein)
    
    def translate_six_frames(self, seq):
        seq = seq.upper()
        dna = sequence("tmp_id",seq)
        reverse_complement_sequence = dna.rev_complement()
        
        results = {}

        for i in range(3):
            results[f"+{i+1}"] = self.translate(seq, i)

        for i in range(3):
            results[f"-{i+1}"] = self.translate(reverse_complement_sequence, i)
        
        return results

        