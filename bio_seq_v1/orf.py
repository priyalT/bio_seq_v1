from bio_seq_v1.translator import Translator

class ORF():
    def __init__(self, seq_id, start, end, frame, strand, protein):
        if start < 0 or end < start:
            raise ValueError("Invalid ORF coordinates")
        if frame not in (0, 1, 2):
            raise ValueError("Frame must be 0, 1 or 2")
        if strand not in ("+", "-"):
            raise ValueError ("Strand must be '+' or '-'")
        self.seq_id = seq_id
        self.start = start
        self.end = end
        self.frame = frame
        self.strand = strand
        self.protein = protein
        self.length = end - start + 1

    def to_dict(self):
        return {
            'seq_id': self.seq_id,
            'start': self.start,
            'end': self.end,
            'frame': self.frame,
            'strand': self.strand,
            'length': self.length,
            'protein': self.protein
        }
    def __eq__(self, other):
        if not isinstance(other, ORF):
            return NotImplemented
        return (
            self.seq_id == other.seq_id and
            self.start == other.start and
            self.end == other.end and
            self.frame == other.frame and
            self.strand == other.strand and
            self.protein == other.protein
        )

class ORFDetector():
    START_CODONS = {"ATG"}
    STOP_CODONS = {"TAA", "TAG", "TGA"}

    def __init__(self, min_length = 0):
        if min_length < 0:
            raise ValueError("min_length must be non-negative")
        self.translator = Translator()
        self.min_length = min_length

    def find_orfs(self, seq):
        orfs = []
        frames = self.translator.translate_six_frames(seq)
        seq_id = seq.id
        seq_len = len(seq.sequence)
        for frame_label, protein in frames.items():
            strand = frame_label[0]
            frame = int(frame_label[1]) - 1
            aa_index = 0
            while aa_index < len(protein):
                if protein[aa_index] == "M":
                    start_aa = aa_index
                    scan = aa_index + 1

                    while scan < len(protein) and protein[scan] != "*":
                        scan += 1

                    if scan < len(protein):  
                        end_aa = scan - 1
                        dna_start = frame + start_aa * 3
                        dna_end = frame + end_aa * 3 + 2
                        prot_seq = protein[start_aa:scan]
                        
                        if strand == "-":
                            rc_start = dna_start
                            rc_end = dna_end
                            dna_start = seq_len - rc_end - 1
                            dna_end = seq_len - rc_start - 1


                        if len(prot_seq) * 3 >= self.min_length:
                            orfs.append(
                                ORF(seq_id, dna_start, dna_end, frame, strand, prot_seq)
                            )

                aa_index += 1

        return orfs
        
    def overlapping_orfs(self, orfs: list):
        if not isinstance(orfs, list):
            raise TypeError("ORFs must be a list of ORF objects")
        
        overlap = []
        for i in range(len(orfs)):
            if not isinstance(orfs[i], ORF):
                raise TypeError("List must contain only ORF objects")
            orf1 = orfs[i]
            for j in range(i+1, len(orfs)):
                if not isinstance(orfs[j], ORF):
                    raise TypeError("List must contain only ORF objects")
                orf2 = orfs[j]
                if orf1.strand != orf2.strand:
                    continue
                if orf1.start <= orf2.end and orf2.start <= orf1.end:
                    overlap.append((orf1, orf2))
        return overlap