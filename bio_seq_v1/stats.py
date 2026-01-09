#this is the file with the statistics involved. we will be calculating length of sequences

class sequence():
    revcomp_dict = {
    "A":"T", "T":"A", "G":"C", "C":"G", "U":"A",
    "R":"Y", "Y":"R", "S":"S", "W":"W",
    "K":"M", "M":"K", "B":"V", "D":"H",
    "H":"D", "V":"B", "N":"N"}
    
    valid = "ACGTUNRYSWKMBDHV-."

    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence.upper()

    def print_sequence_lengths(self):
        return len(self.sequence)
    
    def base_count(self):
            counts = {b:0 for b in self.valid}
            for b in self.sequence:
                if b in counts:
                    counts[b] += 1
            return counts
    
    def gc_content(self):
        g = self.sequence.count("G")
        c = self.sequence.count("C")
        total = sum(self.sequence.count(b) for b in self.valid)
        return ((g+c)/total)*100


    def rev_complement(self):
            reverse = self.sequence[::-1]
            complement = "".join(self.revcomp_dict[b] for b in reverse)
            return complement


