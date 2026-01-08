#this is the file with the statistics involved. we will be calculating length of sequences
valid = "ACGTUNRYSWKMBDHV-."
revcomp_dict = {
    "A":"T", "T":"A", "G":"C", "C":"G", "U":"A",
    "R":"Y", "Y":"R", "S":"S", "W":"W",
    "K":"M", "M":"K", "B":"V", "D":"H",
    "H":"D", "V":"B", "N":"N"
}

def print_sequence_lengths(sequences):
    length_seq = []
    for seq in sequences:
        length_seq.append({'id': seq['id'], "length of sequence": len(seq['sequence'])})
    return length_seq


def base_count(sequences):
    all_counts = []
    for i in sequences:
        counts = {base: 0 for base in valid}
        for j in i['sequence']:
            if j in counts:
                counts[j] += 1
        counts["id"] = i["id"]
        all_counts.append(counts)
    return all_counts

def gc_content(sequences):
    base = base_count(sequences) #outputs a list of dictionaries
    gc_list = []
    for base_counts in base:
        total_bases = sum(base_counts[b] for b in valid)
        g = base_counts["G"]
        c = base_counts["C"]
        gc_percent = ((g+c)/total_bases)*100 if total_bases > 0 else 0
        gc_list.append({"id": base_counts["id"], "GC%": gc_percent})
    return gc_list

def rev_complement(sequences):
    reverse_complement = []
    for seq in sequences:
        reverse = seq['sequence'][::-1]
        complement = "".join(revcomp_dict[d] for d in reverse)
        reverse_complement.append({"id": seq["id"], "reverse complement": complement})
    return reverse_complement


