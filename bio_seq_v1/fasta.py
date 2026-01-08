#defining our parser function
def fasta_parser(filename):
    all_seq =[]
    current_id = None
    current_seq = []
    valid = set("ACGTUNRYSWKMBDHV-.")

    with open(filename, 'r') as fasta: 
        for line in fasta: 
            line = line.strip() 

            if not line:
                continue 

            if line.startswith(">"): 
                if current_id is not None: 
                    all_seq.append({"id": current_id, "sequence": "".join(current_seq)}) 
                current_id = line[1:] 
                current_seq = [] 
            else:
                seq_line = line.upper()
                if all(base in valid for base in seq_line):
                     current_seq.append(seq_line)

        if current_id is not None: 
                    all_seq.append({"id": current_id, "sequence": "".join(current_seq)}) 
                    current_seq = []
    return all_seq            