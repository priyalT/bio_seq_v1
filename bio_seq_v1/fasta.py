# FASTA file parsing and input validation utilities
#importing class sequence from stats.py
#importing Path from pathlib to validate the input path and raise errors
from bio_seq_v1.stats import sequence
from pathlib import Path

#defining the function for validating input path
def validate_input_path(filename):
    path = Path(filename) #setting the filename's path to a variable
    if not path.exists(): #if path does not exist
        raise FileNotFoundError(f"File not found: {path}") #raise the filenotfound error
    if not path.is_file(): #if path is not a file (maybe the path provided is that of a directory)
         raise IsADirectoryError(f"Not a file: {path}") #raise isadirectoryerror
    if path.stat().st_size == 0: #if the file size is zero (empty file)
         raise ValueError(f"File is empty.") #raise Valueerror
    return path 

#defining the fasta parser function with input of the given file
def fasta_parser(filename):
    all_seq =[] #creating an empty list to store sequence objects (id + sequence)
    current_id = None #setting current id to none
    current_seq = [] #creating an empty list for only sequences
    valid = set("ACGTUNRYSWKMBDHV-.") #creating a set for valid bases
    path = validate_input_path(filename) #a variable for the output of validate_input_path function (described above)
    
    with path.open('r') as fasta: #opening the file and reading it
        for linenum, line in enumerate(fasta, start=1): #running the loop on lines and numerating each line
            line = line.strip() #removing whitespace if any

            if not line: #if the line is empty (i.e no whitespaces), skip it    
                continue 

            if line.startswith(">"): #if the line starts with ">", it is probably the header and the indication for starting a new sequence
                if current_id is not None: #if there is information stored in current_id (i.e name of the sequence on which the loop is running)
                    all_seq.append(sequence(current_id, "".join(current_seq))) ##append the completed sequence object to the list
                current_id = line[1:] #then, let current_id store the id (ignoring the ">")
                current_seq = [] #reset sequence accumulator for the new record 
            else:
                if current_id is None: #if current_id is none (i.e loop is at the sequence, yet id is none)
                     raise ValueError(f"Sequence data found before header on line {linenum}") #raise error
                seq_line = line.upper() #set the line of sequence to upper case
                if all(base in valid for base in seq_line): #if all the bases in seq_line also exist in the valid set
                     current_seq.append(seq_line) #then append it to current_seq
                else: 
                     raise ValueError(f"Invalid character in sequence '{current_id}':'{seq_line}' " #otherwise raise error and inform user about the sequence, id, and line
                                      f"on line: {linenum}")

        if current_id is not None: #at the end of the loop, save the last id, sequence
                    all_seq.append(sequence(current_id, "".join(current_seq))) 

    return all_seq  #return the list of sequence objects