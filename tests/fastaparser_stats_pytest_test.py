import pytest
from bio_seq_v1.stats import sequence
from bio_seq_v1.fasta import FASTAParser
fasta_seq = "/Users/priyaltripathi/bio_seq_v1/tests/data/tiny.fasta"
single_seq = "tests/data/single.fasta"

def test_fasta_parser():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    sequences = parser.sequences
    assert isinstance(sequences, list)
    assert all(isinstance(s, sequence) for s in sequences)  
    first_seq = sequences[0].sequence  
    assert first_seq == "ATGCGTACGTAGCTAGTTAGCGATCGGGGCTAGCTAGCTAGCTAG"

#test length of sequence
def test_seq_length():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    seq_obj = parser.sequences[0]
    assert seq_obj.sequence_length() == 45
    

#test gc content     
def test_gc_content():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    seq_obj = parser.sequences[0]
    gc_content = seq_obj.gc_content()
    assert round(gc_content, 2) == 53.33  

#test total bases
def test_total_bases():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    seq_obj = parser.sequences[0]
    base_count = seq_obj.base_count()
    assert base_count["A"] == 10
    assert base_count["C"] == 9
    assert base_count["G"] == 15 
    assert base_count["T"] == 11

#test reverse complement
def test_reverse_complement():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    seq_obj = parser.sequences[0]
    reverse_comp = seq_obj.rev_complement()
    assert reverse_comp == "CTAGCTAGCTAGCTAGCCCCGATCGCTAACTAGCTACGTACGCAT"

#testing edge cases
#test lowercase input
def test_lowercase_input():
    seq = sequence("id", "ctga")
    assert seq.rev_complement() == "TCAG"

def test_single_base():
    seq = sequence("id", "G")
    assert seq.gc_content() == 100.0

def test_multiple_fasta_entries():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    assert len(parser) == 3

def test_single_fasta_entries():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    assert len(parser) == 1

def test_ambig_codes():
    parser = FASTAParser(fasta_seq)
    parser.parse_file()
    seq_obj = parser.sequences[0]
    assert seq_obj.rev_complement() == "KDMAAATTTCCCGGGBYACCCGGGTTTAAACCC"

def test_empty_sequence():
    seq = sequence("id", "ATATAT")
    assert seq.gc_content() == 0.0