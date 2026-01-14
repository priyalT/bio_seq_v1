import pytest
from bio_seq_v1.stats import sequence
from bio_seq_v1.fasta import fasta_parser
fasta_seq = "tests/data/tiny.fasta"
single_seq = "tests/data/single.fasta"

#test fasta parser
def test_fasta_parser():
    parsed = fasta_parser(fasta_seq)
    assert isinstance(parsed, list)
    assert all(isinstance(s, sequence) for s in parsed)
    first_seq = parsed[0]
    assert first_seq.sequence == "ATGCGTACGTAGCTAGTTAGCGATCGGGGCTAGCTAGCTAGCTAG"
    assert first_seq.id == "seq1"

#test length of sequence
def test_seq_length():
    parsed = fasta_parser(fasta_seq)
    assert parsed[0].sequence_length() == 45
    
#test gc content     
def test_gc_content():
    parsed = fasta_parser(fasta_seq)
    gc_content = parsed[0].gc_content()
    assert round(gc_content, 2) == 53.33  

#test total bases
def test_total_bases():
    parsed = fasta_parser(fasta_seq)
    base_count = parsed[0].base_count()
    assert base_count["A"] == 10
    assert base_count["C"] == 9
    assert base_count["G"] == 15 
    assert base_count["T"] == 11

#test reverse complement
def test_reverse_complement():
    parsed = fasta_parser(fasta_seq)
    reverse_comp = parsed[0].rev_complement()
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
    parsed = fasta_parser(fasta_seq)
    assert len(parsed) == 3

def test_single_fasta_entries():
    parsed = fasta_parser(single_seq)
    assert len(parsed) == 1

def test_ambig_codes():
    parsed = fasta_parser(single_seq)
    assert parsed[0].rev_complement() == "KDMAAATTTCCCGGGBYACCCGGGTTTAAACCC"

def test_empty_sequence():
    seq = sequence("id", "ATATAT")
    assert seq.gc_content() == 0.0