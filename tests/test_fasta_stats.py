import pytest
from bio_seq_v1.stats import sequence
from bio_seq_v1.fasta import fasta_parser
fasta_seq = "tests/data/tiny.fasta"

#test fasta parser
def test_fasta_parser():
    parsed = fasta_parser(fasta_seq)
    assert isinstance(parsed, list)
    assert all(isinstance(s, sequence) for s in parsed)
    first_seq = parsed[0]
    assert first_seq.sequence == "ATGCGTACGTAGCTAGTTAGCGATCGGGGCTAGCTAGCTAGCTAG"
    assert first_seq.id == "seq1"

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