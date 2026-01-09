#test reverse complement
#test total bases
import pytest
from bio_seq_v1.stats import sequence
from bio_seq_v1.fasta import fasta_parser
fasta_seq = "tests/data/tiny.fasta"
parsed = fasta_parser(fasta_seq)
#test fasta parser
def test_fasta_parser():
    assert isinstance(parsed, list)
    assert all(isinstance(s, sequence) for s in parsed)
    first_seq = parsed[0]
    assert first_seq.sequence == "ATGCGTACGTAGCTAGTTAGCGATCGGGGCTAGCTAGCTAGCTAG"
    assert first_seq.id == "seq1"

#test gc content     
def test_gc_content():
    gc_content = parsed[0].gc_content()
    assert round(gc_content, 2) == 53.33    