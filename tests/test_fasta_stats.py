#test fasta parser
#test gc content
#test reverse complement
#test total bases
import pytest
from bio_seq_v1.stats import (print_sequence_lengths, base_count, gc_content, rev_complement)

fasta_seq = "tests/data/tiny.fasta"

def test_fasta_parser():
    result = print_sequence_lengths(fasta_seq)
     

