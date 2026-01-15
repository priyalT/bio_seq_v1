import pytest
import os
from bio_seq_v1.stats import sequence
from bio_seq_v1.fasta import FASTAParser
from hypothesis import given, strategies as st

@given(st.text())
def test_fasta_parser_invalid_input_raises(path):
    try:
        parser = FASTAParser(path)
        if path:
            parser.parse_file()
        assert isinstance(parser.sequences, list)
        assert all(isinstance(s, sequence) for s in parser.sequences)
    except(FileNotFoundError, IsADirectoryError, ValueError):
        pass    

@given(st.text(alphabet="ACGTNRYKMSWBDHV", min_size=1))
def test_sequence_length_non_negative(seq_str):
    seq = sequence("id", seq_str)
    assert seq.sequence_length() >= 0

@given(st.text(alphabet="ACGTNRYKMSWBDHV", min_size=1))
def test_gc_content(seq_str):
    seq = sequence("id", seq_str)
    assert 0.0 <= seq.gc_content() <= 100

@given(st.text(alphabet = 'ACGTNRYKMSWBDHV', min_size=1))
def test_base_count(seq_str):
    seq = sequence("id", seq_str)
    counts = seq.base_count()
    length = seq.sequence_length()
    sum(counts.values()) == length

@given(st.text(alphabet='ACGTNRYKMSWBDHV', min_size = 1))
def test_reverse_complement(seq_str):
    seq = sequence("id", seq_str)
    rc1 = seq.rev_complement()
    seq_rc = sequence("id", rc1)
    rc2 = seq_rc.rev_complement()
    assert rc2 == seq_str.upper()

