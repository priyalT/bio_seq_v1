import pytest
from bio_seq_v1.fasta import FASTAParser
from hypothesis import given, strategies as st

@given(st.text(alphabet="ACGTNRYKMSWBDHV", min_size = 1))
def test_fasta_preserves_sequences(seq):
    fasta = ""
    for i, nuc in enumerate(seq):
        fasta += f">seq{i}\n{nuc}\n"
    parser = FASTAParser(strict=False)
    parser.parse_string(fasta)
    expected = list(seq)
    assert [s.sequence for s in parser.sequences] == expected

@given(st.text(alphabet=st.characters(blacklist_characters=">ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.\n"), min_size = 1))
def test_fasta_file_rejection(fasta_str):
    parser = FASTAParser(strict = True)
    with pytest.raises(ValueError):
        parser.parse_string(fasta_str)

@given(st.text(alphabet=st.characters(blacklist_characters=">ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.\n"), min_size = 1))
def test_invalid_fasta_records_errors(fasta_str):
    parser = FASTAParser(strict = False)
    parser.parse_string(fasta_str)
    assert parser.errors

blacklist_chars = "ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.<"
blacklist_chars_lower = blacklist_chars.lower()
blacklist = blacklist_chars + blacklist_chars_lower
@given(st.text(alphabet=st.characters(blacklist_characters=blacklist), min_size = 1))
def test_invalid_nucleotide(fasta_seq):
    fasta = f">seq1\n{fasta_seq}"
    parser = FASTAParser(strict_seq=False)
    parser.parse_string(fasta)
    assert any(
        "invalid" in e.lower() or "empty" in e.lower()
        for e in parser.errors
    )

@given(st.text(alphabet=st.characters(blacklist_characters=blacklist), min_size = 1))
def test_line_reporting(fasta_seq):
    fasta = f">seq1\n{fasta_seq}"
    parser = FASTAParser(strict_seq=False)
    parser.parse_string(fasta)
    assert any(
        ("invalid" in e.lower() or "empty" in e.lower() or "no sequence" in e.lower()) 
        and "line" in e.lower()
        for e in parser.errors)

@given(st.text(alphabet=" ", min_size = 0))
def test_empty_files(fasta_seq):
    parser = FASTAParser(strict = False)
    parser.parse_string(fasta_seq)
    assert parser.errors
    assert any("empty" in e.lower() for e in parser.errors)