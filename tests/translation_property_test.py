import pytest
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from hypothesis import given, strategies as st
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY*X")
stop_codons = {"TAA", "TAG", "TGA"}

@given(st.text(alphabet = "ACGT", min_size = 3))
def test_genetic_code_compliance(seq):
    translator = Translator()
    protein = translator.translate(seq, frame=0)
    for aa in protein:
        assert aa in valid_amino_acids

seq_strategy = st.text(alphabet="ACGT", min_size=9).filter(lambda s: len(set(s)) > 1)
@given(seq_strategy)
def test_frame_offset_correctness(seq):
    translator = Translator()
    p0 = translator.translate(seq, frame = 0)
    p1 = translator.translate(seq, frame = 1)
    p2 = translator.translate(seq, frame = 2)
    assert len({p0, p1, p2}) > 1
    
@given(st.text(alphabet = "ACGT", min_size = 1))
def test_incomplete_codon_handling(seq):
    translator = Translator()
    protein = translator.translate(seq, frame=0)
    expected_len = len(seq) // 3
    assert len(protein) == expected_len


@given(st.text(alphabet = "ACGT", min_size = 3))
def test_stop_codon_representation(seq):
    translator = Translator()
    protein = translator.translate(seq, frame = 0)
    for i in range(0, len(seq)-2, 3):
        codons = seq[i:i+3]
    for codon, aa in zip(codons, protein):
        if codon in stop_codons:
            assert aa == "*"

@given(st.text(alphabet = "ACGTNRYSWKMBDHV", min_size = 3))
def test_ambiguous_code_handling(seq):
    translator = Translator()
    protein = translator.translate(seq, frame = 0)
    codons = [
        seq[i:i+3]
        for i in range(0, len(seq)-2, 3)
    ]
    for codon, aa in zip(codons, protein):
        if codon not in translator.genetic_code:
            assert aa == "X"
        else:
            assert aa in valid_amino_acids

@given(st.text(alphabet = "ACGTNRYSWKMBDHV", min_size = 3))
def test_six_frame_translation_completeness(seq):
    translator = Translator()
    dna = sequence("tmp", seq)
    rev = dna.rev_complement()
    six_frame = translator.translate_six_frames(seq)

    for frame in range(3):
        single = translator.translate(seq, frame = frame)
        assert single == six_frame[f"+{frame+1}"]

    for frame in range(3):
        single = translator.translate(rev, frame=frame)
        assert single == six_frame[f"-{frame+1}"]



