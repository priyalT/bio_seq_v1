import pytest
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector
from hypothesis import given, strategies as st

#ORF detection completeness
@given(st.text(alphabet="ACGTNRYSWKMBDHV", min_size=3))
def test_orf_detection_completeness(seq):
    orf = ORFDetector(min_length=0)
    detected_orfs = orf.find_orfs(seq)
    frames = Translator().translate_six_frames(seq)
    expected = 0