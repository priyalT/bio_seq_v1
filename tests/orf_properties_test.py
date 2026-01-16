import pytest
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector
from hypothesis import given, strategies as st

#ORF detection completeness
@given(st.text(alphabet=""))