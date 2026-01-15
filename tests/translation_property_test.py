import pytest
from bio_seq_v1.stats import sequence
from hypothesis import given, strategies as st

#Genetic code translation compliance
@given(st.text(alphabet = "ACGT", min_size = 1))
def test_genetic_code_compliance(seq):
    

