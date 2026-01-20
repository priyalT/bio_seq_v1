from hypothesis import given, strategies as st
from hypothesis.strategies import sampled_from
from bio_seq_v1.stats import sequence
from bio_seq_v1.motif_search import MotifFinder
from bio_seq_v1.motif_search import Match

#exact motif matching
@given(motif_char=sampled_from(list(MotifFinder.IUPAC.keys())),base=sampled_from(list("ACGT")))
def test_exact_motif_matching(motif_char, base):
    finder = MotifFinder(k=1)
    expected = base in MotifFinder.IUPAC[motif_char]
    assert finder.char_match(motif_char, base) == expected



