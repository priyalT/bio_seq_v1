from hypothesis import given, assume, strategies as st
from hypothesis.strategies import sampled_from
from bio_seq_v1.stats import sequence
from bio_seq_v1.motif_search import MotifFinder
from bio_seq_v1.motif_search import Match

#exact motif matching
@given(k=st.integers(min_value=1, max_value=5), seq = st.text(alphabet="ACGT", min_size=5))
def test_exact_motif_matching(seq, k):
    assume(len(seq) >= k)
    seq_obj = sequence("test_seq", seq)
    finder = MotifFinder(k=k)
    motif = seq[:k]
    motif_search = finder.search_single(seq_obj, motif, max_mismatches=0)
    for m in motif_search:
        assert m.matched_seq == motif

#motif match metadata
@given(seq_char=sampled_from(list("ACGTNRYSWKMBDHV")))
def test_motif_metadata(seq_char):
    seq_obj = sequence("test_seq", seq_char)
    finder = MotifFinder(k=1)
    motif_search = finder.search_both_strands(seq_obj, seq_char, mismatches=0)
    assert all(isinstance(match, Match) for match in motif_search)
    assert all(match.matched_seq == seq_obj.sequence for match in motif_search)
    assert all(match.matched_seq == seq_char for match in motif_search)

#IUPAC motif matching
@given(motif_char=sampled_from(list(MotifFinder.IUPAC.keys())),base=sampled_from(list("ACGT")))
def test_iupac_motif_matching(motif_char, base):
    finder = MotifFinder(k=1)
    expected = base in MotifFinder.IUPAC[motif_char]
    assert finder.char_match(motif_char, base) == expected

#both strands matching
@given(seq_char=sampled_from(list(MotifFinder.IUPAC.keys())))
def test_both_strands(seq_char):
    seq_obj = sequence("test_seq", seq_char)
    finder = MotifFinder(k=1)
    motif_search = finder.search_both_strands(seq_obj, seq_char, mismatches=0)
    for match in motif_search:
        assert match.strand_attributes in {"+", "-"}


