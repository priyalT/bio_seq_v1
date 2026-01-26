import pytest
from hypothesis import given, strategies as st
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector


def count_orfs_in_protein(protein: str) -> int:
    count = 0
    i = 0
    while i < len(protein):
        if protein[i] == "M":
            j = i + 1
            while j < len(protein) and protein[j] != "*":
                j += 1
            if j < len(protein) and protein[j] == "*":
                count += 1
        i += 1
    return count

@given(st.text(alphabet="ACGT", min_size=3))
def test_orf_detection_completeness(dna):
    detector = ORFDetector(min_length=0)
    translator = Translator()

    detected_orfs = detector.find_orfs(dna)
    frames = translator.translate_six_frames(dna)

    expected_count = 0
    for protein in frames.values():
        expected_count += count_orfs_in_protein(protein)

    assert len(detected_orfs) == expected_count

@given(st.text(alphabet="ACGT", min_size=3))
def test_orf_length_filtering(dna):
    detector = ORFDetector(min_length=3)
    detected_orfs = detector.find_orfs(dna)
    for orf in detected_orfs:
        assert orf.length >= detector.min_length

@given(st.text(alphabet="ACGT", min_size=3))
def test_orf_metadata_completeness(dna):
    detector = ORFDetector(min_length=3)
    detected_orfs = detector.find_orfs(dna)
    for orf in detected_orfs:
        hasattr(orf, "start")
        hasattr(orf, "end")
        hasattr(orf, "frame")
        hasattr(orf, "strand")
        hasattr(orf, "protein")
        hasattr(orf, "length")
        
        assert isinstance(orf.start, int)
        assert isinstance(orf.end, int)
        assert isinstance(orf.frame, int)
        assert isinstance(orf.strand, str)
        assert isinstance(orf.protein, str)
        assert isinstance(orf.length, int)

        assert orf.length == orf.end - orf.start + 1
        assert orf.start >= 0
        assert orf.end >= orf.start
        assert orf.frame in {0,1,2}
        assert orf.strand in {"+", "-"}
        assert len(orf.protein) > 0 
        assert "*" not in orf.protein
        assert orf.protein[0] == "M"

        assert (orf.length % 3) == 0

        d = orf.to_dict()
        expected_keys = {"start", "end", "frame", "strand", "protein", "length"}
        assert set(d.keys()) == expected_keys
        assert d["start"] == orf.start
        assert d["length"] == orf.length

@given(st.text(alphabet="ACGT", min_size = 3))
def test_orf_overlap(dna):
    detector = ORFDetector(min_length = 3)
    detected_orfs = detector.find_orfs(dna)
    overlapping_pairs = detector.overlapping_orfs(detected_orfs)
    for orf1, orf2 in overlapping_pairs:
        assert orf1 in detected_orfs
        assert orf2 in detected_orfs
        assert orf1.start <= orf2.end
        assert orf2.start <= orf1.end