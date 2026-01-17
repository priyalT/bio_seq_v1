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