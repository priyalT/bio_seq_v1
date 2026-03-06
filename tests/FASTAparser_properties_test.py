import pytest
from pathlib import Path
from bio_seq_v1.fasta import FASTAParser
from hypothesis import given, strategies as st

TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"


@given(st.text(alphabet="ACGTNRYKMSWBDHV", min_size=1))
def test_fasta_preserves_sequences(seq):
    fasta = ""
    for i, nuc in enumerate(seq):
        fasta += f">seq{i}\n{nuc}\n"
    parser = FASTAParser(strict=False)
    parser.parse_string(fasta)
    expected = list(seq)
    assert [s.sequence for s in parser.sequences] == expected


@given(
    st.text(
        alphabet=st.characters(
            blacklist_characters=">ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.\n"
        ),
        min_size=1,
    )
)
def test_fasta_file_rejection(fasta_str):
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError):
        parser.parse_string(fasta_str)


@given(
    st.text(
        alphabet=st.characters(
            blacklist_characters=">ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.\n"
        ),
        min_size=1,
    )
)
def test_invalid_fasta_records_errors(fasta_str):
    parser = FASTAParser(strict=False)
    parser.parse_string(fasta_str)
    assert parser.errors


blacklist_chars = "ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.<"
blacklist_chars_lower = blacklist_chars.lower()
blacklist = blacklist_chars + blacklist_chars_lower


@given(st.text(alphabet=st.characters(blacklist_characters=blacklist), min_size=1))
def test_invalid_nucleotide(fasta_seq):
    fasta = f">seq1\n{fasta_seq}"
    parser = FASTAParser(strict_seq=False)
    parser.parse_string(fasta)
    assert any("invalid" in e.lower() or "empty" in e.lower() for e in parser.errors)


@given(st.text(alphabet=st.characters(blacklist_characters=blacklist), min_size=1))
def test_line_reporting(fasta_seq):
    fasta = f">seq1\n{fasta_seq}"
    parser = FASTAParser(strict_seq=False)
    parser.parse_string(fasta)
    assert any(
        ("invalid" in e.lower() or "empty" in e.lower() or "no sequence" in e.lower())
        and "line" in e.lower()
        for e in parser.errors
    )


@given(st.text(alphabet=" ", min_size=0))
def test_empty_files(fasta_seq):
    parser = FASTAParser(strict=False)
    parser.parse_string(fasta_seq)
    assert parser.errors
    assert any("empty" in e.lower() for e in parser.errors)


def test_strict_mode_sets_all_flags():
    parser = FASTAParser.strict_mode(str(DATA_DIR / "tiny.fasta"))
    assert parser.strict is True
    assert parser.strict_file is True
    assert parser.strict_seq is True


def test_strict_mode_parses_valid_file():
    parser = FASTAParser.strict_mode(str(DATA_DIR / "tiny.fasta"))
    parser.parse_file()
    assert len(parser.sequences) > 0


def test_strict_mode_rejects_nonexistent_file():
    with pytest.raises(FileNotFoundError):
        FASTAParser.strict_mode("/fake/path/doesnt_exist.fasta")


def test_strict_mode_rejects_empty_file(tmp_path):
    empty_file = tmp_path / "empty.fasta"
    empty_file.write_text("")
    with pytest.raises(ValueError, match="File is empty"):
        FASTAParser.strict_mode(str(empty_file))


def test_strict_file_rejects_directory(tmp_path):
    with pytest.raises(IsADirectoryError):
        FASTAParser(path=str(tmp_path), strict_file=True)


def test_strict_file_no_path():
    with pytest.raises(ValueError, match="requires a file path"):
        FASTAParser(path=None, strict_file=True)


def test_parse_file_no_path_raises():
    parser = FASTAParser()
    with pytest.raises(ValueError, match="No file path"):
        parser.parse_file()


def test_parse_string_no_header_adds_anonymous():
    parser = FASTAParser()
    parser.parse_string("ATGCGTACG")
    assert len(parser.sequences) == 1
    assert parser.sequences[0].id == "anonymous"
    assert parser.sequences[0].sequence == "ATGCGTACG"


def test_parse_string_with_header_keeps_it():
    parser = FASTAParser()
    parser.parse_string(">myseq\nATGCGTACG")
    assert parser.sequences[0].id == "myseq"


def test_empty_header_strict_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="Empty FASTA header"):
        parser.parse_string(">\nATGC")


def test_empty_header_nonstrict_records_error():
    parser = FASTAParser(strict=False)
    parser.parse_string(">\nATGC")
    assert any("empty" in e.lower() and "header" in e.lower() for e in parser.errors)


def test_sequence_before_header_strict_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="before any header"):
        parser.parse_string("ATGC\n>seq1\nGGGG")


def test_sequence_before_header_nonstrict_records_error():
    parser = FASTAParser(strict=False)
    parser.parse_string("ATGC\n>seq1\nGGGG")
    assert any("before any header" in e.lower() for e in parser.errors)
    assert len(parser.sequences) == 1
    assert parser.sequences[0].id == "seq1"


def test_consecutive_headers_strict_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="no sequence"):
        parser.parse_string(">seq1\n>seq2\nATGC")


def test_consecutive_headers_nonstrict_records_error():
    parser = FASTAParser(strict=False)
    parser.parse_string(">seq1\n>seq2\nATGC")
    assert any("no sequence" in e.lower() for e in parser.errors)
    assert len(parser.sequences) == 1
    assert parser.sequences[0].id == "seq2"


def test_last_header_no_sequence_strict_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="no sequence"):
        parser.parse_string(">seq1\nATGC\n>seq2")


def test_last_header_no_sequence_nonstrict():
    parser = FASTAParser(strict=False)
    parser.parse_string(">seq1\nATGC\n>seq2")
    assert any("no sequence" in e.lower() for e in parser.errors)
    assert len(parser.sequences) == 1


def test_strict_seq_invalid_char_raises():
    parser = FASTAParser(strict_seq=True)
    with pytest.raises(ValueError, match="Invalid character"):
        parser.parse_string(">seq1\nATG123XYZ")


def test_strict_seq_valid_chars_passes():
    parser = FASTAParser(strict_seq=True)
    parser.parse_string(">seq1\nATGCNRYSWKM")
    assert len(parser.sequences) == 1


def test_empty_line_before_header_produces_warning():
    parser = FASTAParser(strict=False)
    parser.parse_string("\n\n>seq1\nATGC")
    assert any("ignored" in w.lower() for w in parser.warnings)


def test_empty_line_in_sequence_strict_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="Empty or whitespace"):
        parser.parse_string(">seq1\nATGC\n\nGGGG")


def test_empty_line_in_sequence_nonstrict_records_error():
    parser = FASTAParser(strict=False)
    parser.parse_string(">seq1\nATGC\n\nGGGG")
    assert any("empty" in e.lower() for e in parser.errors)


def test_get_report_success():
    parser = FASTAParser()
    parser.parse_string(">seq1\nATGC")
    report = parser.get_report()
    assert "successful" in report.lower()
    assert "Sequences parsed: 1" in report


def test_get_report_with_errors():
    parser = FASTAParser(strict=False)
    parser.parse_string("ATGC\n>seq1\nGGGG")
    report = parser.get_report()
    assert "errors" in report.lower()
    assert "Errors:" in report


def test_get_report_with_warnings():
    parser = FASTAParser(strict=False)
    parser.parse_string("\n>seq1\nATGC")
    report = parser.get_report()
    if parser.warnings and not parser.errors:
        assert "warnings" in report.lower()


def test_get_report_shows_sequence_count():
    parser = FASTAParser()
    parser.parse_string(">s1\nATGC\n>s2\nGGGG\n>s3\nCCCC")
    report = parser.get_report()
    assert "Sequences parsed: 3" in report


def test_non_ascii_in_sequence_raises():
    parser = FASTAParser(strict=True)
    with pytest.raises(ValueError, match="Invalid character"):
        parser.parse_string(">seq1\nATGC™©")


def test_non_ascii_nonstrict_records_error():
    parser = FASTAParser(strict=False)
    parser.parse_string(">seq1\nATGC™©")
    assert any("invalid" in e.lower() for e in parser.errors)
