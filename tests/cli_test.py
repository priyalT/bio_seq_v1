import pytest
from click.testing import CliRunner
from bio_seq_v1.cli import (
    main,
    print_sequence_lengths_formatted,
    print_gc_content_table,
    print_revcomp,
    print_base_count,
    print_summary,
)
from bio_seq_v1.stats import sequence
from pathlib import Path

TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"
TINY_FASTA = str(DATA_DIR / "tiny.fasta")
SINGLE_FASTA = str(DATA_DIR / "single.fasta")


@pytest.fixture
def runner():
    return CliRunner()


def test_main_help(runner):
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "BioSeq" in result.output


def test_stats_summary_from_file(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA])
    assert result.exit_code == 0
    assert "seq1" in result.output


def test_stats_length_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--length"])
    assert result.exit_code == 0
    assert "Length" in result.output


def test_stats_gc_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--gc"])
    assert result.exit_code == 0
    assert "GC" in result.output or "%" in result.output


def test_stats_revcomp_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--revcomp"])
    assert result.exit_code == 0
    assert "reverse complement" in result.output.lower()


def test_stats_basecount_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--basecount"])
    assert result.exit_code == 0
    assert "A" in result.output


def test_stats_summary_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--summary"])
    assert result.exit_code == 0
    assert "Length" in result.output
    assert "GC" in result.output or "%" in result.output


def test_stats_from_string(runner):
    result = runner.invoke(main, ["stats", "--string", ">seq1\nATGCATGC"])
    assert result.exit_code == 0
    assert "seq1" in result.output


def test_stats_no_input_raises(runner):
    result = runner.invoke(main, ["stats"])
    assert result.exit_code != 0


def test_stats_both_file_and_string_raises(runner):
    result = runner.invoke(
        main, ["stats", "--file", TINY_FASTA, "--string", ">s\nATGC"]
    )
    assert result.exit_code != 0


def test_stats_nonexistent_file(runner):
    result = runner.invoke(main, ["stats", "--file", "/fake/file.fasta"])
    assert result.exit_code != 0 or "failed" in result.output.lower()


def test_stats_export_csv(runner, tmp_path):
    output = str(tmp_path / "out.csv")
    result = runner.invoke(
        main,
        [
            "stats",
            "--file",
            TINY_FASTA,
            "--length",
            "--output",
            output,
            "--format",
            "csv",
        ],
    )
    assert result.exit_code == 0
    assert Path(output).exists()
    assert "saved" in result.output.lower()


def test_stats_export_json(runner, tmp_path):
    output = str(tmp_path / "out.json")
    result = runner.invoke(
        main,
        ["stats", "--file", TINY_FASTA, "--gc", "--output", output, "--format", "json"],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_stats_export_tsv(runner, tmp_path):
    output = str(tmp_path / "out.tsv")
    result = runner.invoke(
        main,
        [
            "stats",
            "--file",
            TINY_FASTA,
            "--basecount",
            "--output",
            output,
            "--format",
            "tsv",
        ],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_stats_strict_flag(runner):
    result = runner.invoke(main, ["stats", "--file", TINY_FASTA, "--strict"])
    assert result.exit_code == 0


def test_stats_single_fasta(runner):
    result = runner.invoke(main, ["stats", "--file", SINGLE_FASTA])
    assert result.exit_code == 0


def test_translate_default_frame(runner):
    result = runner.invoke(main, ["translate", "--file", TINY_FASTA])
    assert result.exit_code == 0
    assert "Frame 0" in result.output


def test_translate_specific_frame(runner):
    result = runner.invoke(main, ["translate", "--file", TINY_FASTA, "--frame", "1"])
    assert result.exit_code == 0
    assert "Frame 1" in result.output


def test_translate_frame_2(runner):
    result = runner.invoke(main, ["translate", "--file", TINY_FASTA, "--frame", "2"])
    assert result.exit_code == 0
    assert "Frame 2" in result.output


def test_translate_six_frames(runner):
    result = runner.invoke(main, ["translate", "--file", TINY_FASTA, "--six-frames"])
    assert result.exit_code == 0
    assert "Six-frame" in result.output


def test_translate_from_string(runner):
    result = runner.invoke(main, ["translate", "--string", ">s1\nATGCATGCATGC"])
    assert result.exit_code == 0


def test_translate_no_input_raises(runner):
    result = runner.invoke(main, ["translate"])
    assert result.exit_code != 0


def test_translate_both_inputs_raises(runner):
    result = runner.invoke(
        main, ["translate", "--file", TINY_FASTA, "--string", ">s\nATGC"]
    )
    assert result.exit_code != 0


def test_translate_export_csv(runner, tmp_path):
    output = str(tmp_path / "proteins.csv")
    result = runner.invoke(
        main,
        ["translate", "--file", TINY_FASTA, "--output", output, "--format", "csv"],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_translate_export_json(runner, tmp_path):
    output = str(tmp_path / "proteins.json")
    result = runner.invoke(
        main,
        ["translate", "--file", TINY_FASTA, "--output", output, "--format", "json"],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_translate_export_fasta(runner, tmp_path):
    output = str(tmp_path / "proteins.fasta")
    result = runner.invoke(
        main,
        ["translate", "--file", TINY_FASTA, "--output", output, "--format", "fasta"],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_translate_six_frames_export_fasta(runner, tmp_path):
    output = str(tmp_path / "six_frame.fasta")
    result = runner.invoke(
        main,
        [
            "translate",
            "--file",
            TINY_FASTA,
            "--six-frames",
            "--output",
            output,
            "--format",
            "fasta",
        ],
    )
    assert result.exit_code == 0
    assert Path(output).exists()


def test_orf_basic(runner):
    result = runner.invoke(main, ["orf", "--file", TINY_FASTA])
    assert result.exit_code == 0
    assert "ORFs found" in result.output


def test_orf_min_length(runner):
    result = runner.invoke(main, ["orf", "--file", TINY_FASTA, "--min-length", "100"])
    assert result.exit_code == 0


def test_orf_overlap(runner):
    result = runner.invoke(main, ["orf", "--file", TINY_FASTA, "--overlap"])
    assert result.exit_code == 0


def test_orf_from_string(runner):
    result = runner.invoke(main, ["orf", "--string", ">s1\nATGAAACCCTTTGGGTAACCC"])
    assert result.exit_code == 0


def test_orf_no_input_raises(runner):
    result = runner.invoke(main, ["orf"])
    assert result.exit_code != 0


def test_orf_both_inputs_raises(runner):
    result = runner.invoke(main, ["orf", "--file", TINY_FASTA, "--string", ">s\nATGC"])
    assert result.exit_code != 0


def test_orf_export_csv(runner, tmp_path):
    output = str(tmp_path / "orfs.csv")
    result = runner.invoke(
        main,
        ["orf", "--file", TINY_FASTA, "--output", output, "--format", "csv"],
    )
    assert result.exit_code == 0


def test_orf_export_json(runner, tmp_path):
    output = str(tmp_path / "orfs.json")
    result = runner.invoke(
        main,
        ["orf", "--file", TINY_FASTA, "--output", output, "--format", "json"],
    )
    assert result.exit_code == 0


def test_orf_export_fasta(runner, tmp_path):
    output = str(tmp_path / "orfs.fasta")
    result = runner.invoke(
        main,
        ["orf", "--file", TINY_FASTA, "--output", output, "--format", "fasta"],
    )
    assert result.exit_code == 0


def test_motif_single_strand(runner):
    result = runner.invoke(main, ["motif", "--file", TINY_FASTA, "--pattern", "ATG"])
    assert result.exit_code == 0
    assert "Motifs" in result.output


def test_motif_both_strands(runner):
    result = runner.invoke(
        main,
        ["motif", "--file", TINY_FASTA, "--pattern", "ATG", "--mode", "both"],
    )
    assert result.exit_code == 0


def test_motif_search_all(runner):
    result = runner.invoke(
        main,
        [
            "motif",
            "--file",
            TINY_FASTA,
            "--pattern",
            "GGG",
            "--mode",
            "search-all",
        ],
    )
    assert result.exit_code == 0


def test_motif_with_mismatch(runner):
    result = runner.invoke(
        main,
        ["motif", "--file", TINY_FASTA, "--pattern", "ATGC", "--mismatch", "1"],
    )
    assert result.exit_code == 0


def test_motif_custom_k(runner):
    result = runner.invoke(
        main,
        ["motif", "--file", TINY_FASTA, "--pattern", "ATG", "--k", "3"],
    )
    assert result.exit_code == 0


def test_motif_from_string(runner):
    result = runner.invoke(
        main, ["motif", "--string", ">s1\nATGCATGCATGC", "--pattern", "ATG"]
    )
    assert result.exit_code == 0


def test_motif_no_input_raises(runner):
    result = runner.invoke(main, ["motif", "--pattern", "ATG"])
    assert result.exit_code != 0


def test_motif_both_inputs_raises(runner):
    result = runner.invoke(
        main,
        [
            "motif",
            "--file",
            TINY_FASTA,
            "--string",
            ">s\nATGC",
            "--pattern",
            "ATG",
        ],
    )
    assert result.exit_code != 0


def test_motif_export_csv(runner, tmp_path):
    output = str(tmp_path / "motifs.csv")
    result = runner.invoke(
        main,
        [
            "motif",
            "--file",
            TINY_FASTA,
            "--pattern",
            "ATG",
            "--output",
            output,
            "--format",
            "csv",
        ],
    )
    assert result.exit_code == 0


def test_motif_export_json(runner, tmp_path):
    output = str(tmp_path / "motifs.json")
    result = runner.invoke(
        main,
        [
            "motif",
            "--file",
            TINY_FASTA,
            "--pattern",
            "GGG",
            "--output",
            output,
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0


def test_config_show(runner):
    result = runner.invoke(main, ["config", "--show"])
    assert result.exit_code == 0
    assert "parsing" in result.output


def test_config_get_existing_key(runner):
    result = runner.invoke(main, ["config", "--get", "export.default_format"])
    assert result.exit_code == 0
    assert "csv" in result.output


def test_config_get_nonexistent_key(runner):
    result = runner.invoke(main, ["config", "--get", "fake.key"])
    assert result.exit_code != 0


def test_config_set_value(runner):
    result = runner.invoke(main, ["config", "--set", "motif.default_k", "10"])
    assert result.exit_code == 0
    assert "Set motif.default_k = 10" in result.output


def test_config_set_boolean_true(runner):
    result = runner.invoke(main, ["config", "--set", "parsing.strict_mode", "true"])
    assert result.exit_code == 0
    assert "True" in result.output


def test_config_set_boolean_false(runner):
    result = runner.invoke(main, ["config", "--set", "parsing.strict_mode", "false"])
    assert result.exit_code == 0
    assert "False" in result.output


def test_config_reset(runner):
    result = runner.invoke(main, ["config", "--reset"])
    assert result.exit_code == 0
    assert "reset" in result.output.lower()


def test_config_no_flags(runner):
    result = runner.invoke(main, ["config"])
    assert result.exit_code == 0
    assert "Use" in result.output


def test_print_sequence_lengths(capsys):
    seqs = [sequence("s1", "ATGCATGC"), sequence("s2", "GGG")]
    print_sequence_lengths_formatted(seqs)
    captured = capsys.readouterr()
    assert "s1" in captured.out
    assert "8" in captured.out
    assert "s2" in captured.out
    assert "3" in captured.out


def test_print_gc_content(capsys):
    seqs = [sequence("s1", "GGCC"), sequence("s2", "AAAA")]
    print_gc_content_table(seqs)
    captured = capsys.readouterr()
    assert "s1" in captured.out
    assert "100" in captured.out
    assert "s2" in captured.out
    assert "0" in captured.out


def test_print_revcomp(capsys):
    seqs = [sequence("s1", "ATGC")]
    print_revcomp(seqs)
    captured = capsys.readouterr()
    assert "reverse complement" in captured.out.lower()
    assert "GCAT" in captured.out


def test_print_base_count(capsys):
    seqs = [sequence("s1", "AAACCCGGGTTT")]
    print_base_count(seqs)
    captured = capsys.readouterr()
    assert "s1" in captured.out
    assert "3" in captured.out


def test_print_summary(capsys):
    seqs = [sequence("s1", "ATGCATGC")]
    print_summary(seqs)
    captured = capsys.readouterr()
    assert "SEQUENCE LENGTHS" in captured.out
    assert "GC CONTENT" in captured.out
    assert "BASE COMPOSITION" in captured.out
