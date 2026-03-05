from tabulate import tabulate
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector
from bio_seq_v1.motif_search import MotifFinder
from bio_seq_v1.export import Exporter
from bio_seq_v1.config import Config
import rich_click as click
import yaml
from pathlib import Path

click.rich_click.USE_RICH_MARKUP = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_COMMANDS_TABLE_COLUMN_WIDTH_RATIO = (1, 2)


@click.group()
@click.pass_context
def main(ctx):
    """BioSeq — A bioinformatics sequence analysis toolkit."""
    ctx.obj = Config()


@main.command()
@click.option(
    "--init", "do_init", is_flag=True, help="Interactively set up your config"
)
@click.option("--show", is_flag=True, help="Show current configuration")
@click.option(
    "--get", "get_key", default=None, help="Get a config value (e.g., motif.default_k)"
)
@click.option(
    "--set",
    "set_key",
    nargs=2,
    default=None,
    help="Set a config value (e.g., --set motif.default_k 6)",
)
@click.option("--reset", is_flag=True, help="Reset config to defaults")
@click.pass_context
def config(ctx, do_init, show, get_key, set_key, reset):
    """View or modify configuration settings."""
    cfg = ctx.obj
    try:
        if do_init:
            cfg.create_config(None)
        elif show:
            click.echo(yaml.dump(cfg.config, default_flow_style=False))
        elif get_key:
            value = cfg.get_config(get_key)
            if value is None:
                raise click.ClickException(f"Key '{get_key}' not found")
            click.echo(f"{get_key} = {value}")
        elif set_key:
            key, value = set_key
            if value.lower() == "true":
                value = True
            elif value.lower() == "false":
                value = False
            elif value.isdigit():
                value = int(value)
            cfg.set_config(key, value)
            output_path = Path.home() / ".bio_seq" / "config.yaml"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                yaml.dump(cfg.config, f, default_flow_style=False)
            click.echo(f"Set {key} = {value}")
            click.echo(f"Saved to {output_path}")
        elif reset:
            cfg.config = cfg._get_defaults()
            click.echo("Config reset to defaults")
        else:
            click.echo("Use --show, --get, --set, --init, or --reset")
    except ValueError as e:
        raise click.ClickException(str(e))


def print_sequence_lengths_formatted(sequences):
    """
    Print a formatted table of sequence IDs and their lengths.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    table = [[s.id, s.sequence_length()] for s in sequences]
    print(tabulate(table, headers=["Sequence ID", "Length"], tablefmt="grid"))


def print_gc_content_table(sequences):
    """
    Print a formatted table of sequence IDs and their GC content percentages.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    table = [[s.id, f"{s.gc_content():.2f}%"] for s in sequences]
    print(tabulate(table, headers=["Sequence", "GC%"], tablefmt="grid"))


def print_revcomp(sequences):
    """
    Print the reverse complement of each sequence in the list.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    for s in sequences:
        print(f">{s.id} reverse complement")
        print(s.rev_complement())
        print("-" * 30)


def print_base_count(sequences):
    """
    Print a table of the counts of each base for each sequence.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    all_counts = [s.base_count() for s in sequences]
    bases_present = [b for b in sequence.valid]
    table = []
    for counts, s in zip(all_counts, sequences):
        row = [s.id] + [counts.get(b, 0) for b in bases_present]
        table.append(row)
    headers = ["Sequence"] + bases_present
    print(tabulate(table, headers=headers, tablefmt="grid"))


def print_summary(sequences):
    """
    Print a full summary of sequences including lengths, GC content, and base composition.

    Args:
        sequences (list of sequence): List of Sequence objects to process.
    """
    print("SEQUENCE LENGTHS")
    print_sequence_lengths_formatted(sequences)
    print()

    print("GC CONTENT")
    print_gc_content_table(sequences)
    print()

    print("BASE COMPOSITION")
    print_base_count(sequences)


@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option("--length", "-l", is_flag=True, help="Compute sequence lengths")
@click.option("--gc", is_flag=True, help="Compute GC content")
@click.option("--revcomp", "-rc", is_flag=True, help="Compute reverse complements")
@click.option("--basecount", "-b", is_flag=True, help="Compute base counts")
@click.option("--summary", is_flag=True, help="Print full summary")
@click.option(
    "--format",
    "export_format",
    type=click.Choice(["csv", "tsv", "json"]),
    default="csv",
    help="Export format: to csv, to tsv, to json, sequences to csv",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (prints to stdout if not given)",
)
def stats(
    file,
    string,
    strict,
    strict_file,
    strict_seq,
    length,
    gc,
    revcomp,
    basecount,
    summary,
    output,
    export_format,
):
    """Analyze sequences — lengths, GC content, base composition."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")
    fasta_parser = FASTAParser(
        path=file, strict=strict, strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err=True)
        return
    click.echo(fasta_parser.get_report())
    click.echo()
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences parsed.")
    if not any([length, gc, revcomp, basecount, summary]):
        summary = True
    export_data = [{"id": s.id} for s in sequences] if output else None
    try:
        if length or summary:
            print_sequence_lengths_formatted(sequences)
            click.echo()
            if export_data:
                for i, s in enumerate(sequences):
                    export_data[i]["length"] = s.sequence_length()
        if gc or summary:
            print_gc_content_table(sequences)
            click.echo()
            if export_data:
                for i, s in enumerate(sequences):
                    export_data[i]["gc_content"] = s.gc_content()
        if revcomp or summary:
            print_revcomp(sequences)
            click.echo()
            if export_data:
                for i, s in enumerate(sequences):
                    export_data[i]["rev_comp"] = s.rev_complement()
        if basecount or summary:
            print_base_count(sequences)
            click.echo()
            if export_data:
                for i, s in enumerate(sequences):
                    export_data[i]["base_count"] = s.base_count()
        if output and export_data:
            if export_format == "csv":
                Exporter.to_csv(export_data, file_path=output)
            elif export_format == "tsv":
                Exporter.to_tsv(export_data, file_path=output)
            elif export_format == "json":
                Exporter.to_json(export_data, file_path=output)
            click.echo(f"Results saved to {output}")
    except ValueError as e:
        raise click.ClickException(str(e))


@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option(
    "--frame",
    type=click.IntRange(0, 2),
    default=0,
    help="Frame for translation (0, 1, or 2)",
)
@click.option(
    "--six-frames",
    is_flag=True,
    help="Translate the parsed sequences in all six frames",
)
@click.option(
    "--format",
    "export_format",
    type=click.Choice(["csv", "tsv", "json", "fasta"]),
    default="csv",
    help="Export format: to csv, to tsv, to json, to fasta",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (prints to stdout if not given)",
)
def translate(
    file,
    string,
    strict,
    strict_file,
    strict_seq,
    frame,
    six_frames,
    export_format,
    output,
):
    """Translate given sequence."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")

    fasta_parser = FASTAParser(
        path=file, strict=strict, strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err=True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")
    translator = Translator()
    export_data = [] if output else None
    try:
        for seq in sequences:
            if six_frames:
                results = translator.translate_six_frames(seq)
                click.echo(f">{seq.id} — Six-frame translation")
                row = {"id": seq.id}
                for frame_label, protein in results.items():
                    click.echo(f" Frame {frame_label}: {protein}")
                    row[f"frame_{frame_label}"] = protein
                if export_data is not None:
                    export_data.append(row)
            else:
                protein = translator.translate(seq, frame)
                click.echo(f">{seq.id} — Frame {frame}")
                click.echo(f"  {protein}")
                if export_data is not None:
                    export_data.append(
                        {"id": seq.id, "frame": frame, "protein": protein}
                    )
            click.echo()
        if output and export_data:
            if export_format == "csv":
                Exporter.to_csv(export_data, file_path=output)
            elif export_format == "tsv":
                Exporter.to_tsv(export_data, file_path=output)
            elif export_format == "json":
                Exporter.to_json(export_data, file_path=output)
            elif export_format == "fasta":
                fasta_seqs = []
                for row in export_data:
                    if "protein" in row:
                        fasta_seqs.append(sequence(row["id"], row["protein"]))
                    else:
                        for key, val in row.items():
                            if key.startswith("frame_"):
                                fasta_seqs.append(sequence(f"{row['id']}_{key}", val))
                Exporter.to_fasta(fasta_seqs, file_path=output)
            click.echo(f"Results saved to {output}")
    except ValueError as e:
        raise click.ClickException(str(e))


@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option("--min-length", default=0, help="Mininmum length for open reading frames")
@click.option("--overlap", is_flag=True, help="Overlapping ORFs")
@click.option(
    "--format",
    "export_format",
    type=click.Choice(["csv", "tsv", "json", "fasta"]),
    default="csv",
    help="Export format: to csv, to tsv, to json, to fasta",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (prints to stdout if not given)",
)
def orf(
    file,
    string,
    strict,
    strict_file,
    strict_seq,
    min_length,
    overlap,
    export_format,
    output,
):
    """Find open reading frames in sequences."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")
    fasta_parser = FASTAParser(
        path=file, strict=strict, strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err=True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")
    detector = ORFDetector(min_length=min_length)
    all_orfs = []
    try:
        for seq in sequences:
            results = detector.find_orfs(seq)
            click.echo(f">{seq.id} — ORFs found: {len(results)}")
            for found_orf in results:
                click.echo(
                    f"  start={found_orf.start}, end={found_orf.end}, protein={found_orf.protein}"
                )
            all_orfs.extend(results)
            if overlap:
                overlaps = detector.overlapping_orfs(results)
                if overlaps:
                    click.echo(f"  Overlapping pairs: {len(overlaps)}")
                    for orf1, orf2 in overlaps:
                        click.echo(
                            f"    {orf1.start}-{orf1.end} ↔ {orf2.start}-{orf2.end}"
                        )
            click.echo()
        if output and all_orfs:
            if export_format == "csv":
                Exporter.orfs_to_csv(all_orfs, file_path=output)
            elif export_format == "tsv":
                data = [o.to_dict() for o in all_orfs]
                Exporter.to_tsv(data, file_path=output)
            elif export_format == "json":
                data = [o.to_dict() for o in all_orfs]
                Exporter.to_json(data, file_path=output)
            elif export_format == "fasta":
                fasta_seqs = [
                    sequence(f"{o.seq_id}_orf_{o.start}_{o.end}", o.protein)
                    for o in all_orfs
                ]
                Exporter.to_fasta(fasta_seqs, file_path=output)
            click.echo(f"Results saved to {output}")
    except ValueError as e:
        raise click.ClickException(str(e))


@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option(
    "--mode",
    type=click.Choice(["single", "both", "search-all"]),
    default="single",
    help="Search mode: single strand, both strands, or all sequences",
)
@click.option("--k", default=3, help="Mininmum length for motif")
@click.option("--mismatch", "-m", default=0, help="Numer of mismatches allowed")
@click.option(
    "--pattern", "-p", required=True, help="Motif pattern to search for (e.g., TATAAA)"
)
@click.option(
    "--format",
    "export_format",
    type=click.Choice(["csv", "tsv", "json"]),
    default="csv",
    help="Export format: to csv, to tsv, to json",
)
@click.option(
    "--output",
    "-o",
    default=None,
    help="Output file path (prints to stdout if not given)",
)
def motif(
    file,
    string,
    strict,
    strict_file,
    strict_seq,
    mode,
    k,
    mismatch,
    pattern,
    export_format,
    output,
):
    """Find motifs in sequences."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")

    fasta_parser = FASTAParser(
        path=file, strict=strict, strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err=True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")

    finder = MotifFinder(k=k)
    all_motifs = []
    try:
        if mode == "single":
            for seq in sequences:
                single_strand = finder.search_single(seq, pattern, mismatch)
                click.echo(
                    f">{seq.id} — Motifs on single strand found: {len(single_strand)}"
                )
                for found_motif in single_strand:
                    click.echo(
                        f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}"
                    )
                all_motifs.extend(single_strand)
        if mode == "both":
            for seq in sequences:
                double_strand = finder.search_both_strands(seq, pattern, mismatch)
                click.echo(
                    f">{seq.id} — Motifs on both strands found: {len(double_strand)}"
                )
                for found_motif in double_strand:
                    click.echo(
                        f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}"
                    )
                all_motifs.extend(double_strand)
        if mode == "search-all":
            fasta_motifs = finder.search_fasta(sequences, pattern, mismatch)
            click.echo(f"Motifs found across all sequences: {len(fasta_motifs)}")
            for found_motif in fasta_motifs:
                click.echo(
                    f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}"
                )
            all_motifs.extend(fasta_motifs)
        if output and all_motifs:
            if export_format == "csv":
                Exporter.motifs_to_csv(all_motifs, file_path=output)
            elif export_format == "tsv":
                data = [o.to_dict() for o in all_motifs]
                Exporter.to_tsv(data, file_path=output)
            elif export_format == "json":
                data = [o.to_dict() for o in all_motifs]
                Exporter.to_json(data, file_path=output)
    except ValueError as e:
        raise click.ClickException(str(e))
