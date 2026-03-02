from tabulate import tabulate
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORF
from bio_seq_v1.orf import ORFDetector
from bio_seq_v1.motif_search import MotifFinder
from bio_seq_v1.export import Exporter
import rich_click as click
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_COMMANDS_TABLE_COLUMN_WIDTH_RATIO = (1, 2)


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


@click.group()
def main():
    """
    Command-line interface entry point.

    Parses arguments to specify FASTA input source and analysis options,
    and prints the requested sequence information.
    """
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
@click.option("--format", "export_format", type=click.Choice(["csv", "tsv", "json"]), default="csv",
              help="Export format: to csv, to tsv, to json, sequences to csv")
@click.option("--output", "-o", default=None, help="Output file path (prints to stdout if not given)")
def stats(file, string, strict, strict_file, strict_seq, length, gc, revcomp, basecount, summary, output, export_format):
    """Analyze sequences — lengths, GC content, base composition."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")
    fasta_parser = FASTAParser(
        path=file,
        strict=strict,
        strict_file=strict_file,
        strict_seq=strict_seq
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



@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option("--frame", type = click.IntRange(0, 2), default=0, help="Frame for translation (0, 1, or 2)")
@click.option("--six-frames", is_flag=True, help="Translate the parsed sequences in all six frames")
@click.option("--format", "export_format", type=click.Choice(["csv", "tsv", "json"]), default="csv",
              help="Export format: to csv, to tsv, to json, sequences to csv")
@click.option("--output", "-o", default=None, help="Output file path (prints to stdout if not given)")
def translate(file, string, strict, strict_file, strict_seq, frame, six_frames, export_format, output):
    """Translate given sequence."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")

    fasta_parser = FASTAParser(
        path=file, strict=strict,
        strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err= True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")
    translator = Translator()
    export_data = [] if output else None


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
                export_data.append({"id": seq.id, "frame": frame, "protein": protein})
        click.echo()
    if output and export_data:
        if export_format == "csv":
            Exporter.to_csv(export_data, file_path=output)
        elif export_format == "tsv":
            Exporter.to_tsv(export_data, file_path=output)
        elif export_format == "json":
            Exporter.to_json(export_data, file_path=output)
        click.echo(f"Results saved to {output}")



@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option("--min-length", default=0, help="Mininmum length for open reading frames")
@click.option("--overlap", is_flag = True, help="Overlapping ORFs")
@click.option("--format", "export_format", type=click.Choice(["orfs-to-csv", "tsv", "json"]), default="csv",
              help="Export format: to csv, to tsv, to json, sequences to csv")
@click.option("--output", "-o", default=None, help="Output file path (prints to stdout if not given)")
def orf(file, string, strict, strict_file, strict_seq, min_length, overlap, export_format, output):
    """Find open reading frames in sequences."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")

    fasta_parser = FASTAParser(
        path=file, strict=strict,
        strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err= True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")
    detector = ORFDetector(min_length=min_length)

    for seq in sequences:
        results = detector.find_orfs(seq)      
        click.echo(f">{seq.id} — ORFs found: {len(results)}")
        for found_orf in results:
            click.echo(f"  start={found_orf.start}, end={found_orf.end}, protein={found_orf.protein}")
        if overlap: 
            overlaps = detector.overlapping_orfs(results)
            if overlaps:
                click.echo(f"  Overlapping pairs: {len(overlaps)}")
                for orf1, orf2 in overlaps:
                    click.echo(f"    {orf1.start}-{orf1.end} ↔ {orf2.start}-{orf2.end}")
        click.echo()

@main.command()
@click.option("--file", "-f", default=None, help="Path to the FASTA file")
@click.option("--string", "-s", default=None, help="FASTA-formatted string")
@click.option("--strict", is_flag=True, help="Enable strict parsing")
@click.option("--strict-file", is_flag=True, help="Enable strict file validation")
@click.option("--strict-seq", is_flag=True, help="Fail on invalid sequence characters")
@click.option("--mode", type=click.Choice(["single", "both", "search-all"]), default="single",
              help="Search mode: single strand, both strands, or all sequences")
@click.option("--k", default=3, help="Mininmum length for motif")
@click.option("--mismatch", "-m", default=0, help="Numer of mismatches allowed")
@click.option("--pattern", "-p", required=True, help="Motif pattern to search for (e.g., TATAAA)")
def motif(file, string, strict, strict_file, strict_seq, mode, k, mismatch, pattern):
    """Find motifs in sequences."""

    if not file and not string:
        raise click.UsageError("Must provide either --file or --string")
    if file and string:
        raise click.UsageError("Cannot use both --file and --string")

    fasta_parser = FASTAParser(
        path=file, strict=strict,
        strict_file=strict_file, strict_seq=strict_seq
    )
    try:
        if file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(string)
    except Exception as e:
        click.echo(f"Parsing failed: {e}", err= True)
        return
    sequences = fasta_parser.sequences
    if not sequences:
        raise click.ClickException("No valid sequences.")
    
    finder = MotifFinder(k=k)

    if mode == 'single':
        for seq in sequences:
            single_strand = finder.search_single(seq, pattern, mismatch)
            click.echo(f">{seq.id} — Motifs on single strand found: {len(single_strand)}")
            for found_motif in single_strand: 
                click.echo(f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}")
    if mode == 'both':
        for seq in sequences:
            double_strand = finder.search_both_strands(seq, pattern, mismatch)    
            click.echo(f">{seq.id} — Motifs on both strands found: {len(double_strand)}")
            for found_motif in double_strand: 
                click.echo(f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}")
    if mode == 'search-all':
        fasta_motifs = finder.search_fasta(sequences, pattern, mismatch)
        click.echo(f"Motifs found across all sequences: {len(fasta_motifs)}")
        for found_motif in fasta_motifs: 
            click.echo(f"  position = {found_motif.position}, matched sequence = {found_motif.matched_seq}, strand attributes = {found_motif.strand_attributes}")

