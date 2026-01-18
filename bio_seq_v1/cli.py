from tabulate import tabulate
import argparse
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence

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

def main():
    """
    Command-line interface entry point.

    Parses arguments to specify FASTA input source and analysis options,
    and prints the requested sequence information.
    """
    arg_parser = argparse.ArgumentParser()

    input_group = arg_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--file", "-f", help="Path to the FASTA file")
    input_group.add_argument("--string", "-s", help="FASTA-formatted string")

    arg_parser.add_argument("--strict", action="store_true",
                            help="Enable strict parsing (fail on structural errors)")
    arg_parser.add_argument("--strict-file", action="store_true",
                            help="Enable strict file validation")
    arg_parser.add_argument("--strict-seq", action="store_true",
                            help="Fail on invalid sequence characters")

    arg_parser.add_argument("--length", "-l", action="store_true",
                            help="Compute sequence length per sequence")
    arg_parser.add_argument("--gc", action="store_true",
                            help="Compute GC content per sequence")
    arg_parser.add_argument("--revcomp", "-rc", action="store_true",
                            help="Compute reverse complements per sequence")
    arg_parser.add_argument("--basecount", "-b", action="store_true",
                            help="Compute base counts per sequence")
    arg_parser.add_argument("--summary", action="store_true",
                            help="Print summary statistics")

    args = arg_parser.parse_args()

    fasta_parser = FASTAParser(
        path=args.file,
        strict=args.strict,
        strict_file=args.strict_file,
        strict_seq=args.strict_seq
    )

    try:
        if args.file:
            fasta_parser.parse_file()
        else:
            fasta_parser.parse_string(args.string)
    except Exception as e:
        print(f"Parsing failed: {e}")
        return

    print(fasta_parser.get_report())
    print()

    sequences = fasta_parser.sequences

    if not sequences:
        raise ValueError("No valid sequences parsed.")

    if not any([args.length, args.gc, args.revcomp, args.basecount, args.summary]):
        print_summary(sequences)
        return

    if args.length:
        print_sequence_lengths_formatted(sequences)
        print()

    if args.gc:
        print_gc_content_table(sequences)
        print()

    if args.revcomp:
        print_revcomp(sequences)
        print()

    if args.basecount:
        print_base_count(sequences)
        print()

    if args.summary:
        print_summary(sequences)