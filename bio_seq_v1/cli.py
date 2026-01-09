from tabulate import tabulate
import argparse
from fasta import fasta_parser
from stats import sequence
valid = "ACGTUNRYSWKMBDHV-."

def print_sequence_lengths_formatted(sequences):
    table = [[s.id, s.print_sequence_lengths()] for s in sequences]
    print(tabulate(table, headers=["Sequence ID", "Length"], tablefmt="grid"))

def print_gc_content_table(sequences):
    table = [[s.id, f"{s.gc_content():.2f}%"] for s in sequences]
    print(tabulate(table, headers=["Sequence", "GC%"], tablefmt="grid"))

def print_revcomp(sequences):
    for s in sequences:
        print(f">{s.id} reverse complement")
        print(s.rev_complement())
        print("-" * 30)

def print_base_count(sequences):
    all_counts = [s.base_count() for s in sequences]
    bases_present = [b for b in sequence.valid]
    table = []
    for counts, s in zip(all_counts, sequences):
        row = [s.id] + [counts.get(b, 0) for b in bases_present]
        table.append(row)
    headers = ["Sequence"] + bases_present
    print(tabulate(table, headers=headers, tablefmt="grid"))

def print_summary(sequences):
    print("SEQUENCE LENGTHS")
    print_sequence_lengths_formatted(sequences)
    print()

    print("GC CONTENT")
    print_gc_content_table(sequences)
    print()

    print("BASE COMPOSITION")
    print_base_count(sequences)


parser = argparse.ArgumentParser()

parser.add_argument("--file", "-f", help="Path to the FASTA file", required=True)
parser.add_argument("--length", "-l", help ="Compute sequence length per sequence", action="store_true")
parser.add_argument("--gc", help ="Compute GC% per sequence", action="store_true")
parser.add_argument("--revcomp", "-rc", help ="Compute reverse complements per sequence", action="store_true")
parser.add_argument("--basecount", "-b", help ="Compute total count for bases per sequence", action="store_true")
parser.add_argument("--summary", help="Print summary statistics", action="store_true")
args = parser.parse_args()

seq_dict = fasta_parser(args.file)
sequences = [sequence(d['id'], d['sequence']) for d in seq_dict]
if not any([args.length, args.gc, args.revcomp, args.basecount, args.summary]):
    print_summary(sequences)
    exit()
if not sequences:
    raise ValueError("No sequences found in FASTA file")

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
    exit()

