from tabulate import tabulate
import argparse
from fasta import fasta_parser
from stats import print_sequence_lengths, base_count, gc_content, rev_complement
valid = "ACGTUNRYSWKMBDHV-."

def print_sequence_lengths_formatted(sequences):
    table = [[s['id'], len(s['sequence'])] for s in sequences]
    print(tabulate(table, headers=["Sequence ID", "Length"], tablefmt="grid"))

def print_gc_content_table(sequences):
    gc_list = gc_content(sequences)
    table = [[item["id"], f"{item['GC%']:.2f}%"] for item in gc_list]
    print(tabulate(table, headers=["Sequence", "GC%"], tablefmt="grid"))

def print_revcomp(sequences):
    rev_list = rev_complement(sequences)
    for item in rev_list:
        print(f">{item['id']}_revcomp")
        print(item["reverse complement"])
        print("-" * 30)

def print_base_count(sequences):
    base_list = base_count(sequences)
    bases_present = [b for b in valid if any(b in item for item in base_list)]
    table = []
    for b in base_list:
        row = [b["id"]] + [b.get(base, 0) for base in bases_present]
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

sequences = fasta_parser(args.file)
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

