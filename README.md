# Bio_seq (version 1.0)

Hi! Welcome to my very first personal project repository! This is bio_seq, a CLI tool for sequence analysis. 


  <img width="500" height="500" alt="Logo of bioseq" src="https://github.com/user-attachments/assets/1ac4ea8a-bdd7-4768-8b86-8b43b46f5016" />

### What does it do?
- Parse FASTA files
- Compute sequence lengths
- Calculate GC content
- Generate reverse complements
- Count nucleotide/base composition (IUPAC supported)
- Command-line interface (CLI)

### How does one install it?
Clone the repository:

```bash
git clone https://github.com/yourusername/bio_seq_v1.git
cd bio_seq_
```
Install dependencies:
```bash
pip install tabulate
```

### Input format
The tool accepts standard FASTA files.
Example:

```fasta
>seq1
ATGCGTACGTAGCTAGTTAGCGATCG
GGGCTAGCTAGCTAGCTAG

>seq2
GGGTTTAAACCCGGGCCCGGGAAATTT
```
### How to use it?
Run the CLI using:
```bash
python cli.py -f test.fasta (your filename here)
```
---
Description of the flags:
| Flag | Description |
|------|-------------|
|-f, --file | Path to FASTA file (required)|
|-l, --length | Print sequence lengths|
|--gc |  Print GC content per sequence|
|-rc, --revcomp| Print reverse complements|
|-b, --basecount| Print base composition|
|--summary| Print all statistics (default)|

### Some example commands and outputs
```bash
python cli.py -f test.fasta --gc
```
<img width="210" height="134" alt="Screenshot 2026-01-09 at 5 30 32 AM" src="https://github.com/user-attachments/assets/b52de267-4820-4fca-b596-bf2d1fae6b8b" />

```bash
python cli.py -f test.fasta --revcomp
```
<img width="368" height="140" alt="Screenshot 2026-01-09 at 5 31 09 AM" src="https://github.com/user-attachments/assets/06952c26-9371-4bfa-83d9-07c69a776262" />

### Supported nucleotide bases
The tool supports standard and IUPAC nucleotide codes, i.e ACGTUNRYSWKMBDHV-.

### Project structure
```
bio_seq_v1/
├── src/
│   ├── cli.py
│   ├── fasta.py
│   └── stats.py
├── docs/
│   └── usage.md
├── tests/
│   ├── test_fasta_parser.py
│   └── test_stats.py
└── README.md
```
### Future Improvements
- Add translation to amino acids
- Add ORF detection
- Package as installable CLI tool
- Add unit tests