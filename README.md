# 🧬 BioSeq

<img width="500" height="500" alt="Logo of bioseq" src="https://github.com/user-attachments/assets/1ac4ea8a-bdd7-4768-8b86-8b43b46f5016" />

**A Python CLI toolkit for bioinformatics sequence analysis.**

![Tests](https://github.com/priyalT/bio_seq-v1-/actions/workflows/test.yml/badge.svg)
![Lint](https://github.com/priyalT/bio_seq-v1-/actions/workflows/lint.yml/badge.svg)
![Docker](https://github.com/priyalT/bio_seq-v1-/actions/workflows/docker.yml/badge.svg)
[![codecov](https://codecov.io/gh/priyalT/bio_seq-v1-/branch/main/graph/badge.svg)](https://codecov.io/gh/priyalT/bio_seq-v1-)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

---

## What is BioSeq?

BioSeq is a command-line toolkit for analyzing biological sequences. It parses FASTA files and provides tools for sequence statistics, translation, ORF detection, and motif searching — all from your terminal.

### Features

- 📊 **Sequence Statistics** — lengths, GC content, base composition, reverse complements
- 🧪 **Translation** — DNA to protein in any reading frame, or all six frames
- 🔍 **ORF Detection** — find open reading frames with configurable minimum length and overlap detection
- 🎯 **Motif Search** — exact and fuzzy pattern matching on single/both strands
- 📁 **Export** — results to CSV, TSV, JSON, or FASTA format
- ⚙️ **Configuration** — YAML-based config system with sensible defaults
- 🐳 **Docker** — containerized and available on Quay.io
- ✅ **Tested** — property-based tests with Hypothesis + pytest

---

## Installation

### Option 1: pip (from source)

```bash
git clone https://github.com/priyalT/bio_seq-v1-.git
cd bio_seq_v1
pip install -e .
```

### Option 2: Docker

```bash
docker pull quay.io/priyal_tripathi/bioseq
docker run quay.io/priyal_tripathi/bioseq --help
```

### Option 3: pip install only

```bash
pip install .
```

### Requirements

- Python 3.9+
- Dependencies (installed automatically): `tabulate`, `pyyaml`, `rich-click`

---

## Quick Start

```bash
# View all available commands
bioseq --help

# Analyze a FASTA file (prints full summary by default)
bioseq stats --file sequences.fasta

# Translate sequences to protein
bioseq translate --file sequences.fasta

# Find open reading frames
bioseq orf --file sequences.fasta

# Search for a motif pattern
bioseq motif --file sequences.fasta --pattern TATAAA
```

---

## Usage

### Input Format

BioSeq accepts standard FASTA files or inline FASTA strings:

```fasta
>seq1 Example sequence
ATGCGTACGTAGCTAGTTAGCGATCG
GGGCTAGCTAGCTAGCTAG

>seq2 Another sequence
GGGTTTAAACCCGGGCCCGGGAAATTT
```

You can provide input via file or string:

```bash
# From a file
bioseq stats --file sequences.fasta

# From a string
bioseq stats --string ">seq1\nATGCGTACG"
```

---

### `bioseq stats` — Sequence Statistics

Analyze sequences for lengths, GC content, base composition, and reverse complements.

```bash
# Full summary (default)
bioseq stats --file sequences.fasta

# Specific analyses
bioseq stats --file sequences.fasta --length
bioseq stats --file sequences.fasta --gc
bioseq stats --file sequences.fasta --revcomp
bioseq stats --file sequences.fasta --basecount

# Export results
bioseq stats --file sequences.fasta --output results.json --format json
```

| Flag | Short | Description |
|------|-------|-------------|
| `--file` | `-f` | Path to FASTA file |
| `--string` | `-s` | FASTA-formatted string |
| `--length` | `-l` | Compute sequence lengths |
| `--gc` | | Compute GC content |
| `--revcomp` | `-rc` | Compute reverse complements |
| `--basecount` | `-b` | Compute base composition |
| `--summary` | | Print all statistics (default) |
| `--strict` | | Enable strict parsing |
| `--format` | | Export format: `csv`, `tsv`, `json` |
| `--output` | `-o` | Output file path |

---

### `bioseq translate` — DNA to Protein Translation

Translate DNA sequences to protein in a specific reading frame or all six frames.

```bash
# Default frame (0)
bioseq translate --file sequences.fasta

# Specific frame
bioseq translate --file sequences.fasta --frame 2

# All six frames
bioseq translate --file sequences.fasta --six-frames

# Export as FASTA
bioseq translate --file sequences.fasta --six-frames --output proteins.fasta --format fasta
```

| Flag | Description |
|------|-------------|
| `--frame` | Reading frame: `0`, `1`, or `2` (default: 0) |
| `--six-frames` | Translate in all six reading frames |
| `--format` | Export format: `csv`, `tsv`, `json`, `fasta` |
| `--output` | Output file path |

---

### `bioseq orf` — Open Reading Frame Detection

Find ORFs in sequences with configurable minimum length and overlap detection.

```bash
# Find all ORFs
bioseq orf --file sequences.fasta

# Minimum length filter
bioseq orf --file sequences.fasta --min-length 100

# Show overlapping ORFs
bioseq orf --file sequences.fasta --overlap

# Export to JSON
bioseq orf --file sequences.fasta --output orfs.json --format json
```

| Flag | Description |
|------|-------------|
| `--min-length` | Minimum ORF length (default: 0) |
| `--overlap` | Show overlapping ORF pairs |
| `--format` | Export format: `csv`, `tsv`, `json`, `fasta` |
| `--output` | Output file path |

---

### `bioseq motif` — Motif Search

Search for sequence motifs with exact or fuzzy matching.

```bash
# Search single strand
bioseq motif --file sequences.fasta --pattern TATAAA

# Search both strands
bioseq motif --file sequences.fasta --pattern TATAAA --mode both

# Fuzzy matching (allow 1 mismatch)
bioseq motif --file sequences.fasta --pattern TATAAA --mismatch 1

# Search across all sequences
bioseq motif --file sequences.fasta --pattern TATAAA --mode search-all
```

| Flag | Short | Description |
|------|-------|-------------|
| `--pattern` | `-p` | Motif pattern to search for (required) |
| `--mode` | | `single`, `both`, or `search-all` (default: single) |
| `--mismatch` | `-m` | Number of mismatches allowed (default: 0) |
| `--k` | | Minimum motif length (default: 3) |
| `--format` | | Export format: `csv`, `tsv`, `json` |
| `--output` | | Output file path |

---

### `bioseq config` — Configuration

Manage BioSeq configuration with a YAML-based config system.

```bash
# Initialize config interactively
bioseq config --init

# Show current config
bioseq config --show

# Get a specific value
bioseq config --get motif.default_k

# Set a value
bioseq config --set motif.default_k 6

# Reset to defaults
bioseq config --reset
```

---

## Docker Usage

### Pull and run

```bash
docker pull quay.io/priyal_tripathi/bioseq
docker run quay.io/priyal_tripathi/bioseq --help
```

### Analyze local files using volume mounts

```bash
# Mount a local directory into the container
docker run -v $(pwd)/data:/data quay.io/priyal_tripathi/bioseq stats --file /data/sequences.fasta

# Export results to your local machine
docker run -v $(pwd)/data:/data -v $(pwd)/output:/output \
  quay.io/priyal_tripathi/bioseq stats --file /data/sequences.fasta --output /output/results.json --format json
```

---

## Supported Nucleotide Bases

BioSeq supports standard and IUPAC nucleotide codes:

`A`, `C`, `G`, `T`, `U`, `N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `-`, `.`

---

## Test Data

Sample FASTA files are included in [`tests/data/`](tests/data/) for testing:

| File | Description |
|------|-------------|
| `tiny.fasta` | 3 short sequences for quick testing |
| `single.fasta` | Single sequence for edge case testing |

---

## Project Structure

```
bio_seq_v1/
├── bio_seq_v1/              # Source package
│   ├── cli.py               # CLI commands (Click)
│   ├── config.py            # YAML configuration system
│   ├── export.py            # Export to CSV/TSV/JSON/FASTA
│   ├── fasta.py             # FASTA parser
│   ├── motif_search.py      # Motif pattern matching
│   ├── orf.py               # Open reading frame detection
│   ├── stats.py             # Sequence statistics
│   └── translator.py        # DNA to protein translation
├── tests/                   # Test suite
│   ├── data/                # Test FASTA files
│   └── *_test.py            # Property-based & unit tests
├── examples/                # Example scripts
├── .github/workflows/       # CI/CD pipelines
│   ├── test.yml             # Pytest across Python 3.9–3.11
│   ├── lint.yml             # Black + Ruff
│   └── docker.yml           # Build & push to Quay.io
├── Dockerfile               # Container image
├── pyproject.toml           # Package configuration
└── README.md
```

---

## Development

### Prerequisites

- Python 3.9+
- pip

### Setup

```bash
git clone https://github.com/priyalT/bio_seq-v1-.git
cd bio_seq_v1
pip install -e .
pip install pytest hypothesis pytest-cov black ruff
```

### Run Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=bio_seq_v1 --cov-branch --cov-report=html

# Open coverage report
open htmlcov/index.html
```

### Lint & Format

```bash
# Check formatting
black --check bio_seq_v1/ tests/

# Auto-format
black bio_seq_v1/ tests/

# Lint
ruff check bio_seq_v1/ tests/

# Auto-fix lint issues
ruff check bio_seq_v1/ tests/ --fix
```

### Build Docker Image

```bash
docker build -t bioseq .
docker run bioseq --help
```

---

## CI/CD

This project uses GitHub Actions for continuous integration:

| Workflow | What it does |
|----------|--------------|
| **Tests** | Runs pytest on Python 3.9, 3.10, 3.11 |
| **Lint** | Checks formatting (Black) and linting (Ruff) |
| **Docker** | Builds image, smoke tests, pushes to Quay.io |

All workflows run automatically on every push to `main` and on pull requests.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Author

**Priyal Tripathi** — [priyaltripathi2910@gmail.com](mailto:priyaltripathi2910@gmail.com)
