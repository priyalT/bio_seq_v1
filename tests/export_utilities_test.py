import pytest
import csv
import io
import json
import sys
import os
import tempfile
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.export import Exporter
from hypothesis import given, strategies as st
from bio_seq_v1.stats import sequence


@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size=1))
def test_csv_output(fasta):
    seq = sequence("id", fasta)
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv") as tmp_file:
        tmp_filepath = tmp_file.name
    try:
        Exporter.sequences_to_csv([seq], tmp_filepath)
        with open(tmp_filepath, newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            loaded_data = list(reader)

        assert loaded_data[0]["id"] == seq.id
        assert loaded_data[0]["sequence"] == seq.sequence

    finally:
        os.unlink(tmp_filepath)


@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size=1))
def test_json_output(fasta):
    seq = sequence("id", fasta)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".json"
    ) as tmp_file:
        tmp_filepath = tmp_file.name
    try:
        Exporter.to_json([seq.to_dict()], tmp_filepath)
        with open(tmp_filepath) as jsonfile:
            loaded_data = json.load(jsonfile)

        assert loaded_data[0]["id"] == seq.id
        assert loaded_data[0]["sequence"] == seq.sequence

    finally:
        os.unlink(tmp_filepath)


@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size=1))
def test_tsv_output(fasta):
    seq = sequence("id", fasta)
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tsv") as tmp_file:
        tmp_filepath = tmp_file.name
    try:
        Exporter.to_tsv([seq.to_dict()], tmp_filepath)
        with open(tmp_filepath, newline="") as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            loaded_data = list(reader)

        assert loaded_data[0]["id"] == seq.id
        assert loaded_data[0]["sequence"] == seq.sequence

    finally:
        os.unlink(tmp_filepath)


@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size=1))
def test_fasta_output(fasta):
    seq = sequence("id", fasta)
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".fasta"
    ) as tmp_file:
        tmp_filepath = tmp_file.name
    try:
        Exporter.to_fasta([seq], tmp_filepath)
        parser = FASTAParser(path=tmp_filepath)
        parser.parse_file()
        loaded_data = parser.sequences

        assert loaded_data[0].id == seq.id
        assert loaded_data[0].sequence == seq.sequence

    finally:
        os.unlink(tmp_filepath)
