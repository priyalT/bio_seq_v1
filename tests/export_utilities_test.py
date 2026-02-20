import pytest
import csv
import io
import json
import sys
import os
import tempfile
from bio_seq_v1.stats import sequence
from bio_seq_v1.export import Exporter
from hypothesis import given, strategies as st
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector
from bio_seq_v1.fasta import FASTAParser

@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size= 1))
def test_output_file_writing(fasta):
    seq = sequence("id", fasta)
    with tempfile.NamedTemporaryFile(mode = 'w', delete = False, suffix='.csv') as tmp_file:
        tmp_filepath = tmp_file.name
    try:
        Exporter.sequences_to_csv([seq], tmp_filepath)
        with open(tmp_filepath, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            loaded_data = list(reader)
        
        assert loaded_data[0]["id"] == seq.id
        assert loaded_data[0]["sequence"] == seq.sequence

        
    finally:
        os.unlink(tmp_filepath)