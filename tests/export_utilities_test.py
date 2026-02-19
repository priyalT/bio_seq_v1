import pytest
import csv
import io
import json
import sys
import os
import tempfile
from bio_seq_v1.export import Exporter
from hypothesis import given, strategies as st
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector
from bio_seq_v1.fasta import FASTAParser

# #roundtrip = csv se export kia fir import karne pe data remains the same

@given(st.text(alphabet="ACGTNRUYKMSWBDHV", min_size= 1))
def test_csv_export_roundtrip(fasta):
    export = Exporter()
    seq = sequence("id", fasta)
    length = seq.sequence_length()
    counts = seq.base_count

    with tempfile.NamedTemporaryFile(mode = 'w', delete = False, suffix='.csv') as tmp_file:
        tmp_filepath = tmp_file.name
    
    try:
        export.to_csv(length, tmp_filepath)
    
        with open(tmp_filepath, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            loaded_data = list(reader)

        assert loaded_data == seq
    
    finally:
        os.unlink(tmp_filepath)

