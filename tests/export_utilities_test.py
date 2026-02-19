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


def test_csv_export():
    export = Exporter()
    sample_data = [
        {'name': 'Alice', 'age' : '30', 'city' : 'London'},
        {'name': 'Bob', 'age' : '25', 'city': 'NYC'}
    ]
    headers = ['name', 'age', 'city']
    with tempfile.NamedTemporaryFile(mode = 'w', delete = False, suffix='.csv') as tmp_file:
        tmp_filepath = tmp_file.name
    
    try:
        export.to_csv(sample_data, tmp_filepath)

        with open(tmp_filepath, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            rows = list(reader)
        
        assert len(rows) == 2
        assert rows[0]['name'] == 'Alice'
        assert rows[1]['city'] == 'NYC'
        assert reader.fieldnames == headers
    
    finally:
        os.unlink(tmp_filepath)
        