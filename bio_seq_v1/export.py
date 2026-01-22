from pathlib import Path
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector, ORF
from bio_seq_v1.motif_search import MotifFinder, Match
import csv
import io
import json

class Exporter:

    @staticmethod
    def sequences_to_csv(sequence) -> str:
        if not sequence:
            return ""
        fieldnames = ["id", "length", "gc_content", "A", "C", "G", "T", "reverse complement"]
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=fieldnames)
        writer.writeheader()

        for seq in sequence:
            base = seq.base_count()
            writer.writerow({
                "id": seq.id,
                "length": seq.sequence_length(),
                "gc_content": seq.gc_content(),
                "A": base.get("A", 0),
                "C": base.get("C", 0),
                "G": base.get("G", 0),
                "T": base.get("T", 0),
            })
        return buffer.getvalue()
    
    @staticmethod
    def 

        
    