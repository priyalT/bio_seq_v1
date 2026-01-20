from pathlib import Path
from bio_seq_v1.fasta import FASTAParser
from bio_seq_v1.stats import sequence
from bio_seq_v1.translator import Translator
from bio_seq_v1.orf import ORFDetector, ORF
from bio_seq_v1.motif_search import MotifFinder, Match

class Exporter():
    def __init__(self, path):
        if path:
            self.path = Path(path) 
        else:
            self.path = Path(".")
    
    def to_csv(self):

    