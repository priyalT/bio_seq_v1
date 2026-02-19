from pathlib import Path
import csv
import io
import json
import sys

class Exporter:
    
    @staticmethod
    def _write_or_print(content: str, file_path = None):
        if file_path:
            Path(file_path).parent.mkdir(parents=True, exist_ok=True)
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(content)
        else:
            sys.stdout.write(content)
    
    @staticmethod
    def to_csv(data, file_path=None, delimiter=","):
        if not data:
            Exporter._write_or_print("", file_path)
            return

        if isinstance(data[0], dict) and 'fieldnames' in data[0]:
            fieldnames = data[0]['fieldnames']
            rows = data[1:]
        else:
            fieldnames = list(data[0].keys()) 
            rows = data

        buffer = io.StringIO()
        writer = csv.DictWriter(
            buffer,
            fieldnames=fieldnames,
            delimiter=delimiter
        )
        writer.writeheader()
        writer.writerows(rows)
        Exporter._write_or_print(buffer.getvalue(), file_path)

    @staticmethod
    def to_tsv(data, file_path=None):
        Exporter.to_csv(data, file_path, delimiter='\t')

    @staticmethod
    def to_json(data, file_path=None):
        content = json.dumps(data, indent = 2)
        Exporter._write_or_print(content, file_path)

    @staticmethod
    def sequences_to_csv(sequences, file_path=None):
        rows = [
            {
                "id": seq.id,
                "length": seq.sequence_length(),
                "gc content": seq.gc_content(),
                "sequence": seq.sequence,
                "reverse complement": seq.rev_complement()
            }
            for seq in sequences
        ]
        Exporter.to_csv(rows, file_path)

    @staticmethod
    def orfs_to_csv(orfs, file_path=None):
        rows = [
            {
                "sequence id": orf.seq_id,
                'start': orf.start,
                'end': orf.end,
                'frame': orf.frame,
                'strand': orf.strand,
                'length': orf.length,
                'protein': orf.protein
            }
            for orf in orfs
        ]
        Exporter.to_csv(rows, file_path)

    @staticmethod
    def motifs_to_csv(matches, file_path=None):
        rows=[
            {
                "id" : m.id,
                "position" : m.position,
                'matched_seq' : m.matched_seq,
                'strand_attributes' : m.strand_attributes
            }
            for m in matches
        ]
        Exporter.to_csv(rows, file_path)

    @staticmethod
    def to_fasta(sequences, file_path=None):
        lines = []
        for seq in sequences:
            lines.append(f">{seq.id}")
            lines.append(seq.sequence)
        content = "\n".join(lines) + "\n"
        Exporter._write_or_print(content,file_path)

    