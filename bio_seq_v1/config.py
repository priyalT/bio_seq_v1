from pathlib import Path
import io

class Config:
    def __init__(self):
        config_path = Path("/bio_seq_v1/config.yaml")
        if config_path.exists():
            with config_path.open("r", encoding="utf-8") as config:
                content = config.read()
                print("Pre-existing configuration found, loading saved settings.")
                load_config()
        else:
            pass

    def load_config():
        config = {
            'parsing': {
                'strict_mode' : False,
                'strict_file' : False,
                'strict_sequence' : False,
                'max_errors' : 100,
                'case_sensitive' : False,
                'allowed_characters' : 'ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.'
            },
            'output' : {
                    'decimal_places': 2,
                    'cli_table_format' : 'grid',
                    'export_formats' : {
                        'markdown' : 'pipe'
                    },
                    'use_rich' : True,
                    'color_scheme' : 'auto',
                    'html_report' : {
                        'template' : 'default',
                        'interactive' : True,
                        'include_plots' : True
                    },
                    'verbose' : False,
                    'line_width' : 80,
                    'show_line_numbers' : False,
                    'uppercase_output' : True
            },
            'export' : {
                'default_format' : 'csv',
                'default_delimiter' : ',',
                'include_headers' : True,           
                'pretty_json' : True,               
                'json_indent': 2                  
            },
            'statistics': {
                'gc_include_ambigcodes' : False,
                'count_gaps' : True,              
                'count_unknown': True,
                'preserve_case': False            
            },
            'translation' : {
                'genetic_code' : 'standard',
                'start_codons': ['ATG'],
                'stop_codons': ['TAA', 'TAG', 'TGA'],
                'unknown_codon_char': 'X',
                'include_stop_in_protein': False
            },
            'orf' : {
                'min_length' : 0,                   
                'min_protein_length' : 0,           
                'search_both_strands': True,       
                'allow_alternative_starts': False, 
                'nested_orfs': True               
            }
            }
                  

        }
    
    def create_config():
        pass

    def get_config(name):
        pass

 