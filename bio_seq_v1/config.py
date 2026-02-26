from pathlib import Path
import yaml


class Config:
    def __init__(self):
        self.config = self._get_defaults()
        config_path = self._find_config_file()
        if config_path:
                print(f"Pre-existing configuration found, loading saved settings from {config_path}")
                self.load_config(config_path)
        else:
                print("No configuration found, using defaults.")

    def _get_defaults(self):
        return {
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
            },
            'motif' : {
                'default_k' : 6,
                'max_mismatches': 0,
                'case_sensitive': False,
                'search_reverse_strand': True,
                'iupac_support': True
            },
            'performance' : {
                'max_file_size_mb' : 1000,
                'chunk_size' : 10000,
                'parallel_processing' : False,      
            },
            'cli' : {
                'default_action': 'summary',
                'auto_detect_format': True,
                'color_output': True,
                'progress_bar': False          
            },
            'files': {
                'auto_decompress': True,
                'encoding': "utf-8",
                'line_endings': "auto",
                'backup_on_export': True
            },
            'logging': {
                'enabled': False,
                'log_file': "bio_seq.log",
                'log_level': "INFO",
                'log_format': "%(asctime)s - %(levelname)s - %(message)s"
            }}
    
    def _find_config_file(self):
        project_config = Path("./bio_seq_config.yaml")
        if project_config.exists():
            return project_config
        user_config = Path.home() / '.bio_seq' / 'config.yaml'
        if user_config.exists():
            return user_config
        return None
    
    def _deep_merge(self, base, override):
        for key, value in override.items():
            if key in base and isinstance(base[key], dict) and isinstance(value,dict):
                self._merge_(base[key], value)
            else:
                base[key] = value
    
    def load_config(self, config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f)
            if user_config:
                self._deep_merge(self.config, user_config)
        except FileNotFoundError:
             print(f"Config file not found: {config_path}")
        except yaml.YAMLError as e:
             print(f"Error parsing YAML: {e}")
             print("Using default configuration.")
        except Exception as e:
             print(f"Unexpected error loading config: {e}")
             print("Using default configuration.")

    
    def create_config(self):
         

    def get_config(self, key_path, default=None):
         keys = key_path.split('.')
         value = self.config
         for key in keys:
              if isinstance(value, dict) and key in value:
                   value = value[key]
              else:
                   return default
         return value
    
    def set_config(self, key_path, value, default=None):
         pass

 