from pathlib import Path
import yaml


class Config:
    def __init__(self):
        self.config = self._get_defaults()
        self._config_source = "defaults"
        config_path = self._find_config_file()

        if config_path:
            print(
                f"Pre-existing configuration found, loading saved settings from {config_path}"
            )
            self.load_config(config_path)
            self._config_source = str(config_path)

        # else:
        #         self.create_config(config_path)
        # self._validate_config()

    def create_config(self, output_path):
        if output_path is None:
            output_path = Path(__file__).parent.parent / "config.yaml"
        if Path(output_path).exists():
            user_in = input(
                "A config file already exists. Would you like to overwrite it? (Y/N) "
            )
            if user_in != "Y":
                print("Keeping existing config.")
                return
        user_in = input(
            "Since no config files were found, the following options will allow you to set some general defaults according to your liking. \n"
            "Would you like to input these defaults? If not, the config file formed will contain defaults set beforehand. (Y/N) "
        )
        if user_in == "Y":
            defaults = self._get_defaults()

            strict_mode = input("Strict mode on (T/F): ")
            if strict_mode not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            strict_file = input("Strict file parsing (T/F): ")
            if strict_file not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            strict_sequence = input("Strict sequence parsing (T/F): ")
            if strict_sequence not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            case_sensitive = input("Case sensitive parsing (T/F): ")
            if case_sensitive not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            default_format = input("Default format for export (csv, tsv or json): ")
            if default_format not in ("csv", "tsv", "json"):
                raise ValueError(
                    "Value other than csv, tsv or json entered. You can export currently in csv, tsv or json only."
                )
            include_headers = input("Include headers in export? (T/F) ")
            if include_headers not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            gc_include_ambigcodes = input(
                "Include ambiguous codes during sequence statistics? (T/F) "
            )
            if gc_include_ambigcodes not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            count_unknown = input(
                "Count unknown characters during sequence statistics? (T/F) "
            )
            if count_unknown not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            start_codons_input = input(
                "Enter start codons separated by commas (e.g., ATG, GTG, TTG): "
            )
            start_codons = [
                codon.strip().upper()
                for codon in start_codons_input.split(",")
                if codon.strip()
            ]
            stop_codons_input = input(
                "Enter stop codons separated by commas (e.g., TAA, TAG, TGA): "
            )
            stop_codons = [
                codon.strip().upper()
                for codon in stop_codons_input.split(",")
                if codon.strip()
            ]
            unknown_codon_char = str(
                input(
                    "Please enter the character used to deal with unknown codons (e.g., X): "
                )
            )

            search_both_strands = input("Allow orf search on both strands? (T/F) ")
            if search_both_strands not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            allow_alternative_starts = input(
                "Allow orf search with alternate starts? (T/F) "
            )
            if allow_alternative_starts not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            nested_orfs = input("Allow nested orf search? (T/F) ")
            if nested_orfs not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            search_reverse_strand = input(
                "Allow search on reverse strand for motif? (T/F) "
            )
            if search_reverse_strand not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )
            iupac_support = input("Allow IUPAC support? (T/F) ")
            if iupac_support not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            enable_logging = input("Enable logging? (T/F) ")
            if enable_logging not in ("T", "F"):
                raise ValueError(
                    "Value other than T or F entered. Please enter only T or F."
                )

            user_inputs = {
                "parsing": {
                    "strict_mode": strict_mode,
                    "strict_file": strict_file,
                    "strict_sequence": strict_sequence,
                    "case_sensitive": case_sensitive,
                },
                "export": {
                    "default_format": default_format,
                    "include_headers": include_headers,
                },
                "statistics": {
                    "gc_include_ambigcodes": gc_include_ambigcodes,
                    "count_unknown": count_unknown,
                },
                "translation": {
                    "start_codons": start_codons,
                    "stop_codons": stop_codons,
                    "unknown_codon_char": unknown_codon_char,
                },
                "orf": {
                    "search_both_strands": search_both_strands,
                    "allow_alternative_starts": allow_alternative_starts,
                    "nested_orfs": nested_orfs,
                },
                "motif": {
                    "search_reverse_strand": search_reverse_strand,
                    "iupac_support": iupac_support,
                },
                "logging": {
                    "enabled": enable_logging,
                },
            }

            for inp in user_inputs.values():
                for key, val in inp.items():
                    if val == "T":
                        inp[key] = True
                    elif val == "F":
                        inp[key] = False

            self._deep_merge(defaults, user_inputs)
            self.config = defaults
            if output_path is None:
                output_path = Path(__file__).parent.parent / "config.yaml"
            output_path = Path(output_path)
            output_path.parent.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                yaml.dump(defaults, f, default_flow_style=False, sort_keys=False)
            print(f"Config file created at: {output_path}.")

        elif user_in == "N":
            if output_path is None:
                output_path = Path(__file__).parent.parent / "config.yaml"
            output_path = Path(output_path)
            output_path.parent.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                yaml.dump(
                    self._get_defaults(), f, default_flow_style=False, sort_keys=False
                )
            print(f"Config file created at: {output_path} using general defaults.")

        else:
            raise ValueError(
                "Value other than Y or N passed. Please enter either Y or N."
            )

    def _get_defaults(self):
        return {
            "parsing": {
                "strict_mode": False,
                "strict_file": False,
                "strict_sequence": False,
                "max_errors": 100,
                "case_sensitive": False,
                "allowed_characters": "ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.",
            },
            "output": {
                "decimal_places": 2,
                "cli_table_format": "grid",
                "export_formats": {"markdown": "pipe"},
                "use_rich": True,
                "color_scheme": "auto",
                "html_report": {
                    "template": "default",
                    "interactive": True,
                    "include_plots": True,
                },
                "verbose": False,
                "line_width": 80,
                "show_line_numbers": False,
                "uppercase_output": True,
            },
            "export": {
                "default_format": "csv",
                "default_delimiter": ",",
                "include_headers": True,
                "pretty_json": True,
                "json_indent": 2,
            },
            "statistics": {
                "gc_include_ambigcodes": False,
                "count_gaps": True,
                "count_unknown": True,
                "preserve_case": False,
            },
            "translation": {
                "genetic_code": "standard",
                "start_codons": ["ATG"],
                "stop_codons": ["TAA", "TAG", "TGA"],
                "unknown_codon_char": "X",
                "include_stop_in_protein": False,
            },
            "orf": {
                "min_length": 0,
                "min_protein_length": 0,
                "search_both_strands": True,
                "allow_alternative_starts": False,
                "nested_orfs": True,
            },
            "motif": {
                "default_k": 6,
                "max_mismatches": 0,
                "case_sensitive": False,
                "search_reverse_strand": True,
                "iupac_support": True,
            },
            "performance": {
                "max_file_size_mb": 1000,
                "chunk_size": 10000,
                "parallel_processing": False,
            },
            "cli": {
                "default_action": "summary",
                "auto_detect_format": True,
                "color_output": True,
                "progress_bar": False,
            },
            "files": {
                "auto_decompress": True,
                "encoding": "utf-8",
                "line_endings": "auto",
                "backup_on_export": True,
            },
            "logging": {
                "enabled": False,
                "log_file": "bio_seq.log",
                "log_level": "INFO",
                "log_format": "%(asctime)s - %(levelname)s - %(message)s",
            },
        }

    def _find_config_file(self):
        project_config = Path("./bio_seq_config.yaml")
        if project_config.exists():
            return project_config
        user_config = Path(__file__).parent.parent / "config.yaml"
        if user_config.exists():
            return user_config
        return None

    def _deep_merge(self, base, override):
        for key, value in override.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                self._deep_merge(base[key], value)
            else:
                base[key] = value

    def __repr__(self):
        return f"Config(loaded_from={self._config_source})"

    def load_config(self, config_path):
        try:
            with open(config_path, "r") as f:
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

    def get_config(self, key_path, default=None):
        keys = key_path.split(".")
        value = self.config
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        return value

    def set_config(self, key_path, value):
        keys = key_path.split(".")
        config = self.config
        for key in keys[:-1]:
            if key not in config:
                config[key] = {}
            config = config[key]
        config[keys[-1]] = value

    def _validate_config(self):
        """Validate config values are within acceptable ranges"""
        if self.config["parsing"]["max_errors"] < 0:
            print("Warning: max_errors cannot be negative, using 0")
            self.config["parsing"]["max_errors"] = 0
        if self.config["output"]["decimal_places"] < 0:
            print("Warning: decimal_places cannot be negative, using 2")
            self.config["output"]["decimal_places"] = 2
        if self.config["orf"]["min_length"] < 0:
            print("Warning: min_length cannot be negative, using 0")
            self.config["orf"]["min_length"] = 0
        if self.config["orf"]["min_protein_length"] < 0:
            print("Warning: min_protein_length cannot be negative, using 0")
            self.config["orf"]["min_protein_length"] = 0
        if self.config["motif"]["default_k"] < 0:
            print("Warning: default_k cannot be negative, using 0")
            self.config["motif"]["default_k"] = 0
