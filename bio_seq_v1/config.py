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
        
    
    def create_config():
        pass

    def get_config(name):
        pass

 