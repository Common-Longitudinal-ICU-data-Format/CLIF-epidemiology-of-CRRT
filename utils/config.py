import json
import os
import sys

def load_config():
    # Get the directory where this file (config.py) is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up one level to project root, then into config folder
    project_root = os.path.dirname(current_dir)
    json_path = os.path.join(project_root, "config", "config.json")

    if os.path.exists(json_path):
        with open(json_path, 'r') as file:
            config = json.load(file)
        print("Loaded configuration from config.json")
    else:
        raise FileNotFoundError("Configuration file not found.",
                                "Please create config.json based on the config_template.")

    # Validate config fields
    code_dir = os.path.join(project_root, "code")
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)
    from pipeline_helpers import validate_config
    config = validate_config(config)

    return config
# Load the configuration
config = load_config()
