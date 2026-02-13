import json
import os

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
    
    return config
# Load the configuration
config = load_config()
