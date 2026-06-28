import os
import sys

def load_config():
    # Delegate to the single shared loader in pipeline_helpers, which honors the
    # CLIF_CONFIG env var (defaulting to config/config.json). Kept as a thin
    # wrapper so existing `from utils.config import config` imports still work.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    code_dir = os.path.join(project_root, "code")
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)
    from pipeline_helpers import load_config as _load_config
    config = _load_config()
    print(f"Loaded configuration from {config.get('_config_path', 'config.json')}")
    return config
# Load the configuration
config = load_config()
