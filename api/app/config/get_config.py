import os
import yaml
from easydict import EasyDict


def read_config(what="import_columns", file_name="config.yml"):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    config_path = os.path.join(script_dir, file_name)

    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    return EasyDict(config[what])
