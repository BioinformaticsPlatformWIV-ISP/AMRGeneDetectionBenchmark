from importlib.resources import files
from pathlib import Path

import yaml


def _load_config() -> dict:
    """
    Loads the main configuration.
    :return: Config dictionary
    """
    _path_config = Path(str(files('detection').joinpath('config/config.yml')))
    if not _path_config.exists():
        raise FileNotFoundError(
            f"No config file found, update the 'config.yml.sample' file and move it to: {_path_config}")
    with open(_path_config) as handle:
        return yaml.load(handle, Loader=yaml.SafeLoader)


def load_tools_config() -> dict:
    """
    Loads the tools configuration.
    :return: Tools config dictionary
    """
    config_dict = _load_config()
    _path_config = Path(str(files('detection').joinpath('config/tools.yml')))
    if not _path_config.exists():
        raise FileNotFoundError(f"No tools config file found: {_path_config}")
    with open(_path_config) as handle:
        data_db = yaml.load(handle.read().format(
            db_directory=config_dict['database_directory'],
            argo_conda_env=config_dict['conda_env_names']['argo'],
            argpore_conda_env=config_dict['conda_env_names']['argpore'],
            deeparg_conda_env=config_dict['conda_env_names']['deeparg'],
            shortbred_conda_env=config_dict['conda_env_names']['shortbred'],
        ), Loader=yaml.SafeLoader)
    return data_db


def load_dataset_config() -> dict:
    """
    Loads the dataset configuration.
    :return: Dataset config dictionary
    """
    _path_config = Path(str(files('detection').joinpath('config/datasets.yml')))
    if not _path_config.exists():
        raise FileNotFoundError(f"No dataset config file found: {_path_config}")
    with open(_path_config) as handle:
        data_db = yaml.load(handle, Loader=yaml.SafeLoader)
    return data_db

config = _load_config()
