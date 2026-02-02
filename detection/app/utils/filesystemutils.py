import logging
import shutil
import tempfile
from pathlib import Path

import humanize
import requests

from detection.config import config


def create_temp_dir(prefix: str) -> tempfile.TemporaryDirectory:
    """
    Creates a temporary
    :param prefix: Directory prefix
    :return: Path to temporary directory
    """
    return tempfile.TemporaryDirectory(prefix=prefix, dir=config['dir_temp'])


def move_directory_contents(source: Path, target: Path) -> None:
    """
    Moves the content of the source directory to the target directory.
    :param source: Source directory
    :param target: Target directory
    :return: None
    """
    if not target.exists():
        target.mkdir(parents=True)
    for path in source.iterdir():
        if path.name.startswith('.'):
            continue
        if path.name.endswith('.fasta'):
            logging.debug(f"Moving: {path.name}")
        if path.is_file():
            shutil.move(str(path), str(target / path.name))
        elif path.is_dir():
            move_directory_contents(path, target / path.name)
