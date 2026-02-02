import logging
from pathlib import Path

from detection.config import config


def initialize_logging() -> None:
    """
    Initializes the logging.
    :return: None
    """
    logging.basicConfig(
        level=logging.DEBUG,
        format=config['logging']['format']
    )


def add_filehandler(dir_logs: Path, tool_key: str) -> logging.Handler:
    """
    Adds a file handler to the default logger.
    :param dir_logs: Directory to store logs
    :param tool_key: Tool key
    :return: Handler
    """
    handler = logging.FileHandler(str(dir_logs / f'{tool_key}.log'))
    handler.setFormatter(logging.Formatter(config['logging']['format']))
    logging.getLogger().addHandler(handler)
    return handler
