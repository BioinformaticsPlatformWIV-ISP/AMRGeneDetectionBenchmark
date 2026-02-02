from pathlib import Path


def relative_to_absolute_path(path: str) -> Path:
    """
    Takes a relative path and returns the absolute path.
    :param path: Relative or absolute path
    :return: The resolved absolute Path object
    """
    return Path(path).expanduser().resolve()
