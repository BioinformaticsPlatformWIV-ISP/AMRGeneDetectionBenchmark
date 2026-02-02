import binascii
from pathlib import Path


class FileSystemHelper(object):
    """
    This class contains utility function to work with the file system.
    """

    @staticmethod
    def is_gzipped(path: Path) -> bool:
        """
        Checks if the given file is compressed with gzip.
        :param path: Path
        :return: True if gzipped, False otherwise
        """
        with path.open('rb') as handle:
            magic_number = binascii.hexlify(handle.read(2))
        return magic_number == b'1f8b'

