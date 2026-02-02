import abc
import json
import logging
from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.utils.command import Command


class BaseHandlerError(RuntimeError):
    """
    Error that can be raised when a base method fails.
    """
    pass


class BaseHandler(metaclass=abc.ABCMeta):
    """
    Base class for detection handlers.
    """

    def __init__(self, key: str, config: dict, tool_params: dict, version: str | None = None) -> None:
        """
        Initializes a tool handler.
        :param key: Tool key
        :param config: Tool configuration
        :param version: Tool version
        :return: None
        """
        self._key = key
        self._config = config
        self._tool_params = tool_params
        self._version: str = version if version is not None else self.get_version()
        logging.debug(f"{self.__class__.__name__} initialized for Tool '{key}'")

    @abc.abstractmethod
    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        return self._version

    @property
    def key(self) -> str:
        """
        Returns the tool key.
        :return: Tool key
        """
        return self._key

    def clean_read_header(self, input_fastq: Path, output_fastq_name: str, dir_running: Path) -> Path:
        """
        This method replaces the read headers from the input FASTQ file by numbered read names.
        Some tools have issues processing reads with special characters.
        :param input_fastq: input FASTQ file
        :param output_fastq_name: output FASTQ file name
        :param dir_running: Temporary directory
        :return: Path to the FASTQ file
        """
        output_fastq = Path(dir_running) / output_fastq_name
        command = Command(' '.join([
            f'zcat {input_fastq} | awk',
            ''''{print (NR%4 == 1) ? "@read_" ++i : $0}' | gzip -c > ''',
            f'{output_fastq}'
        ]))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise BaseHandlerError(f"Error running '{self.key}': {command.stderr}")
        return output_fastq

    def convert_fastq_to_fasta(self, input_fastq: Path, filename: str, dir_running: Path, **kwargs: Any) -> Path:
        """
        Converts the FASTQ input to FASTA.
        :param input_fastq: Input FASTQ file
        :param filename: Output FASTA file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: Path to the output FASTA file
        """
        output_fasta_file = Path(dir_running) / filename
        command = Command(' '.join([
            'seqkit fq2fa',
            f'--out-file {output_fasta_file}',
            f'{input_fastq}'
        ]), kwargs.get('dependencies', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise BaseHandlerError(f"Error running '{self.key}': {command.stderr}")
        return output_fasta_file

    def export_tsv_pr_file(self, input_tsv: pd.DataFrame, dir_out: Path, header: bool = True) -> None:
        """
        Creates a TSV output file with the format necessary for the Precision-Recall curves.
        :param header: output file contains a header
        :param input_tsv: input TSV file to write to dir_out
        :param dir_out: Output directory
        :return: None
        """
        logging.debug(f"Exporting TSV file to output directory {dir_out}")
        output_name = dir_out / f"{self._key}_for_pr.tsv"
        with open(output_name, "w") as handle:
            input_tsv.to_csv(handle, sep='\t', index=False, header=header)

    def export_json_file(self, input_json: list[dict], dir_out: Path) -> None:
        """
        Creates a JSON output file with the results.
        :param input_json: input JSON file to write to dir_out
        :param dir_out: Output directory
        :return: None
        """
        logging.debug(f"Exporting JSON file to output directory {dir_out}")
        output_name = dir_out / f"{self._key}_cleaned.json"
        with open(output_name, "w") as handle:
            json.dump(input_json, handle, indent=2)
