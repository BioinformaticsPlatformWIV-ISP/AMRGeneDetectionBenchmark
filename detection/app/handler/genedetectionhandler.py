import abc
import json
from pathlib import Path
from typing import Any

from detection.app.handler.basehandler import BaseHandler
from detection.app.utils.command import Command


class GeneDetectionError(RuntimeError):
    """
    Error that can be raised when a gene detection fails.
    """
    pass


class GeneDetectionHandler(BaseHandler):

    @abc.abstractmethod
    def run_gene_detection(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any) -> Path:
        """
        Runs the gene detection module.
        :param input_fastq: input FASTQ file
        :param dir_out: Directory to store the gene detection results
        :param dir_temp: Directory to store files temporarily (is automatically deleted afterward)
        :param kwargs: Additional keyword arguments
        :return: Tool output file path
        """
        pass

    def run_hamronizer(self, tool_key: str, output_tool_file: Path, dir_running: Path, **kwargs: Any) -> list[dict]:
        """
        Runs hAMRonizer.
        :param output_tool_file: Output file to hAMRonize
        :param tool_key: tool name
        :param dir_running: Running directory
        :param kwargs: Additional arguments
        :return: None
        """
        db_version = self._tool_params['db_version']
        tool_version = self._version
        input_file_name_for_hamronizer = f'{tool_key}_hamronized'
        output_file_name = Path(dir_running) / 'output_hamronized.json'
        command = Command(' '.join([
            f'hamronize {tool_key}',
            '--format json',
            f'--analysis_software_version {tool_version}',
            f'--reference_database_version {db_version}',
            f'--input_file_name {input_file_name_for_hamronizer}',
            f'--output {output_file_name}',
            f'{output_tool_file}'
        ]), self._tool_params.get('dependencies', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        output_pandas_table = json.load(open(output_file_name))
        return output_pandas_table
