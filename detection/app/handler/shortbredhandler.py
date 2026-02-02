from importlib.resources import files
from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.handler.genedetectionhandler import GeneDetectionError, GeneDetectionHandler
from detection.app.utils.command import Command


class ShortbredHandler(GeneDetectionHandler):
    """
    Class to manage the AMR gene detection by ShortBRED.
    """

    key = 'shortbred'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command('shortbred_identify.py --version',
                          self._tool_params.get('dependencies', None),
                          self._tool_params.get('conda', None))
        command.execute(Path.cwd())
        return command.stderr.split('v')[1].strip()

    def run_gene_detection(self, input_fastq: Path, dir_running: Path, **kwargs: Any) -> Path:
        """
        Runs the ShortBRED gene detection.
        :param input_fastq: Input FASTQ file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: Path to the output file
        """
        nb_threads = self._config['resources']['threads_per_job']
        shortbred_table = Path(dir_running) / 'shortbred.tsv'
        command = Command(' '.join([
            'shortbred_quantify.py',
            f'--wgs {input_fastq}',
            f'--markers {self._tool_params["markers_path"]}',
            f'--results {shortbred_table}',
            f'--tmp {dir_running}',
            f'--threads {nb_threads}'
        ]), self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        return shortbred_table

    @staticmethod
    def clean_shortbred_like_hamronizer(input_table: Path, **kwargs: Any) -> list[dict]:
        """
        Manually cleans the ShortBRED output (as it is not available on hAMRonizer).
        :param input_table: shortbred output table
        :param kwargs: Additional arguments
        :return: List of records
        """
        output_table = pd.read_table(input_table, sep="\t", header=0,
                                     names=['reference_accession', 'coverage_depth', 'hits', 'length_amr_gene'])
        output_table['coverage_depth'] = output_table['coverage_depth'].astype(float)
        output_table['hits'] = output_table['hits'].astype(int)
        output_table['length_amr_gene'] = output_table['length_amr_gene'].astype(int)
        correspondence_file = Path(str(files('detection').joinpath('resources/shortbred_correspondence.csv')))
        correspondence_dict = {line.strip().split(',')[0]: line.strip().split(',')[1]
                               for line in correspondence_file.open().readlines()}
        output_table['reference_accession'] = output_table['reference_accession'].map(correspondence_dict)
        output_table = output_table[output_table['coverage_depth'] > 0]
        output_table['input_sequence_id'] = 'sequence_id'
        output_table['coverage_percentage'] = 95  # Placeholder
        output_table['sequence_identity'] = 70  # Placeholder
        output_pandas_table = output_table.to_dict(orient="records")
        return output_pandas_table[1:]

    def run(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any) -> None:
        """
        Runs ShortBRED and export the TSV and JSON output to the output directory.
        :param input_fastq: input fastq file
        :param dir_out: output directory
        :param dir_temp: running directory
        :param kwargs: additional arguments
        :return: None
        """
        cleaned_reads = self.clean_read_header(input_fastq, 'reads.fastq.gz', dir_temp)
        shortbred_output = self.run_gene_detection(cleaned_reads, dir_temp)
        cleaned_table = self.clean_shortbred_like_hamronizer(shortbred_output, **kwargs)
        self.export_json_file(cleaned_table, dir_out)
