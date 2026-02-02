import glob
from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.handler.genedetectionhandler import GeneDetectionError, GeneDetectionHandler
from detection.app.utils.command import Command


class ArgoHandler(GeneDetectionHandler):
    """
    Class to manage the AMR gene detection by Argo
    """

    key = 'argo'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command('argo --version',
                          self._tool_params.get('dependencies', None),
                          self._tool_params.get('conda', None))
        command.execute(Path.cwd())
        return command.stdout.strip()

    def run_gene_detection(self, input_fastq: Path, dir_running: Path, **kwargs: Any) -> Path:
        """
        Runs the Argo gene detection.
        :param input_fastq: Input FASTQ file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: Argo output file path
        """
        nb_threads = self._config['resources']['threads_per_job']
        output_prefix = Path(dir_running) / 'argo'
        command = Command(' '.join([
            'argo',
            f'--db {self._tool_params["db"]}',
            f'--output {output_prefix}',
            f'--threads {nb_threads}',
            f'{input_fastq}',
        ]), self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        argo_output_file = [Path(f) for f in glob.glob(f"{dir_running}/argo/*sarg.tsv")][0]
        return argo_output_file

    @staticmethod
    def clean_argo_like_hamronizer(input_file: Path, **kwargs: Any) -> list[dict]:
        """
        Manually cleans Argo output (as it is not available on hAMRonizer).
        :param input_file: Output file from Argo
        :param kwargs: Additional arguments
        :return: List of records
        """
        output_table = pd.read_table(input_file, sep="\t", header=0,
                                     names=['lineage', 'type', 'reference_accession',
                                            'coverage_depth', 'genome', 'abundance'])

        output_table['sequence_identity'] = 0.9  # Placeholder
        output_table['coverage_percentage'] = 0.9  # Placeholder
        output_table['input_sequence_id'] = 'sequence_id'

        output_pandas_table = output_table.to_dict(orient="records")
        return output_pandas_table[1:]

    def run(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any) -> None:
        """
        Runs Argo and exports the clean table.
        :param input_fastq: input fastq file
        :param dir_out: output directory
        :param dir_temp: running directory
        :param kwargs: additional arguments
        :return: None
        """
        argo_output_table = self.run_gene_detection(input_fastq, dir_temp)
        cleaned_table = self.clean_argo_like_hamronizer(argo_output_table, **kwargs)
        self.export_json_file(cleaned_table, dir_out)
