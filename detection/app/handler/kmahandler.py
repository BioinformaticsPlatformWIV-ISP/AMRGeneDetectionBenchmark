from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.handler.genedetectionhandler import GeneDetectionError, GeneDetectionHandler
from detection.app.utils.command import Command


class KMAHandler(GeneDetectionHandler):
    """
    Class to manage the AMR gene detection by KMA
    """

    key = 'kma'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command('kma -v', self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(Path.cwd())
        return command.stdout.split('-')[1].strip()

    def run_gene_detection(self, input_fastq: Path, dir_running: Path, **kwargs: Any) -> Path:
        """
        Runs the KMA gene detection.
        :param input_fastq: Input FASTQ file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: The KMA output file
        """
        nb_threads = self._config['resources']['threads_per_job']
        output_prefix = Path(dir_running) / 'output'
        command = Command(' '.join([
            f'kma -i {input_fastq}',
            f'-o {output_prefix}',
            f'-t_db {self._tool_params["db"]}',
            f'-t {nb_threads}',
            '-bcNano -ont -nc -nf -na'
        ]), self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        return Path(dir_running) / 'output.res'

    @staticmethod
    def filter_kma_output(kma_table: Path, dir_running: Path, min_identity: float = 0,
                          min_length_coverage: float = 0) -> Path:
        """
        Cleans the gene names and filters the KMA output based on identity and coverage (placeholder).
        :param kma_table: Original KMA output file
        :param dir_running: Temporary directory
        :param min_identity: Minimum identity to filter
        :param min_length_coverage: Minimum coverage to filter
        :return: Path to the output KMA table
        """
        pd_table = pd.read_table(kma_table)
        pd_table['#Template'] = [x[5] for x in pd_table['#Template'].str.split('|', expand=False)]
        filtered_pd_table = pd_table.loc[(pd_table['Query_Identity'] >= min_identity)
                                         & (pd_table['Template_Coverage'] >= min_length_coverage)]
        output_filtered_table = Path(dir_running) / 'output_filtered.res'
        filtered_pd_table.to_csv(output_filtered_table, sep='\t', index=False)
        return output_filtered_table

    @staticmethod
    def make_pr_dataframe(kma_table: Path) -> pd.DataFrame:
        """
        Generates a table necessary for Precision-Recall curves.
        :param kma_table: Original KMA output file
        :return: Pandas dataframe
        """
        pd_table = pd.read_table(kma_table, usecols=['#Template', 'Query_Identity', 'Template_Coverage', 'Depth'])
        pd_table.columns = ['amr_gene', 'identity', 'length_cov', 'depth']
        pd_table['amr_gene'] = [x[5] for x in pd_table['amr_gene'].str.split('|', expand=False)]
        return pd_table

    def run(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any):
        """
        Runs KMA and export the TSV and JSON output to the output directory.
        :param input_fastq: Input fastq file
        :param dir_out: Output directory
        :param dir_temp: Temporary directory
        :param kwargs: Additional arguments
        :return: None
        """
        kma_output_table = self.run_gene_detection(input_fastq, dir_temp)
        kma_filtered_table = self.filter_kma_output(kma_output_table, dir_temp)
        cleaned_table = self.run_hamronizer('kmerresistance', kma_filtered_table, dir_temp, **kwargs)
        self.export_json_file(cleaned_table, dir_out)
        pr_dataframe = self.make_pr_dataframe(kma_output_table)
        self.export_tsv_pr_file(pr_dataframe, dir_out)
