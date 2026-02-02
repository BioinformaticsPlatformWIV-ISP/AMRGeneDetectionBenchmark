import glob
from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.handler.genedetectionhandler import GeneDetectionError, GeneDetectionHandler
from detection.app.utils.command import Command


class Argpore2Handler(GeneDetectionHandler):
    """
    Class to manage the AMR gene detection by ARGpore2
    """

    key = 'argpore2'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        command = Command('argpore.sh -h',
                          self._tool_params.get('dependencies', None),
                          self._tool_params.get('conda', None))
        command.execute(Path.cwd())
        return command.stdout.splitlines()[0].split(' ')[1].strip()

    def run_gene_detection(self, input_fasta: Path, dir_running: Path, **kwargs: Any) -> Path:
        """
        Runs the ARGpore2 gene detection.
        :param input_fasta: Input FASTA file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: output file path
        """
        nb_threads = self._config['resources']['threads_per_job']
        command = Command(' '.join([
            'argpore.sh',
            f'-f {input_fasta}',
            f'-t {nb_threads}',
        ]), self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        argpore2_output_file = [Path(f) for f in glob.glob(f"{dir_running}/*summary*")][-1]
        return argpore2_output_file

    def clean_argpore2_like_hamronizer(self, input_file: Path, **kwargs: Any) -> list[dict]:
        """
        Manually cleans ARGpore2 output (as it is not available on hAMRonizer).
        :param input_file: output file from ARGpore2
        :param kwargs: Additional arguments
        :return: List of records
        """
        output_table = pd.read_table(input_file, sep="\t", header=None,
                                     names=['reference_accession', 'type', 'coverage_depth'])

        output_table['sequence_identity'] = self._config['alignment_filtering']['identity'] / 100
        output_table['sequence_identity'] = 0.9  # Placeholder
        output_table['coverage_percentage'] = 0.75  # Placeholder
        output_table['input_sequence_id'] = 'sequence_id'

        output_pandas_table = output_table.to_dict(orient="records")
        return output_pandas_table[1:]

    @staticmethod
    def make_pr_dataframe(dir_running: Path) -> pd.DataFrame:
        """
        Generates a table necessary for Precision-Recall curves.
        :param dir_running: Running directory
        :return: Pandas dataframe
        """
        last_arg = glob.glob(f"{dir_running}/*arg.pr.tab")[-1]
        last_arg_df = pd.read_table(last_arg, sep="\t", usecols=[2, 3, 4, 13], header=0,
                                    names=['amr_gene', 'identity', 'length_align', 'length_amr_gene'])
        last_arg_df.apply(pd.to_numeric, errors='ignore')
        last_arg_df['length_cov'] = last_arg_df['length_align'] / last_arg_df['length_amr_gene']
        df = last_arg_df[(last_arg_df['identity'] >= 50) & (last_arg_df['length_cov'] >= 0.5)]
        df['length_cov'] = df['length_cov'].apply(lambda x: f'{x * 100:.2f}')
        df = df[['amr_gene', 'identity', 'length_cov']]
        df['amr_gene'] = [x[5] for x in df['amr_gene'].str.split('|', expand=False)]
        return df

    def run(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any) -> None:
        """
        Runs ARGpore2 and exports the clean table.
        :param input_fastq: input fastq file
        :param dir_out: output directory
        :param dir_temp: running directory
        :param kwargs: additional arguments
        :return: None
        """
        fasta_filename = 'reads.fasta'
        fasta_reads = self.convert_fastq_to_fasta(input_fastq, fasta_filename, dir_temp, **self._config['read_filtering'])
        argpore2_file = self.run_gene_detection(fasta_reads, dir_temp)
        cleaned_table = self.clean_argpore2_like_hamronizer(argpore2_file)
        self.export_json_file(cleaned_table, dir_out)
        pr_dataframe = self.make_pr_dataframe(dir_temp)
        self.export_tsv_pr_file(pr_dataframe, dir_out)
