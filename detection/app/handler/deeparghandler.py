from pathlib import Path
from typing import Any

import pandas as pd

from detection.app.handler.genedetectionhandler import GeneDetectionError, GeneDetectionHandler
from detection.app.utils.command import Command


class DeepargHandler(GeneDetectionHandler):
    """
    Class to manage the AMR gene detection by DeepARG.
    """

    key = 'deeparg'

    def get_version(self) -> str:
        """
        Retrieves the tool version.
        :return: Tool version
        """
        # DeepARG does not have a --version command, and it does not look like it will be updated anytime soon
        command = Command('deeparg -h',
                          self._tool_params.get('dependencies', None),
                          self._tool_params.get('conda', None))
        command.execute(Path.cwd())
        if 'command not found' in command.stderr:
            raise IndexError('Tool version not found')
        return '1.0.2'

    def run_gene_detection(self, input_fasta: Path, dir_running: Path, **kwargs: Any) -> Path:
        """
        Runs the DeepARG gene detection.
        :param input_fasta: Input FASTA file
        :param dir_running: Temporary directory
        :param kwargs: Additional arguments
        :return: Path to the output file
        """
        output_prefix = Path(dir_running) / 'deeparg'
        command = Command(' '.join([
            'deeparg predict',
            '--model SS',
            '--type nucl',
            f'--input {input_fasta}',
            f'--data-path {self._tool_params["model_path"]}',
            f'--output-file {output_prefix}'
        ]), self._tool_params.get('dependencies', None), self._tool_params.get('conda', None))
        command.execute(dir_running)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running '{self.key}': {command.stderr}")
        deeparg_table_path = str(Path(output_prefix))+'.mapping.ARG'
        deeparg_table = pd.read_table(deeparg_table_path, sep='\t')
        deeparg_table['best-hit'] = [x[4] for x in deeparg_table['best-hit'].str.split('|', expand=False)]

        output_file = str(Path(deeparg_table_path))+'_gene_name_cleaned'
        deeparg_table.to_csv(output_file, sep='\t', index=False)
        return Path(output_file)

    @staticmethod
    def account_for_depth_in_hamronizer(hamronized_dict: list[dict]) -> list:
        """
        Add depth from the total counts of genes in hAMRonized table.
        :param hamronized_dict: hAMRonized deeparg dictionary
        :return: hAMRonized dictionary with depth
        """
        all_amr_detected = [k['gene_name'] for k in hamronized_dict]
        all_amr_counts = {k: all_amr_detected.count(k) for k in set(all_amr_detected)}
        output_hamronized = []
        for entry in hamronized_dict:
            current_gene = entry['gene_name']
            if current_gene in all_amr_counts:
                entry['coverage_depth'] = all_amr_counts[current_gene]
                output_hamronized.append(entry)
                del all_amr_counts[current_gene]
        return output_hamronized

    def make_pr_dataframe(self, deeparg_cleaned_file: Path, **kwargs: Any) -> pd.DataFrame:
        """
        Generates a table necessary for Precision-Recall curves.
        :param deeparg_cleaned_file: Path to the deeparg cleaned file
        :param kwargs: Additional arguments
        :return: Pandas dataframe
        """
        gene_to_length_file = f'{self._tool_params["model_path"]}/database/v2/features.gene.length'
        gene_to_length_dictionary = {line.strip().split()[0].split('|')[4]: float(line.strip().split()[1]) for line
                                     in open(gene_to_length_file).readlines()}
        alignment_df = pd.read_table(deeparg_cleaned_file, sep='\t', usecols=[5, 7, 8], skiprows=1,
                                     names=['amr_gene', 'identity', 'align_length'])
        alignment_df['length_cov'] = [float(row['align_length']) / float(gene_to_length_dictionary[row['amr_gene']])
                                      for index, row in alignment_df.iterrows()]
        df = alignment_df[(alignment_df['identity'] >= 50) & (alignment_df['length_cov'] >= 0.5)]
        df = df[['amr_gene', 'identity', 'length_cov']]
        df['length_cov'] = df['length_cov'].apply(lambda x: f'{x * 100:.2f}')
        return df

    def run(self, input_fastq: Path, dir_out: Path, dir_temp: Path, **kwargs: Any) -> None:
        """
        Runs DeepARG and export the TSV and JSON output to the output directory.
        :param input_fastq: input fastq file
        :param dir_out: output directory
        :param dir_temp: running directory
        :param kwargs: additional arguments
        :return: None
        """
        fasta_filename = 'reads.fasta'
        fasta_reads = self.convert_fastq_to_fasta(input_fastq, fasta_filename, dir_temp, **self._config['read_filtering'])
        deeparg_output = self.run_gene_detection(fasta_reads, dir_temp)
        cleaned_table = self.run_hamronizer(DeepargHandler.key, deeparg_output, dir_temp)
        final_table = self.account_for_depth_in_hamronizer(cleaned_table)
        self.export_json_file(final_table, dir_out)
        pr_dataframe = self.make_pr_dataframe(deeparg_output)
        self.export_tsv_pr_file(pr_dataframe, dir_out)
