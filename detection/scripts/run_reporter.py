#!/usr/bin/env python
import argparse
import glob
import json
import logging
from importlib.resources import files
from pathlib import Path
from typing import Any, Union

import pandas as pd
import yaml
from jinja2 import Environment, FileSystemLoader

from detection.app.reporter.genedetectionreporter import GeneDetectionReporter
from detection.app.utils.loggingutils import initialize_logging
from detection.config import load_tools_config
from detection.scripts import relative_to_absolute_path

REPORT_GR_TEMPLATE = Path(__file__).parent.parent / 'templates'


def _parse_arguments() -> argparse.Namespace:
    """
    Parses the command line arguments.
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir-in', type=relative_to_absolute_path, help='Input directory')
    parser.add_argument('--dir-working', type=relative_to_absolute_path, help='Working directory',
                        default=Path(Path.cwd(), 'reporter_working'))
    parser.add_argument('--dir-out', type=relative_to_absolute_path, help='Output directory',
                        default=Path(Path.cwd(), 'reporter_out'))
    parser.add_argument('--sample-name', type=str, help='Sample name', required=True)
    parser.add_argument('--ground-truth', type=relative_to_absolute_path, help='Ground truth')
    parser.add_argument('--configfile', type=relative_to_absolute_path, help='config file used by the main detection')
    parser.add_argument('--coverage-amr-genes', type=relative_to_absolute_path, help='Coverage amr genes file')
    parser.add_argument('--blast-ground-truth', action='store_true',
                        help='tick if you have blast-based ground truth [True]', default=True)
    parser.add_argument('--convert-gene-name', action='store_true', help='convert gene allele to gene name [True]',
                        default=True)
    return parser.parse_args()


def _check_arguments(arguments: argparse.Namespace) -> argparse.Namespace:
    """
    Checks if the provided command line arguments are valid, raises an exception otherwise.
    :param arguments: Parsed arguments
    :return: None
    """
    if not arguments.ground_truth:
        raise ValueError("Ground truth file not given.")
    if not arguments.coverage_amr_genes:
        raise ValueError("Coverage of the AMR genes file not given.")

    if not Path(arguments.ground_truth).exists():
        raise ValueError(f"Ground truth file {arguments.ground_truth} not found.")
    if not Path(arguments.coverage_amr_genes).exists():
        raise ValueError(f"Coverage AMR genes file {arguments.coverage_amr_genes} not found.")
    return arguments


def _format_gene_id(arguments: argparse.Namespace, gene: str, correspondence_dict: dict = None) -> str:
    """
    Formats the gene name nicely.
    :param arguments: Command line arguments
    :param gene: gene name
    :param correspondence_dict: allele to gene correspondence dictionary
    :return: formatted gene name
    """
    if arguments.convert_gene_name:
        try:
            gene = correspondence_dict[gene]['gene']
        except KeyError:
            pass
    else:
        try:
            gene = correspondence_dict[gene]['allele']
        except KeyError:
            pass
    return gene


def _retrieve_amr(arguments: argparse.Namespace, correspondence_dict: dict = None) -> dict[str, list[Any]]:
    """
    Retrieves the clean AMR output files.
    :param arguments: Command line arguments
    :param correspondence_dict: allele to gene correspondence dictionary
    :return: None
    """
    result = {}
    tools_config = load_tools_config()
    tool_correspondence = {tool: tools_config[tool]['clean_name'] for tool in tools_config.keys()}
    for json_output in glob.glob(f"{arguments.dir_in}/*/*_cleaned.json"):
        tool = tool_correspondence[Path(json_output).parent.name]
        result[tool] = []
        json_dictionary = json.load(open(json_output))
        for gene in json_dictionary:
            the_gene_id = _format_gene_id(arguments, gene['reference_accession'], correspondence_dict)
            result[tool].append(
                (the_gene_id, gene['input_sequence_id'], gene['coverage_percentage'],
                 gene['sequence_identity'], gene['coverage_depth']))
    return result


def _retrieve_blast_ground_truth(arguments: argparse.Namespace, genes_covered_df: pd.DataFrame, output_dir: Path,
                                 correspondence_dict: dict = None) -> pd.DataFrame:
    """
    Retrieves ground truth from blast-based approach.
    :param arguments: Command line arguments
    :param genes_covered_df: genes coverage in sequencing data
    :param output_dir: output directory
    :param correspondence_dict: NDARO allele to gene correspondence dictionary
    :return: Dataframe of ground truth AMR genes
    """
    genes_covered_dict = dict(zip(genes_covered_df['type_gene'], genes_covered_df['coverage']))
    with open(arguments.ground_truth) as handle:
        ground_truth_dataframe = pd.read_csv(handle, sep='\t', header=0, usecols=[1, 2, 3, 4],
                                             names=['gene', 'identity', 'coverage', 'genome'])
        ground_truth_dataframe_covered = ground_truth_dataframe[
            ground_truth_dataframe['gene'].isin(genes_covered_dict.keys())]
        ground_truth_dataframe_covered['seq_coverage'] = ground_truth_dataframe_covered['gene'].map(genes_covered_dict)
        ground_truth_dataframe_covered['coverage'] = list(
            map(lambda x: eval(x) * 100, ground_truth_dataframe_covered['coverage']))
        new_gene_id = []
        if arguments.convert_gene_name:
            for gene in ground_truth_dataframe_covered['gene']:
                new_gene_id.append([v['gene'] for k, v in correspondence_dict.items() if gene == v['allele']][0])
        ground_truth_dataframe_covered['gene'] = new_gene_id
    output_file = Path(output_dir, 'updated_ground_truth.tsv')
    ground_truth_dataframe_covered.to_csv(output_file, sep='\t', index=False)
    return ground_truth_dataframe_covered


def _retrieve_summary_for_precision_recall(arguments: argparse.Namespace,
                                           correspondence_dict: dict = None) -> dict[str, pd.DataFrame]:
    """
    Retrieves the summary file for precision and recall curves.
    :param arguments: Command line arguments
    :param correspondence_dict: AMRFinder correspondence dictionary
    :return: None
    """
    result = {}
    tools_config = load_tools_config()
    tool_correspondence = {tool: tools_config[tool]['clean_name'] for tool in tools_config.keys()}
    for tsv_output in glob.glob(f"{arguments.dir_in}/*/*_for_pr.tsv"):
        tool = tool_correspondence[Path(tsv_output).parent.name]
        pr_table = pd.read_table(tsv_output, sep='\t', header=0)
        pr_table['amr_gene'] = pr_table['amr_gene'].apply(
            lambda x: _format_gene_id(arguments, x, correspondence_dict))
        result[tool] = pr_table
    return result


def _generate_per_sample_report(graph_dict: dict[str, Union[Path, list[dict[str, Path]]]],
                                output_html: Path) -> None:
    """
    Generates the per sample report based on the jinja template.
    :param graph_dict: statistics dictionary
    :param output_html: Path to the HTML output report file
    :return: None
    """
    env = Environment(loader=FileSystemLoader(REPORT_GR_TEMPLATE))
    template = env.get_template('ground_truth_report.jinja')
    html_path = Path(output_html.parent.resolve())
    text_report = template.render(title="Per sample report",
                                  graph_PA=graph_dict['graph_PA'].relative_to(html_path),
                                  graph_stats=graph_dict['graph_stats'].relative_to(html_path),
                                  graph_dictionary=[{'src': p['src'].relative_to(html_path, walk_up=True)}
                                                    for p in graph_dict['graph_depth']],
                                  graph_precision_recall=[{'src': p['src'].relative_to(html_path, walk_up=True),
                                                           'src_reverse': p['src_reverse'].relative_to(html_path, walk_up=True)}
                                                          for p in graph_dict['precision_recall']])
    with open(output_html, 'w') as handle:
        handle.write(text_report)


def _generate_summary(stats_dictionary: dict, output_tsv: Path, config_dict: dict) -> None:
    """
    Generates a TSV summary file.
    :param stats_dictionary: statistics dictionary
    :param output_tsv: output tsv file
    :param config_dict: config dictionary
    :return: None
    """
    with open(output_tsv, 'w') as handle:
        for tool in stats_dictionary:
            to_write = (f'{tool}\tidentity_filter\t{config_dict["alignment_filtering"]["identity"]}\n'
                        f'{tool}\tlength_filter\t{config_dict["alignment_filtering"]["length"]}\n')
            handle.write(to_write)
            for stat in stats_dictionary[tool]:
                if isinstance(stats_dictionary[tool][stat], float):
                    to_write = f'{tool}\t{stat}\t{stats_dictionary[tool][stat]}\n'
                else:
                    to_write = f'{tool}\t{stat}\t{",".join(list(stats_dictionary[tool][stat]))}\n'
                handle.write(to_write)


def main() -> None:
    """
    Runs the single sample reporter.
    :return: None
    """
    initialize_logging()
    logging.info(f'Starting reporter script')

    # Parse the command line arguments
    args = _check_arguments(_parse_arguments())
    config_dictionary = yaml.safe_load(open(args.configfile))

    genes_covered_table = pd.read_table(args.coverage_amr_genes, usecols=['type_gene', 'coverage'])
    level = 'gene' if args.convert_gene_name else 'allele'
    output_directory = Path(f'{args.dir_out}/{args.sample_name}/reporter_{level}')
    report_directory = Path(f'{args.dir_out}/{args.sample_name}/reporter_{level}/report')
    if not report_directory.exists():
        report_directory.mkdir(parents=True)
    if args.blast_ground_truth:
        correspondence_file = Path(str(files('detection').joinpath('resources/ncbi_amr_correspondence.csv')))
        correspondence_dictionary = {line.strip().split(',')[0]: {'allele': line.strip().split(',')[1],
                                                                  'gene': line.strip().split(',')[2]}
                                     for line in correspondence_file.open().readlines()}
        ground_truth = _retrieve_blast_ground_truth(args, genes_covered_table, output_directory,
                                                    correspondence_dictionary)
    else:
        logging.warning('AMRFinder ground truth is not supported at the moment.')
        exit(1)

    amr_detected = _retrieve_amr(args, correspondence_dictionary)
    pr_values = _retrieve_summary_for_precision_recall(args, correspondence_dictionary)

    reporter = GeneDetectionReporter(output_directory, ground_truth, amr_detected, pr_values, args.dir_in, level)
    reporter.run()

    graph_dictionary = reporter.output_dictionary
    output_report = output_directory / 'report' / 'per_sample_report.html'
    _generate_per_sample_report(graph_dictionary, output_report)
    output_summary = output_directory / 'statistics_summary.tsv'
    _generate_summary(reporter.output_stats, output_summary, config_dictionary)


if __name__ == '__main__':
    main()
