import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotnine

from detection.app.reporter.basereporter import BaseReporter
from detection.app.utils.statutils import calculate_metrics


class GeneDetectionReporter(BaseReporter):
    """
    Base class for the gene detection reporter.
    """

    RANGE_THRESHOLDS = list(range(50, 105, 5))

    def detected_genes_to_readable_table(self, tool_id: str, ground_truth: pd.DataFrame,
                                         detected_amr: pd.DataFrame) -> None:
        """
        For the purpose of visual exploration, creates a readable csv.
        :param tool_id: name of the tool
        :param ground_truth: ground truth dataframe
        :param detected_amr: detected amr dataframe
        :return: None
        """
        output_directory = self.output_dir / f'{tool_id}'
        if not output_directory.exists():
            output_directory.mkdir(parents=True, exist_ok=True)
        merged_tables = pd.merge(ground_truth, detected_amr, how='outer', on='gene')
        merged_tables = merged_tables.sort_values(by=['depth'], ascending=False, na_position='last')
        merged_tables['depth'] = pd.to_numeric(merged_tables['depth'])
        values = []
        for index, row in merged_tables.iterrows():
            if (not isinstance(row['genome_x'], str)) and (row['depth'] > 0):
                values.append('FP')
            elif (isinstance(row['genome_x'], str)) and (row['depth'] > 0):
                values.append('TP')
            else:
                values.append('FN')
        merged_tables['classification'] = values
        output_name = self.output_dir / f'{tool_id}' / f'{tool_id}_for_comparison.csv'
        with open(output_name, "w") as handle:
            merged_tables.to_csv(handle, sep=';', index=False, header=list(merged_tables.columns))

    def retrieve_statistics(self) -> dict[str, dict]:
        """
        Generates the GT genes, TP genes, FP genes, FN genes; precision, recall, F1 score; genes
        associated with their depth
        :return: Dictionary with aforementioned entries
        """
        gt_genes = set(self.ground_truth['gene'])
        result = {}
        for tool in self.detected_amr:
            result_json = self.detected_amr[tool]
            pandas_detected_amr = pd.DataFrame.from_records(result_json,
                                                            columns=['gene', 'genome', 'coverage', 'identity', 'depth'])
            pandas_detected_amr[['coverage', 'identity', 'depth']] = (
                pandas_detected_amr[['coverage', 'identity', 'depth']].apply(pd.to_numeric))
            genes_found = set(pandas_detected_amr['gene'])
            genes_by_depth = dict(zip(pandas_detected_amr.gene, pandas_detected_amr.depth))
            precision, recall, f1_score = calculate_metrics(gt_genes, genes_found)
            true_positives = genes_found.intersection(gt_genes)
            false_positives = genes_found.difference(gt_genes)
            false_negatives = gt_genes.difference(genes_found)
            result[tool] = {'GT': gt_genes, 'TP': true_positives, 'FN': false_negatives, 'FP': false_positives,
                            'Precision': precision, 'Recall': recall, 'F1 score': f1_score,
                            'genes_by_depth': genes_by_depth}

            # Export results to nice CSV table
            self.detected_genes_to_readable_table(tool, self.ground_truth, pandas_detected_amr)
        self.output_stats = result
        return result

    def plot_presence_absence(self) -> None:
        """
        Plots the presence and absence as a white and black PAP matrix.
        :return: None
        """
        output_figure_path = self.output_dir / 'report' / 'presence_absence.png'
        p = plotnine.ggplot(self.ground_truth, plotnine.aes(x='genome', y='gene'))
        p += plotnine.geom_tile()
        p += plotnine.coord_equal(expand=True)
        p += plotnine.theme_bw()
        p += plotnine.theme(
            axis_text_x=plotnine.element_text(angle=45, hjust=1),
            axis_text_y=plotnine.element_text(style='italic'))
        p += plotnine.labs(title='Presence-Absence plot.')
        p.draw(show=True)
        p.save(output_figure_path, verbose=False, dpi=300, bbox_inches='tight', pad_inches=0)
        self.output_dictionary['graph_PA'] = output_figure_path

    def plot_statistics(self) -> None:
        """
        Plots the gene detection statistics (Precision, Recall, F1 score).
        :return: None
        """
        output_figure_path = self.output_dir / 'report' / 'statistics.png'
        stats = self.retrieve_statistics()
        stats_to_keep = ['Precision', 'Recall', 'F1 score']
        combined_stats = False
        for tool in stats:
            stats_filtered = {key: stats[tool][key] for key in stats_to_keep}
            pandas_stats = pd.DataFrame.from_dict(stats_filtered, orient='index', columns=["value"])
            pandas_stats.insert(1, 'statistic', pandas_stats.index.tolist())
            pandas_stats.insert(2, 'tool', str(tool))
            pandas_stats['value'] = pd.to_numeric(pandas_stats['value'])
            pandas_stats['statistic'] = pd.Categorical(pandas_stats['statistic'], stats_to_keep)

            if combined_stats is False:
                combined_stats = pandas_stats.copy()
            else:
                combined_stats = pd.concat([combined_stats, pandas_stats], ignore_index=True)
        combined_stats["statistic"] = combined_stats["statistic"].cat.reorder_categories(
            ["Precision", "Recall", "F1 score"])
        p = plotnine.ggplot(combined_stats, plotnine.aes(x='statistic', y='value', fill='tool'))
        p += plotnine.geom_col(position=plotnine.position_dodge(preserve='single'), color='black')
        p += plotnine.ylim(0, 1)
        p += plotnine.labs(x='Statistic', y='Value', fill='Tool')
        p += plotnine.scale_fill_brewer(type="qual", palette="Set2")
        p.draw(show=True)
        p.save(output_figure_path, verbose=False, dpi=300)
        self.output_dictionary['graph_stats'] = output_figure_path

    @staticmethod
    def text_in_tile(data_table: pd.DataFrame) -> list:
        """
        Automated function to generate text inside a tile from a data table.
        Adapted from the tutorial webpage of plotnine.
        :param data_table: pandas table with genes and depth to plot
        :return: list of layers of geom_text
        """
        layers = [
            plotnine.geom_text(data_table, plotnine.aes(label="gene"), nudge_y=0.2, size=6.5),
            plotnine.geom_text(data_table, plotnine.aes(label="depth"), nudge_y=-0.225, fontweight="normal", size=5.5),
        ]
        return layers

    @staticmethod
    def format_depth_table(data_table: pd.DataFrame, table_type: str = 'FP', tool_used: str = None) -> pd.DataFrame:
        """
        Create the depth table from the merged data frame.
        :param data_table: merged df
        :param table_type: FP or TP table
        :param tool_used: tool used to generate data
        :return: Depth table
        """
        if table_type == 'FP':
            depth_table = data_table[data_table['genome'].isnull()]
            depth_table['type'] = table_type
        else:
            depth_table = data_table
        depth_table['y_position'] = 1
        depth_table['x_position'] = 1
        count_position = 1
        number_of_columns = 20
        for i in range(0, len(depth_table.index), number_of_columns):
            depth_table['y_position'][i:i + number_of_columns] = count_position
            depth_table['x_position'][i:i + number_of_columns] = (
                list(range(1, number_of_columns + 1)))[:len(depth_table[i:i + number_of_columns].index)]
            count_position += 1
        if tool_used == 'ShortBRED':
            try:
                depth_table['depth'] = depth_table['depth'].apply(lambda x: f'{x:.2f}')
            except ValueError:
                pass
        return depth_table

    def create_dataframe_for_plotting(self, df: dict[str, Any], tool) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
        """
        Creates a pandas dataframe required for plotting TP/FN/FP square plots.
        :param df: input pandas dataframe
        :param tool: tool used to generate data
        :return: dataframes for the TP and FP + the counts in dictionary format
        """
        pandas_stats = pd.DataFrame(df['genes_by_depth'], index=['depth', ])
        pandas_stats_transposed = pandas_stats.transpose()
        pandas_stats_transposed.insert(1, 'gene', pandas_stats_transposed.index.tolist())
        merged_df = pd.merge(self.ground_truth, pandas_stats_transposed, on='gene', how='outer')

        # Create a DF with depth values from the detection
        genes_from_ground_truth = merged_df.dropna(subset=['genome'])
        genes_from_ground_truth['depth'] = genes_from_ground_truth['depth'].replace(np.nan, 0)
        genes_from_ground_truth['presence'] = 'TP'
        genes_from_ground_truth.presence = np.where(genes_from_ground_truth.depth.eq(0), 'FN',
                                                    genes_from_ground_truth.presence)
        genes_from_ground_truth_subset = genes_from_ground_truth.drop_duplicates(subset=['gene'])
        counts_tp_fn = genes_from_ground_truth_subset['presence'].value_counts().to_dict()
        genes_only_fp = self.format_depth_table(merged_df, tool_used=tool)
        genes_from_ground_truth = self.format_depth_table(genes_from_ground_truth_subset, 'TP', tool_used=tool)

        return genes_from_ground_truth, genes_only_fp, counts_tp_fn

    def plot_statistics_per_gene(self, tile_width: float = 0.95, tile_height: float = 0.95) -> None:
        """
        Plots the TP/FN for retrieved genes based on depth.
        :param tile_width: tile width for geom_tile
        :param tile_height: tile height for geom_tile
        :return: None
        """
        stats = self.retrieve_statistics()
        stats_to_keep = ['TP', 'FN', 'FP', 'genes_by_depth']
        self.output_dictionary['graph_depth'] = []
        for tool in stats:
            output_directory = self.output_dir / tool
            output_plot_tp_fn = output_directory / 'depth_by_ground_truth.png'
            output_plot_fp = output_directory / 'depth_by_false_positive.png'

            stats_filtered = {key: stats[tool][key] for key in stats_to_keep}
            genes_tp, genes_fp, counts_tp_fn = self.create_dataframe_for_plotting(stats_filtered, tool)

            # Starting the plot - TP/FN
            p = plotnine.ggplot(genes_tp,
                                plotnine.aes(x='x_position', y='y_position', fill='factor(presence)'))
            p += plotnine.geom_tile(plotnine.aes(width=tile_width, height=tile_height), color='black')
            p += self.text_in_tile(genes_tp)
            p += plotnine.coord_equal(expand=False)
            p += plotnine.theme(axis_text_x=plotnine.element_blank(),
                                axis_title_x=plotnine.element_blank(),
                                axis_text_y=plotnine.element_blank(),
                                axis_title_y=plotnine.element_blank(),
                                axis_ticks_y=plotnine.element_blank(),
                                axis_ticks_x=plotnine.element_blank())
            p += plotnine.scale_fill_manual(values={'TP': '#228B22', 'FN': '#FF0000'})
            p += plotnine.labs(fill='presence in sample', title=f'{tool} - TP/FN AMR gene coverage',
                               caption=f'#TP = {counts_tp_fn.get("TP", 0)}, #FN = {counts_tp_fn.get("FN", 0)}')
            p += plotnine.theme(plot_caption=plotnine.element_text(margin={"t": 6, "units": "lines"}),
                                figure_size=(8, 2))
            p.draw(show=True)
            p.save(output_plot_tp_fn, verbose=False, dpi=300)

            # Starting the plot - FP
            number_of_fps = len(genes_fp.index)
            try:
                p = plotnine.ggplot(genes_fp, plotnine.aes(x='x_position', y='y_position', fill='type'))
                p += plotnine.geom_tile(plotnine.aes(width=tile_width, height=tile_height), color='black')
                p += self.text_in_tile(genes_fp)
                p += plotnine.coord_equal(expand=False)
                p += plotnine.theme(axis_text_x=plotnine.element_blank(),
                                    axis_title_x=plotnine.element_blank(),
                                    axis_text_y=plotnine.element_blank(),
                                    axis_title_y=plotnine.element_blank(),
                                    axis_ticks_y=plotnine.element_blank(),
                                    axis_ticks_x=plotnine.element_blank(),
                                    legend_position='none')
                p += plotnine.scale_fill_manual(values={'FP': '#ADD8E6'})
                p += plotnine.labs(fill='type', title=f'{tool} - FP AMR gene coverage',
                                   caption=f'#FP = {number_of_fps}')
                p += plotnine.theme(plot_caption=plotnine.element_text(margin={"t": 6, "units": "lines"}),
                                    figure_size=(8, 2 * max(genes_fp['y_position'])))
                p.draw(show=True)
                p.save(output_plot_fp, verbose=False, limitsize=False, dpi=300)

            except (ZeroDivisionError, ValueError):
                logging.warning(f'{tool} - no false positives detected.')
                p = plotnine.ggplot(genes_fp, plotnine.aes(x='gene', y='1', fill='type'))
                p += plotnine.annotate("text", x=1, y=1, label="no FP")
                p.draw(show=True)
                p.save(output_plot_fp, verbose=False, dpi=300)

            self.output_dictionary['graph_depth'].append({'src': output_plot_tp_fn})
            self.output_dictionary['graph_depth'].append({'src': output_plot_fp})

    def create_precision_recall_table(self) -> dict[str, pd.DataFrame]:
        """
        Calculate precision and recall for multiple runs.
        :return: dictionary of precision, recall dataframes
        """
        gt_genes = set(self.ground_truth['gene'])
        output_dict = {}
        for tool, tsv in self.pr_json.items():
            tsv[['length_cov', 'identity']] = (
                tsv[['length_cov', 'identity']].apply(pd.to_numeric))
            dict_to_convert = {'Precision': [], 'Recall': [], 'Coverage_threshold': [],
                               'Identity_threshold': [], 'level': []}
            for ident_thresh in GeneDetectionReporter.RANGE_THRESHOLDS:
                for cov_thresh in GeneDetectionReporter.RANGE_THRESHOLDS:
                    sub_amr_cov_identity_file = tsv[(tsv['identity'] >= ident_thresh) &
                                                    (tsv['length_cov'] >= cov_thresh)]
                    genes_found = set(sub_amr_cov_identity_file['amr_gene'].tolist())
                    precision, recall, f1_score = calculate_metrics(gt_genes, genes_found)
                    dict_to_convert['Precision'].append(precision)
                    dict_to_convert['Recall'].append(recall)
                    dict_to_convert['Identity_threshold'].append(f'{ident_thresh}')
                    dict_to_convert['Coverage_threshold'].append(f'{cov_thresh}')
                    dict_to_convert['level'].append(self.level)
            tool_result_dataframe = pd.DataFrame.from_dict(dict_to_convert)
            tool_result_dataframe.apply(pd.to_numeric, errors='coerce')
            output_dict[tool] = tool_result_dataframe
            output_name = self.output_dir / f'{tool}' / f'{tool}_precision_recall_scores.csv'
            with open(output_name, "w") as handle:
                output_dict[tool].to_csv(handle, sep=';', index=False, header=list(output_dict[tool].columns))
        return output_dict

    @staticmethod
    def plot_precision_recall_curve(input_dataframe: pd.DataFrame, filename: Path, tool_name: str) -> None:
        """
        Plots the coverage and identity threshold in function of precision and recall
        :param input_dataframe: dataframe containing precision and recall values
        :param filename: filename to store the generated plot
        :param tool_name: tool used to generate the values
        :return: None
        """
        p = plotnine.ggplot(input_dataframe, plotnine.aes(x='Recall', y='Precision',
                                                          color='factor(Coverage_threshold)',
                                                          group='Coverage_threshold'))
        p += plotnine.geom_line(position=plotnine.position_dodge(.03))
        p += plotnine.ylim(0, 1)
        p += plotnine.xlim(0, 1)
        p += plotnine.geom_text(plotnine.aes(label='Identity_threshold'))
        p += plotnine.theme(figure_size=(9, 5))
        p += plotnine.labs(title=f'Precision-Recall curves for {tool_name}', color='Coverage')
        p += plotnine.scale_color_manual(
            values=['#543005', '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#f5f5f5', '#c7eae5', '#80cdc1',
                    '#35978f', '#01665e', '#003c30'])
        p += plotnine.labs(x='Recall', y='Precision')
        p.draw(show=True)
        p.save(filename, verbose=False, dpi=300)

    @staticmethod
    def plot_precision_recall_circle(input_dataframe: pd.DataFrame, filename: Path, tool_name: str) -> None:
        """
        Plots the precision and recall values in function of metrics
        :param input_dataframe: dataframe containing precision and recall values
        :param filename: filename to store the generated plot
        :param tool_name: tool used to generate the values
        :return: None
        """
        p = plotnine.ggplot(input_dataframe,
                            plotnine.aes(x='Coverage_threshold', y='Identity_threshold', size='Precision'))
        p += plotnine.geom_point(plotnine.aes(x='Coverage_threshold', y='Identity_threshold'), size=10,
                                 color='black', fill="lightgrey", stroke=1.0)
        p += plotnine.geom_point(plotnine.aes(color='Recall'))
        p += plotnine.scale_color_gradient(low="#f7fbff", high="#08306b", limits=(0, 1))
        p += plotnine.scale_size(limits=(0, 1))
        p += plotnine.theme(figure_size=(9, 5))
        p += plotnine.labs(title=f'Precision-Recall round plots for {tool_name}', color='Recall')
        p += plotnine.labs(x='Coverage (%)', y='Identity (%)')
        p.draw(show=True)
        p.save(filename, verbose=False, dpi=300)

    def plot_precision_recall_values(self) -> None:
        """
        Plots the PR curves and circle plots.
        :return: None
        """
        pr_json = self.create_precision_recall_table()
        self.output_dictionary['precision_recall'] = []
        for tool, pr_dataframe in pr_json.items():
            output_directory = self.output_dir / tool
            if not output_directory.exists():
                output_directory.mkdir(parents=True, exist_ok=True)

            pr_filename = output_directory / 'precision_recall.png'
            reversed_pr_filename = output_directory / 'precision_recall_reversed.png'

            pr_dataframe.sort_values('Precision', ascending=False, inplace=True)
            pr_dataframe.sort_values('Coverage_threshold', ascending=True, inplace=True)

            # Plots the standard precision-recall curve
            self.plot_precision_recall_curve(pr_dataframe, pr_filename, tool)

            # Now the other way around as circle plot
            self.plot_precision_recall_circle(pr_dataframe, reversed_pr_filename, tool)

            # Save the plots for the HTML report
            self.output_dictionary['precision_recall'].append({'src_reverse': reversed_pr_filename,
                                                               'src': pr_filename})

    def run(self) -> None:
        """
        Runs the reporter.
        :return: None
        """
        self.plot_presence_absence()
        self.plot_statistics()
        self.plot_statistics_per_gene()
        self.plot_precision_recall_values()
