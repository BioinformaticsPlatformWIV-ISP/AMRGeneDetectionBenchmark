import abc
import logging
from pathlib import Path

import pandas as pd


class BaseReporter(object, metaclass=abc.ABCMeta):
    """
    Base class for report generation.
    """

    def __init__(self, output_dir: Path, ground_truth: pd.DataFrame, detected_amr: dict, pr_json: dict,
                 tool_dir: Path, level: str) -> None:
        """
        Initializes the base reporter.
        :param output_dir: output directory
        :param ground_truth: ground truth amr genes
        :param detected_amr: detected amr dictionary
        :param pr_json: precision-recall json
        :param tool_dir: tool directory
        :param level: report level (allele or gene)
        :return: None
        """
        self.tool_dir = tool_dir
        self.output_dir = output_dir
        self.ground_truth = ground_truth
        self.detected_amr = detected_amr
        self.pr_json = pr_json
        self.output_dictionary = dict()
        self.output_stats = dict()
        self.level = level
        logging.debug(f"{self.__class__.__name__} initialized.")
