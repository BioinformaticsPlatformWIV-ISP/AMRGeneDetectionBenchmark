import logging
from pathlib import Path
from typing import Optional

import yaml

from detection.app.utils.command import Command
from detection.config import config, load_tools_config


class ParallelToolHandler(object):
    """
    This class handles tools by running them in parallel using Snakemake.
    """

    def __init__(self, nb_jobs: int = 1, threads_per_job: int = 1,
                 min_identity: float = 0.9, min_coverage: float = 0.8) -> None:
        """
        Initializes the tool handler.
        :param nb_jobs: Number of jobs
        :param threads_per_job: Threads per job
        :param min_identity: Minimum identity
        :param min_coverage: Minimum coverage
        :return: None
        """
        self._nb_jobs = nb_jobs
        self._threads_per_job = threads_per_job
        self._min_identity = min_identity  # Placeholder
        self._min_coverage = min_coverage  # Placeholder
        self._path_config = None

    def run_tools(self, path_snakefile: Path, dir_out: Path, dir_logs: Path, dir_working: Path,
                  keys: Optional[list[str]] = None, tools_config: Optional[dict] = None) -> None:
        """
        Runs the tools in parallel.
        :param dir_working: working directory
        :param path_snakefile: Path to the main snakefile
        :param dir_out: Output directory
        :param dir_logs: Directory to store logs
        :param keys: (Optional) limit the updates to the tools with these keys
        :param tools_config: Tools configuration (optional), defaults to the content of the db.yml config file
        :return: None
        """
        # Update the keys
        tools_config = tools_config if tools_config is not None else load_tools_config()
        if keys is None:
            keys = tools_config.keys()
            logging.debug(f"'--keys' not specified, running all {len(keys)} tools")
        tools_config = {key: entry for key, entry in tools_config.items() if key in keys}
        ParallelToolHandler._validate_tool_configuration(tools_config)

        logging.debug(f'Running tools in working directory: {dir_working}')

        # Create working directory if it doesn't exist yet
        dir_working.mkdir(exist_ok=True, parents=True)

        # Create config file and snakefile
        self._path_config = self._create_config_file(dir_out, dir_working, dir_logs, tools_config, keys)

        # Run Snakefile
        try:
            command = self._run_snakefile(path_snakefile, self._path_config, dir_working)
            self._save_snakemake_logs(command, dir_logs)
        except RuntimeError as err:
            logging.error(f'Error occurred during Snakemake execution: {err}')
            raise err

    def _create_config_file(self, dir_out: Path, dir_working: Path, dir_logs: Path, tools_config: dict,
                            keys: Optional[list[str]] = None) -> Path:
        """
        Creates the temporary Snakemake config file to run the database updates.
        :param dir_out: Output directory
        :param dir_working: Working directory
        :param dir_logs: Directory to store logs
        :param tools_config: Tools list configuration
        :param keys: Limit the run to these tools
        :return: Path to config file
        """
        # Collect the config data
        config_data = {
            'tools': {key_tool: data for key_tool, data in tools_config.items() if key_tool in keys},
            'resources': config['resources'],
            'read_filtering': config['read_filtering'],
            'alignment_filtering': config['alignment_filtering'],
            'dir_out': str(dir_out),
            'dir_logs': str(dir_logs),
            'working_dir': str(dir_working)
        }
        path_config = dir_working / 'config.yml'
        config_data['resources']['threads_per_job'] = self._threads_per_job

        config_data['alignment_filtering']['coverage'] = self._min_coverage
        config_data['alignment_filtering']['identity'] = self._min_identity

        # Store in new file
        with path_config.open('w', encoding='utf-8') as handle:
            yaml.safe_dump(config_data, handle)
        logging.info(f'Snakemake config file created: {path_config}')
        return path_config

    def _run_snakefile(self, path_snakefile: Path, path_config: Path, dir_: Path) -> Command:
        """
        Runs the Snakefile with the given input
        :param path_snakefile: Path to the snakefile
        :param path_config: Path to the config file
        :param dir_: Directory to run the command
        :return: Executed Snakemake command
        """
        nb_cores = self._nb_jobs * self._threads_per_job
        logging.debug(f'Total nb. cores: {nb_cores} ({self._nb_jobs} jobs, {self._threads_per_job} cores per job)')
        command = Command("snakemake "
                          f'--snakefile {path_snakefile} '
                          f'--configfile {path_config} '
                          f'--cores {nb_cores}')
        command.execute(dir_, silent=True)
        logging.debug(f'Command working directory: {dir_}')
        if not command.returncode == 0:
            raise RuntimeError(f"Error executing Snakemake: {command.stderr}")
        return command

    @staticmethod
    def _save_snakemake_logs(command: Command, dir_logs: Path) -> None:
        """
        Saves the Snakemake logs in the output directory.
        :param command: Snakemake command
        :param dir_logs: Directory to store logs
        :return: None
        """
        logging.info(f'Storing Snakemake logs in: {dir_logs}')
        with (dir_logs / 'snakemake_stdout.txt').open('w') as handle:
            handle.write(command.stdout)
        with (dir_logs / 'snakemake_stderr.txt').open('w') as handle:
            handle.write(command.stderr)

    @staticmethod
    def _validate_tool_configuration(tools_config: dict) -> bool:
        """
        Checks if the given configuration is valid, raises an exception otherwise
        :param tools_config: Tool configuration to check
        :return: True if valid
        """
        pass

    def export_config_path(self) -> Path:
        """
        Export the path to the config file.
        :return: Path
        """
        return self._path_config
