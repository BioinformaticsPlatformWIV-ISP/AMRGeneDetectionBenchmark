#!/usr/bin/env python
import argparse
import logging
import shutil
from datetime import datetime
from importlib.resources import files
from multiprocessing import cpu_count
from pathlib import Path
from typing import Any, Union, Callable

from detection.app.paralleltoolhandler import ParallelToolHandler
from detection.app.utils.command import Command
from detection.app.utils.filesystemhelper import FileSystemHelper
from detection.app.utils.loggingutils import initialize_logging
from detection.config import config, load_tools_config
from detection.scripts import relative_to_absolute_path


def _parse_arguments() -> argparse.Namespace:
    """
    Parses the command line arguments.
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq-se', type=relative_to_absolute_path, help='Input ONT FASTQ file')
    parser.add_argument('--sample-name', type=str, help='Sample name')
    parser.add_argument('--dir-out', type=relative_to_absolute_path, help='Directory to save results',
                        default=Path(Path.cwd(), 'main_out'))
    parser.add_argument('--dir-working', type=relative_to_absolute_path, help='Working directory',
                        default=Path(Path.cwd(), 'main_working'))
    parser.add_argument('--min-identity', type=float, help='Minimum identity percentage [0-100]',
                        choices=range(0, 101), default=90, metavar='')
    parser.add_argument('--min-coverage', type=float, help='Minimum length coverage [0-100]',
                        choices=range(0, 101), default=80, metavar='')
    parser.add_argument('--keys-list', type=relative_to_absolute_path,
                        help='TXT file containing a list of tools to run')
    parser.add_argument('--keys', help='Limit to the tools with these keys', nargs='+')
    parser.add_argument('--list', action='store_true', help='Lists the available tools, does NOT run')
    parser.add_argument('--nb-jobs', type=int, default=config['resources']['nb_jobs'],
                        help='Number of jobs to start in parallel')
    parser.add_argument('--threads-per-job', type=int, default=config['resources']['threads_per_job'],
                        help='Number of threads per job')
    parser.add_argument('--debug', action='store_true', help='Run in debug mode, does not delete intermediary files')
    parser.add_argument('--test', action='store_true', help='Run a test to check the current installation')
    return parser.parse_args()


def _check_arguments(arguments: argparse.Namespace) -> None:
    """
    Checks if the provided command line arguments are valid, raises an exception otherwise.
    :param arguments: Command line arguments
    :return: None
    """
    # Check if output directory is specified
    if (arguments.fastq_se is None) and (arguments.list is False):
        raise ValueError("FASTQ argument '--fastq-se' is required")

    # Check number of jobs
    if arguments.nb_jobs > cpu_count():
        logging.warning(f'Number of jobs ({arguments.nb_jobs}) exceeds the available CPUs ({cpu_count()})')
        logging.warning(f'Lowering number of parallel jobs to {cpu_count()}')
        arguments.nb_jobs = cpu_count()

    # Check the keys
    if arguments.keys is not None:
        if arguments.keys_list is not None:
            raise ValueError("arguments keys and keys-list are mutually exclusive")
        tools_config = load_tools_config()
        all_keys = arguments.keys.split(',')
        for key in all_keys:
            if key not in tools_config:
                raise ValueError(f"Key '{key}' not found in tools config")
        logging.info(f'Running {len(all_keys)} tools')


def _get_keys(arguments: argparse.Namespace) -> Union[list[str], None]:
    """
    Returns the keys of the tools that need to be run.
    :param arguments: Command line arguments
    :return: List of tool keys (or None when all tools should be run)
    """
    if arguments.keys is None and arguments.keys_list is None:
        # No keys specified -> return None
        return None
    elif arguments.keys is not None:
        keys_ = arguments.keys
    else:
        with arguments.keys_list.open() as handle:
            keys_ = [line.strip() for line in handle.readlines() if len(line.strip()) > 0]

    # Check if are present in the tools configuration
    tools_config = load_tools_config()
    for key in keys_:
        if key not in tools_config:
            raise ValueError(f"Key '{key}' not found in tools config")
    return keys_


def _list_keys_from_yml(load_function: Callable, info: str) -> None:
    """
    List the keys from a yml file.
    :param load_function: which load function to use
    :param info: type of entries to display
    :return: None
    """
    config_content = load_function()
    logging.info(f'Available {info}:')
    for key, db_data in sorted(config_content.items()):
        print(f'{key}')


def _symlink_input(arguments: argparse.Namespace, dir_out: Path) -> dict[str, Any]:
    """
    Symlinks the input files.
    :param arguments: Command line arguments
    :param dir_out: Directory to symlink to
    :return: Dictionary of FASTQ input
    """
    # Determine link names
    gzipped = FileSystemHelper.is_gzipped(arguments.fastq_se)
    links = {'input': arguments.fastq_se, 'symlink': f"input.fastq{'.gz' if gzipped else ''}"}

    # Create directory
    dir_links = dir_out / 'input'
    if not dir_links.exists():
        dir_links.mkdir(parents=True)

    # Link files
    paths_new = []
    path_orig, link_name = links['input'], links['symlink']
    path_new = dir_links / link_name
    logging.debug(f"Symlinking input file: {path_orig} -> {link_name}")
    if not path_new.is_symlink():
        path_new.symlink_to(path_orig)
    paths_new.append(path_new)

    # Return output dictionary
    return {'name': path_new.name, 'path': path_new}


def _test_installation(arguments: argparse.Namespace) -> None:
    """
    Test whether the required tools are installed.
    :param arguments: Command line arguments
    :return: None
    """
    from detection.app.handler import handler_by_key
    tool_config = load_tools_config()
    # First testing seqkit
    command_seqkit = Command('seqkit version')
    command_seqkit.execute(Path.cwd())
    if not command_seqkit.stdout.strip():
        raise RuntimeError(f"SeqKit is not available. Please verify your dependencies.")
    for tool in arguments.keys:
        try:
            handler_by_key[tool](key=tool, config=config, tool_params=tool_config[tool].get('params', {})).get_version()
        except IndexError:
            raise RuntimeError(f"Tool '{tool}' is not available. Please verify your dependencies.")
        if tool == 'shortbred':
            command_usearch = Command('usearch --version')
            command_usearch.execute(Path.cwd())
            if not command_usearch.stdout.strip():
                raise RuntimeError(f"Dependency usearch for ShortBRED is not available. Please verify your dependencies.")

    logging.info('Selected tools are seemingly installed. You can proceed forward.')


def main() -> None:
    """
    Runs the main detection.
    :return: None
    """
    initialize_logging()
    logging.info(f'Starting main script')

    # Parse the command line arguments
    args = _parse_arguments()

    # List the available tools if specified
    if args.list:
        _list_keys_from_yml(load_tools_config, 'tools')
        exit(0)

    keys = _get_keys(args)

    if args.test:
        _test_installation(args)
        exit(0)

    snakefile = Path(str(files('detection').joinpath('resources/main_gene_detection.smk')))

    # Create the logging directory
    dir_logs = Path(config['logging']['dir_logs'], datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    dir_logs.mkdir(exist_ok=True)
    logging.info(f'Storing logs in: {dir_logs}')

    output_directory = Path(args.dir_out) / args.sample_name

    fastq_symlinked = _symlink_input(args, output_directory)['path']

    # Run the database updates in parallel
    t_start = datetime.now()
    handler = ParallelToolHandler(args.nb_jobs, args.threads_per_job, args.min_identity, args.min_coverage)
    try:
        handler.run_tools(snakefile, output_directory, dir_logs, args.dir_working, keys)
    except RuntimeError as err:
        logging.error(f'Snakemake execution failed')
        raise err

    # Format the total duration as a string
    total_duration = datetime.now() - t_start
    total_duration_str = f'{int(total_duration.total_seconds()) // 60}m {int(total_duration.total_seconds()) % 60}s'
    logging.info(f'Total duration: {total_duration_str}')
    logging.info(f'Exporting the config file to the output directory')
    shutil.copy(handler.export_config_path(), output_directory)

    if args.debug:
        output_dir_logs = Path(output_directory) / 'logs'
        if not output_dir_logs.exists():
            output_dir_logs.mkdir(exist_ok=True)
        shutil.copytree(dir_logs, output_dir_logs, dirs_exist_ok=True)


if __name__ == '__main__':
    main()
