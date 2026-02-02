import logging
from datetime import datetime
from pathlib import Path

from detection.app.handler import handler_by_key
from detection.app.handler.genedetectionhandler import GeneDetectionError
from detection.app.utils import loggingutils
from detection.app.utils.command import Command
from detection.app.utils.loggingutils import initialize_logging

initialize_logging()

rule run_all:
    """
    Rule to ensure that all jobs are ran.
    """
    input:
        TSV = expand(Path(config['dir_out'], '{tool}', '{tool}_cleaned.json'), tool=config['tools'].keys())


rule filter_reads:
    """
    Filters reads based on quality and length.
    """
    input:
        FASTQ = Path(config['dir_out']) / 'input' / 'input.fastq.gz'
    output:
        FASTQ = Path(config['dir_out']) / 'input' / 'fastq_filtered.fastq.gz'
    params:
        length_filter = config['read_filtering']['length'],
        quality_filter = config['read_filtering']['quality'],
        output_directory = Path(config['dir_out']) / 'input'
    threads: 36
    run:
        command = Command(' '.join([
            'seqkit seq',
            f'-m {params.length_filter}',
            f'-Q {params.quality_filter}',
            f'-j {threads}',
            f'--out-file {output.FASTQ}',
            f"{input.FASTQ}"
        ]))
        command.execute(params.output_directory)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running seqkit seq: {command.stderr}")

rule reads_statistics:
    """
    Computes statistics on input reads.
    """
    input:
        FASTQ = rules.filter_reads.output.FASTQ
    output:
        TSV = Path(config['dir_out']) / 'input' / 'fastq_filtered_stats.tsv'
    params:
        output_directory = Path(config['dir_out']) / 'input'
    run:
        command = Command(' '.join([
            'seqkit stats -aT',
            f'{input.FASTQ}',
            f'> {output.TSV}'
        ]))
        command.execute(params.output_directory)
        if not command.returncode == 0:
            raise GeneDetectionError(f"Error running seqkit stats: {command.stderr}")


rule run_tool:
    """
    Runs a single tool and export a clean output file.
    """
    input:
        FASTQ = rules.filter_reads.output.FASTQ,
        FASTQ_STATS = rules.reads_statistics.output.TSV
    output:
        JSON = Path(config['dir_out']) / '{tool}' / '{tool}_cleaned.json'
    threads: 72
    params:
        tool = lambda wildcards: wildcards.tool,
        tool_params = lambda wildcards: config['tools'][wildcards.tool].get('params', {}),
        handler = lambda wildcards: wildcards.tool,
        dir_out = lambda wildcards: Path(config['dir_out'], wildcards.tool),
        dir_logs = config['dir_logs'],
        working_directory = config['working_dir']
    run:
        # Attach a logger
        log_handler = loggingutils.add_filehandler(Path(params.dir_logs), str(params.tool))

        # Save the start time
        t0 = datetime.now()

        # Initializes tool handler
        handler = handler_by_key[params.handler](str(params.tool), config, params.tool_params)

        try:
            # Attempt to run the gene detection tool handler
            handler.run(Path(input.FASTQ), Path(str(params.dir_out)), params.working_directory)
        except BaseException as err:
            message = str(err)
            logging.error(err, exc_info=True)
        finally:
            # Remove logger
            logging.getLogger().removeHandler(log_handler)
