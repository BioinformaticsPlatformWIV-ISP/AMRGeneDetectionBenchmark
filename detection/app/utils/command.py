import logging
import subprocess
from pathlib import Path

from detection.config import config


class Command(object):
    """
    This class can be used to execute command line calls.
    """

    def __init__(self, command: str = None, dependencies: list[str] = None, conda: str = None) -> None:
        """
        Initializes a command.
        :param command: Command
        :param dependencies: Command dependencies
        :return: None
        """
        self._command = command
        self._stdout = None
        self._stderr = None
        self._returncode = None
        self._dependencies = dependencies
        self._conda = conda

    @property
    def command(self) -> str:
        """
        Returns the command line call.
        :return: Command line call
        """
        return self._command

    @command.setter
    def command(self, command: str) -> None:
        """
        Sets the command.
        :param command: Command
        :return: None
        """
        self._command = command

    @property
    def stdout(self) -> str:
        """
        Returns the stdout.
        :return: Stdout
        """
        return self._stdout

    @property
    def stderr(self) -> str:
        """
        Returns the stderr.
        :return: Stderr
        """
        return self._stderr

    @property
    def returncode(self) -> int:
        """
        Returns the return code of the command line call.
        :return: Return code
        """
        return self._returncode

    def execute(self, dir_working: Path, silent: bool = False) -> None:
        """
        Executes the command
        :param dir_working: Working directory
        :param silent: If True, stderr and stdout are not logged
        :return: None
        """
        cmd = self.command
        if self._dependencies is not None:
            cmd = f"module load {' '.join(self._dependencies)}; {cmd}"
        if self._conda is not None:
            cmd = f"{self._build_dependencies()}; {cmd}"
        logging.debug(f"Executing command '{cmd}'")
        procedure = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            executable='/bin/bash',
            cwd=str(dir_working))
        self._stdout = procedure.stdout.decode('utf-8')
        self._stderr = procedure.stderr.decode('utf-8')
        self._returncode = procedure.returncode

        if (not silent) and (self.stdout != ''):
            logging.debug(f'stdout: {self.stdout}')
        if (not silent) and (self.stderr != ''):
            logging.debug(f'stderr: {self.stderr}')

    def _build_dependencies(self) -> str:
        """
        Builds the dependencies.
        :return: Command to load dependencies
        """
        conda_root = config.get('conda_root')
        return '; '.join([
            f". {Path(conda_root) / 'etc' / 'profile.d' / 'conda.sh'}",
            f'conda activate {self._conda[0]}'
        ])
