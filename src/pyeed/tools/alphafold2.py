"""
AlphaFold2 runner module.

This module provides functionality to run AlphaFold2 predictions on protein sequences.
AlphaFold2 must be installed in a conda environment named 'alphafold_env' as described in:
https://github.com/google-deepmind/alphafold/tree/main
"""

import logging
import os
import subprocess
from pathlib import Path

import torch

logger = logging.getLogger(__name__)


class AlphaFoldRunner:
    """Class to manage and execute AlphaFold2 protein structure predictions."""

    def __init__(self, data_dir: str, output_dir: str) -> None:
        """
        Initialize the AlphaFold runner with required directories.

        Args:
            data_dir: Path to the directory containing AlphaFold model data
            output_dir: Path where prediction results will be stored

        Raises:
            FileNotFoundError: If required paths or files are not found
            EnvironmentError: If no GPU is detected
        """
        # Get the base directory of the pyeed project
        self.base_dir = Path(os.path.dirname(os.path.abspath(__file__)))
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)

        # Set path to the docker run script
        self.docker_script = self.base_dir / "resources/alphafold/docker_run.py"

        # Validate required paths
        if not self.docker_script.exists():
            raise FileNotFoundError(
                f"Docker run script not found: {self.docker_script}"
            )
        if not self.data_dir.exists():
            raise FileNotFoundError(
                f"AlphaFold data directory not found: {self.data_dir}"
            )
        if not self.output_dir.exists():
            self.output_dir.mkdir(parents=True)

        # Verify GPU availability
        if not torch.cuda.is_available():
            raise EnvironmentError("No GPU detected. AlphaFold requires a GPU to run.")

        logger.info("GPU detected. AlphaFold will run on GPU.")

    def run_alphafold(
        self, sequence: str, sequence_id: str, max_template_date: str = "2022-01-01"
    ) -> dict[str, str]:
        """
        Run AlphaFold prediction on a given protein sequence.

        Args:
            sequence: The protein sequence to predict structure for
            sequence_id: Unique identifier for the sequence
            max_template_date: Latest date allowed for template structures (YYYY-MM-DD format)

        Returns:
            dict[str, str]: Contains AlphaFold confidence score and structure path with
                keys 'confidence_score' and 'structure_path'

        Raises:
            RuntimeError: If AlphaFold execution fails
        """
        # Sanitize sequence ID by replacing dots with underscores
        sequence_id = sequence_id.replace(".", "_")

        # Create FASTA file for the sequence
        fasta_path = self.output_dir / f"{sequence_id}.fasta"
        with open(fasta_path, "w", encoding="utf-8") as fasta_file:
            fasta_file.write(f">{sequence_id}\n{sequence}")
        logger.info("Created FASTA file at: %s", fasta_path)

        # Construct AlphaFold command
        cmd = [
            "source ~/anaconda3/etc/profile.d/conda.sh && "
            "conda activate alphafold_env && "
            f"python {self.docker_script} "
            f"--fasta_paths={fasta_path} "
            f"--max_template_date={max_template_date} "
            f"--data_dir={self.data_dir} "
            f"--output_dir={self.output_dir}"
        ]

        logger.info("Running AlphaFold with command: %s", " ".join(cmd))

        # Execute AlphaFold command
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            shell=True,
            executable="/bin/bash",
            check=False,
        )

        # Log output streams
        if process.stdout:
            logger.info("AlphaFold stdout:\n%s", process.stdout)
        if process.stderr:
            logger.error("AlphaFold stderr:\n%s", process.stderr)

        if process.returncode != 0:
            error_msg = f"AlphaFold execution failed: {process.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return {}  # TODO: Implement return value with confidence score and structure path
