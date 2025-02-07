import logging
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Annotated

from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import RedirectResponse
from pydantic import BaseModel, Field

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# Initialize FastAPI application
app = FastAPI(title="MMSeqs2 API", version="1.0.0")


@app.get("/")
async def root() -> RedirectResponse:
    """Redirect to API documentation."""
    logger.debug("Accessed root endpoint.")
    return RedirectResponse(url="/docs")


@app.get("/help")
def get_help() -> str:
    """Retrieve help information from MMSeqs2."""
    logger.debug("Accessed /help endpoint.")
    command = ["mmseqs", "-h"]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"MMSeqs2 help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


class ClusterParams(BaseModel):
    """Parameters for MMSeqs2 clustering."""

    min_seq_id: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Minimum sequence identity (0-1)"
    )
    coverage: float = Field(
        default=0.8, ge=0.0, le=1.0, description="Minimum coverage (0-1)"
    )
    cov_mode: int = Field(
        default=0,
        ge=0,
        le=2,
        description="Coverage mode (0: bidirectional, 1: query, 2: target)",
    )
    threads: int = Field(
        default=4,
        ge=1,
        description="Number of CPU threads to use",
    )
    cluster_mode: int = Field(
        default=0,
        ge=0,
        le=2,
        description="Cluster mode (0: set-cover, 1: connected-component, 2: greedy)",
    )
    sensitivity: float = Field(
        default=7.5,
        ge=1.0,
        le=9.0,
        description="Sensitivity: 1.0 faster; 7.5 default; 9.0 more sensitive",
    )
    seq_id_mode: int = Field(
        default=0,
        ge=0,
        le=1,
        description="Sequence identity definition (0: alignment length, 1: shorter sequence)",
    )
    rescore_mode: int = Field(
        default=0,
        ge=0,
        le=1,
        description="Rescore overlapping alignments (0: no, 1: yes)",
    )


class MMSeqs2Runner:
    """Handles MMSeqs2 clustering operations."""

    def __init__(self, base_dir: Path = Path("/app")):
        self.base_dir = base_dir
        self.query_file = base_dir / "in.fasta"
        self.result_file = base_dir / "output"
        self.tmp_dir = base_dir / "tmp"
        self.cluster_file = base_dir / "output_cluster.tsv"
        self.rep_seq_file = base_dir / "output_rep_seq.fasta"
        self.all_seqs_file = base_dir / "output_all_seqs.fasta"
        self.tmp_dir.mkdir(exist_ok=True)

    async def run_clustering(self, query: str, params: ClusterParams) -> str:
        """Run MMSeqs2 clustering on input sequences."""
        # Write query to file
        self.query_file.write_text(query)
        logger.debug(f"Query file content:\n{self.query_file.read_text()}")

        # Run clustering
        cmd = [
            "mmseqs",
            "easy-cluster",
            str(self.query_file),
            str(self.result_file),
            str(self.tmp_dir),
            "--min-seq-id",
            str(params.min_seq_id),
            "-c",
            str(params.coverage),
            "--cov-mode",
            str(params.cov_mode),
            "--threads",
            str(params.threads),
            "--cluster-mode",
            str(params.cluster_mode),
            "-s",
            str(params.sensitivity),
            "--seq-id-mode",
            str(params.seq_id_mode),
            "--rescore-mode",
            str(params.rescore_mode),
        ]

        try:
            logger.debug(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.debug(f"Cluster file content:\n{self.cluster_file.read_text()}")
            time.sleep(0.1)
            content = self.cluster_file.read_text()
            logger.debug(f"Cluster file content (raw): {repr(content)}")
            self.remove_files_and_folders()
            return content
        except subprocess.CalledProcessError as e:
            self.remove_files_and_folders()
            raise HTTPException(
                status_code=500, detail=f"MMSeqs2 clustering failed: {e.stderr}"
            )

    def remove_files_and_folders(self) -> None:
        """Deletes the temporary directory and all related files."""
        # Remove the entire temporary directory and its contents
        if self.tmp_dir.exists():
            shutil.rmtree(self.tmp_dir)

        # Remove additional files outside tmp_dir
        for file in [
            self.cluster_file,
            self.query_file,
            self.rep_seq_file,
            self.all_seqs_file,
            self.result_file,
        ]:
            if file.exists():
                file.unlink()


@app.post("/cluster")
async def cluster(
    query: Annotated[str, Body(description="FASTA formatted sequences")],
    params: Annotated[ClusterParams, Body()] = ClusterParams(),
) -> str:
    """Cluster sequences using MMSeqs2.

    Args:
        query: FASTA formatted sequences
        params: MMSeqs2 clustering parameters

    Returns:
        Clustering results in TSV format
    """
    runner = MMSeqs2Runner()
    return await runner.run_clustering(query, params)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app:app",
        host="0.0.0.0",
        port=8001,
        reload=True,
    )
