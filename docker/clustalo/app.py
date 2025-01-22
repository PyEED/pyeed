import logging
import subprocess
import sys
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import RedirectResponse

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# Initialize FastAPI application
app = FastAPI(title="Clustal Omega API", version="1.0.0")


@app.get("/")
async def root() -> RedirectResponse:
    """Redirect to API documentation."""
    return RedirectResponse(url="/docs")


@app.get("/help")
def get_help() -> str:
    """Retrieve help information from Clustal Omega."""
    command = ["clustalo", "--help"]

    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Clustal Omega help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


class ClustalORunner:
    """Handles Clustal Omega alignment operations."""

    def __init__(self, base_dir: Path = Path("/app/data")):
        self.base_dir = base_dir
        self.input_file = base_dir / "input.fasta"
        self.output_file = base_dir / "output.clustal"
        base_dir.mkdir(parents=True, exist_ok=True)

    async def run_alignment(self, sequences: str) -> str:
        """Run Clustal Omega alignment on input sequences."""
        # Write sequences to file
        self.input_file.write_text(sequences)
        logger.debug(f"Input file content:\n{sequences}")

        # Build command
        cmd = [
            "clustalo",
            "-i",
            str(self.input_file),
            "-o",
            str(self.output_file),
            "--outfmt=fasta",
            "--force",  # Overwrite output file if it exists
        ]

        try:
            logger.debug(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Read and return results
            alignment = self.output_file.read_text()
            logger.debug(f"Alignment result:\n{alignment}")
            self._cleanup()
            return alignment

        except subprocess.CalledProcessError as e:
            self._cleanup()
            raise HTTPException(
                status_code=500, detail=f"Clustal Omega alignment failed: {e.stderr}"
            )

    def _cleanup(self) -> None:
        """Remove files from the service directory."""
        self.input_file.unlink()
        self.output_file.unlink()


@app.post("/align")
async def align(
    file: UploadFile = File(...),
) -> str:
    """Align sequences using Clustal Omega.

    Args:
        file: FASTA file containing sequences to align

    Returns:
        Alignment in Clustal format
    """
    content = await file.read()
    sequences = content.decode()

    runner = ClustalORunner()
    return await runner.run_alignment(sequences)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app:app",
        host="0.0.0.0",
        port=5001,
        reload=True,
    )
