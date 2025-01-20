import json
import logging
import os
import subprocess
import sys

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import RedirectResponse

app = FastAPI()

logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def create_fastas_file_from_seq(query_string: str, filename: str) -> None:
    """Creates a FASTA file from a single string containing FASTA-formatted sequences."""
    logger.debug(f"Creating FASTA file: {filename}")

    def validate_sequence(sequence: str) -> bool:
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*X")
        sequence = sequence.upper().strip().replace("\n", "")
        return all(char in valid_chars for char in sequence)

    try:
        lines = query_string.strip().split("\n")
        multifasta: list[str] = []
        current_header: str | None = None
        current_sequence: list[str] = []

        for line in lines:
            if line.startswith(">"):
                if current_header:
                    sequence = "".join(current_sequence)
                    if not validate_sequence(sequence):
                        raise ValueError(
                            f"Invalid characters in sequence under {current_header}"
                        )
                    multifasta.append(f"{current_header}\n{sequence}")
                current_header = line.strip()
                current_sequence = []
            else:
                current_sequence.append(line.strip())

        if current_header and current_sequence:
            sequence = "".join(current_sequence)
            if not validate_sequence(sequence):
                raise ValueError(
                    f"Invalid characters in sequence under {current_header}"
                )
            multifasta.append(f"{current_header}\n{sequence}")

        with open(filename, "w", encoding="utf-8") as f:
            f.write("\n".join(multifasta) + "\n")

        logger.debug("FASTA file created successfully")

    except Exception as e:
        logger.error(f"Error creating FASTA file: {e}")
        raise


@app.get("/")
async def read_root() -> None:
    logger.debug("Entering root endpoint")
    return RedirectResponse(url="/docs")  # type: ignore


@app.get("/help")
def help() -> str:
    logger.debug("Entering /help endpoint")

    command = ["mmseqs", "-h"]
    logger.debug(f"Running command: {command}")

    try:
        result = subprocess.run(command, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"MMSeqs2 help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


@app.post("/easycluster")
async def easycluster(request: Request) -> str:
    logger.debug("Entering /easycluster endpoint")

    try:
        data = await request.json()
        logger.debug(f"Received request data: {data}")

        BASE_DIR = "/app"
        query_filename = os.path.join(BASE_DIR, "in.fasta")
        result_filename = os.path.join(BASE_DIR, "output")
        tmp_dir = os.path.join(BASE_DIR, "tmp")

        logger.debug("Creating directories and files")
        os.makedirs(tmp_dir, exist_ok=True)
        open(result_filename, "w").close()

        create_fastas_file_from_seq(data["query"], query_filename)

        command = [
            "mmseqs",
            "easy-cluster",
            query_filename,
            result_filename,
            "--min-seq-id",
            str(data["min_seq_id"]),
            "-c",
            str(data["coverage"]),
            "--cov-mode",
            str(data["cov_mode"]),
            tmp_dir,
        ]
        logger.debug(f"Running command: {command}")

        result = subprocess.run(command, capture_output=True, text=True, check=True)
        logger.debug(f"Command output: {result.stdout}")

        logger.debug("Reading result file")
        with open("/app/output_all_seqs.fasta", "r") as file:
            result_data = file.read()

        return result_data

    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse request JSON: {e}")
        raise HTTPException(status_code=400, detail="Invalid JSON format")
    except subprocess.CalledProcessError as e:
        logger.error(f"MMSeqs2 command failed: {e.stderr}")
        raise HTTPException(status_code=500, detail=f"Command failed: {e.stderr}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="0.0.0.0", port=8001, reload=True)
