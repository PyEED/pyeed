import json
import logging
import os
import subprocess
import sys

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import RedirectResponse

app = FastAPI()

# Configure logging to output to stdout without buffering
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def to_fasta(seq: str) -> str:
    return f">query_sequence\n{seq}"


@app.get("/")
async def read_root() -> None:
    logger.debug("Entering root endpoint")
    return RedirectResponse(url="/docs")  # type: ignore


# endpoiunt to upload fasta file (query)
@app.post("/upload_fasta")
async def upload_fasta(fasta_string: str) -> None:
    logger.debug("Entering upload_fasta endpoint")
    logger.debug(f"Received request data: {fasta_string}")

    # upload file to /usr/local/bin/data#
    with open("/usr/local/bin/data/query.fasta", "w") as file:
        file.write(fasta_string)

    logger.debug("Fasta file uploaded successfully")


# read fasta file
@app.get("/read_fasta")
async def read_fasta() -> str:
    logger.debug("Entering read_fasta endpoint")
    with open("/usr/local/bin/data/query.fasta", "r") as file:
        return file.read()


@app.get("/blastp_help")
def blastp_help() -> str:
    logger.debug("Entering /blastp_help endpoint")

    command = ["blastp", "-help"]
    logger.debug(f"Running command: {command}")

    try:
        result = subprocess.run(command, capture_output=True, text=True)

        # Return the help text
        return result.stdout
    except subprocess.CalledProcessError as e:
        # Log and raise an HTTP exception if the subprocess fails
        logger.error(f"blastp help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


@app.get("/blastn_help")
def blastn_help() -> str:
    logger.debug("Entering /blastn_help endpoint")

    command = ["blastn", "-help"]
    logger.debug(f"Running command: {command}")

    try:
        result = subprocess.run(command, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"blastn help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


@app.post("/blast")
async def run_blast(request: Request) -> dict[str, str]:
    """Run BLAST search with provided parameters."""
    try:
        data = await request.json()
        logger.debug(f"Received request data: {data}")

        mode = data["mode"]
        sequence = data["sequence"]
        db_path = data["db_path"]
        db_name = data["db_name"]
        evalue = float(data["evalue"])
        max_target_seqs = int(data["max_target_seqs"])
        num_threads = int(data["num_threads"])

        query_path = "/usr/local/bin/data/query.fasta"
        result_path = "/usr/local/bin/data/result.out"

        # Create FASTA file
        with open(query_path, "w") as file:
            file.write(to_fasta(sequence))

        # Run BLAST
        command = [
            mode,
            "-query",
            query_path,
            "-db",
            f"{db_path}/{db_name}",
            "-evalue",
            str(evalue),
            "-outfmt",
            "0",
            "-num_threads",
            str(num_threads),
            "-out",
            result_path,
            "-max_target_seqs",
            str(max_target_seqs),
        ]

        logger.debug(f"Running command: {command}")
        subprocess.run(command, capture_output=True, check=True, text=True)

        # Read results
        with open(result_path, "r") as file:
            result_data = file.read()

        # Cleanup
        os.remove(query_path)
        os.remove(result_path)

        return {"result": result_data}

    except Exception as e:
        logger.error(f"Error running BLAST: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/blastn")
async def run_blastn(request: Request) -> dict[str, str]:
    logger.debug("Entering /blastn endpoint")

    try:
        data = await request.json()
        logger.debug(f"Received request data: {data}")

        query_filename = "in.fasta"
        result_filename = "out.out"

        logger.debug("Creating input and output files")
        open(result_filename, "w").close()
        with open(query_filename, "w") as file:
            file.write(to_fasta(data["query"]))

        command = [
            "blastn",
            "-query",
            query_filename,
            "-db",
            data["db"],
            "-evalue",
            str(data["evalue"]),
            "-outfmt",
            str(data["outfmt"]),
            "-num_threads",
            str(data["num_threads"]),
            "-out",
            result_filename,
            "-max_target_seqs",
            str(data["max_target_seqs"]),
        ]
        logger.debug(f"Running command: {command}")

        subprocess.run(command, capture_output=True, check=True, text=True)

        logger.debug("Reading results file")
        with open(result_filename, "r") as file:
            result_data = file.read()

        return {"result": result_data}

    except json.JSONDecodeError as e:
        logger.error(f"Failed to parse request JSON: {e}")
        raise HTTPException(status_code=400, detail="Invalid JSON format")
    except subprocess.CalledProcessError as e:
        logger.error(f"BLASTN command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="0.0.0.0", port=6001, reload=True)
