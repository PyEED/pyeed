from fastapi import FastAPI, HTTPException, Request
from starlette.responses import FileResponse
import logging

import subprocess
import os
import shutil

app = FastAPI()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.info("FastAPI server is running...")


def create_fastas_file_from_seq(query_string, filename):
    """
    Creates a FASTA file from a single string containing FASTA-formatted sequences.

    Args:
        query_string (str): String containing FASTA-formatted sequences.
        filename (str): Path to the output FASTA file.

    Raises:
        ValueError: If any sequence contains invalid characters.
    """
    def validate_sequence(sequence: str) -> bool:
        """Validate that a sequence contains only valid amino acid characters."""
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*X")  # Allow amino acids + stop codon (*), unknown (X)
        sequence = sequence.upper().strip().replace("\n", "")  # Remove whitespace and newlines
        return all(char in valid_chars for char in sequence)

    # Split query string into lines
    lines = query_string.strip().split("\n")

    # Parse headers and sequences
    multifasta = []
    current_header = None
    current_sequence = []

    for line in lines:
        if line.startswith(">"):  # Header line
            if current_header:  # Save the previous sequence
                sequence = "".join(current_sequence)
                if not validate_sequence(sequence):
                    raise ValueError(f"Invalid characters in sequence under {current_header}")
                multifasta.append(f"{current_header}\n{sequence}")
            current_header = line.strip()  # Update header
            current_sequence = []  # Reset sequence buffer
        else:  # Sequence line
            current_sequence.append(line.strip())

    # Add the last sequence
    if current_header and current_sequence:
        sequence = "".join(current_sequence)
        if not validate_sequence(sequence):
            raise ValueError(f"Invalid characters in sequence under {current_header}")
        multifasta.append(f"{current_header}\n{sequence}")

    # Write to file
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("\n".join(multifasta) + "\n")  # Ensure newline at end of file

    print(f"FASTA file created: {filename}")
    
@app.get("/")
async def read_root():
    return {"message": "Welcome to the MMSeqs2 API!"}

@app.get("/help")
def help():
    try: 
        results = subprocess.run(
            ["mmseqs", "-h"],
            capture_output=True,
            text=True,
        )
        return {"help": results.stdout}
    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=400, detail=f"Command failed {e.stderr}")

@app.post("/easycluster")
async def easycluster(request: Request):
    data = await request.json()
    logger.info(f"Received request data: {data}")
    
    BASE_DIR = "/app"
    query_filename = os.path.join(BASE_DIR, "in.fasta")
    result_filename = os.path.join(BASE_DIR, "output")
    tmp_dir = os.path.join(BASE_DIR, "tmp")

    os.makedirs(tmp_dir, exist_ok=True)
    open(result_filename, 'w').close()  # Clear or create result file

    # Create the FASTA file from the query string
    create_fastas_file_from_seq(data['query'], query_filename)

    # Run the mmseqs2 command
    command = [
        "mmseqs", 
        "easy-cluster", 
        query_filename, 
        result_filename, 
        '--min-seq-id', str(data['min_seq_id']),
        '-c', str(data['coverage']),
        '--cov-mode', str(data['cov_mode']),
        tmp_dir]
    logger.info(f"Running command: {' '.join(command)}")
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        logger.info(f"Command output: {result.stdout}")

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with return code {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        raise HTTPException(status_code=500, detail=f"Command failed: {e.stderr}")

    with open("/app/output_all_seqs.fasta", 'r') as file:
        logger.info(f"Reading result file: /app/output_all_seqs.fasta")
        result = file.read()
        
    return result

if __name__ == '__main__':
    import uvicorn
    
    uvicorn.run("app:app", host="0.0.0.0", port=8001, reload=True)