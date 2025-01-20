import logging
import os
import subprocess
import sys

from fastapi import FastAPI, HTTPException, UploadFile
from fastapi.responses import RedirectResponse

app = FastAPI()

logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s - %(message)s",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


@app.get("/")
async def read_root() -> None:
    logger.debug("Entering root endpoint, redirecting to /docs")
    return RedirectResponse(url="/docs")  # type: ignore


@app.get("/help")
def clustalo_help() -> str:
    logger.debug("Entering /help endpoint")

    command = ["clustalo", "--help"]
    logger.debug(f"Running command: {command}")

    try:
        result = subprocess.run(command, capture_output=True, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"ClustalO help command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")


@app.post("/align")
async def align(file: UploadFile) -> str:
    logger.debug("Entering /align endpoint")

    try:
        path = f"/app/data/{file.filename}"
        logger.debug(f"Writing uploaded file to: {path}")

        content = await file.read()
        with open(path, "wb") as f:
            f.write(content)

        command = ["clustalo", "-i", path, "--outfmt=clu"]
        logger.debug(f"Running command: {command}")

        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, command, result.stderr
            )

        logger.debug("Cleaning up temporary file")
        os.remove(path)

        return result.stdout

    except subprocess.CalledProcessError as e:
        logger.error(f"ClustalO command failed: {e.stderr}")
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="0.0.0.0", port=5001, reload=True)
