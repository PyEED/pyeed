import os
import subprocess

from fastapi import FastAPI, UploadFile

app = FastAPI()


@app.get("/help")
def clustalo_help():
    print("Received request for help")
    try:
        result = subprocess.run(
            ["clustalo", "--help"],
            capture_output=True,
            text=True,
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        return {"error": "Command failed", "stderr": e.stderr}, 400


@app.post("/align")
def align(file: UploadFile):
    path = f"/app/data/{file.filename}"
    with open(path, "wb") as f:
        f.write(file.file.read())

    try:
        command = ["clustalo", "-i", path, "--outfmt=clu"]
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )
        os.remove(path)

        return result.stdout
    except subprocess.CalledProcessError as e:
        return {"error": "Command failed", "stderr": e.stderr}, 400


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="0.0.0.0", port=5001, reload=True)
