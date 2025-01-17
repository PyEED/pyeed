from fastapi import FastAPI, HTTPException, Request, UploadFile
import subprocess
import os

app = FastAPI()

def create_fasta_files_from_seq(seq, filename):
    with open(filename, 'w') as file: 
        file.write(f">seq\n{seq}\n")

@app.get("/")
async def read_root(): 
    return {"message": "Welcome to the BLAST API"}

@app.get("/blastp_help")
def blastp_help():
    try:
        result = subprocess.run(
            ["blastp", "--help"],
            capture_output=True,
            text=True,
        )
        return {"help": result.stdout}
    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")

@app.get("/blastn_help")
def blastn_help():
    try:
        result = subprocess.run(
            ["blastn", "--help"],
            capture_output=True,
            text=True,
        )
        return {"help": result.stdout}
    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")

@app.post("/blastp")
async def run_blastp(request: Request):
    data = await request.json()
    
    print(f" Received request to run blastp with data: {data}")
    
    query_filename = f"in.fasta"
    result_filename = f"out.out"
    
    # Clear or create result file
    open(result_filename, 'w').close()
    
    # Create the fasta file from the query string
    create_fasta_files_from_seq(data['query'], query_filename)
    
    try:
        command = [
            "blastp",
            '-query', query_filename,
            '-db', data['db'],
            '-evalue', str(data['evalue']),
            '-outfmt', str(data['outfmt']),
            '-num_threads', str(data['num_threads']),
            '-out', result_filename,
            '-max_target_seqs', str(data['max_target_seqs'])
        ]
        
        result = subprocess.run(
            command, 
            capture_output=True, 
            check=True,
            text=True
        )
        
        with open(result_filename, 'r') as file:
            result_data = file.read()
        
        return {"result": result_data}
    
    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")
    

@app.post("/blastn")
async def run_blastn(request: Request):
    data = await request.json()
    
    print(f" Received request to run blastn with data: {data}")
    
    query_filename = f"in.fasta"
    result_filename = f"out.out"
    
    # Clear or create result file
    open(result_filename, 'w').close()
    
    # Create the fasta file from the query string
    create_fasta_files_from_seq(data['query'], query_filename)
    
    try:
        command = [
            "blastn",
            '-query', query_filename,
            '-db', data['db'],
            '-evalue', str(data['evalue']),
            '-outfmt', str(data['outfmt']),
            '-num_threads', str(data['num_threads']),
            '-out', result_filename,
            '-max_target_seqs', str(data['max_target_seqs'])
        ]
        
        result = subprocess.run(
            command, 
            capture_output=True, 
            check=True,
            text=True
        )
        
        with open(result_filename, 'r') as file:
            result_data = file.read()
        
        return {"result": result_data}
    
    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=400, detail=f"Command failed: {e.stderr}")
    
    
if __name__ == "__main__":
    import uvicorn

    uvicorn.run("app:app", host="0.0.0.0", port=6001, reload=True)

