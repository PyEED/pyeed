# pyeed Docker Services

This directory contains Docker configurations for PyEED's bioinformatics services:
- BLAST for sequence similarity searches
- Clustal Omega for multiple sequence alignment
- MMSeqs2 for fast sequence clustering

## Setup

1. Install Docker and Docker Compose on your system
2. Clone this repository
3. Navigate to the docker directory
4. Run the services:
```bash
docker-compose up -d
```

## Service Details

### BLAST Service (port 6001)
The BLAST service requires a local BLAST database. The service mounts the local database directory to the container:

```yaml
volumes:
  - ./blast:/usr/local/bin/data
```

#### Setting up BLAST Database
1. Place your FASTA sequences in `./blast/test_files/`
2. Create BLAST database in `./blast/test_db/`:
```bash
cd blast/test_db
makeblastdb -in ../test_files/protein_sequences.fasta -dbtype prot -out protein_db
```

The database files (*.pdb, *.phr, etc.) will be available to the container through the mounted volume.

### Clustal Omega Service (port 5001)
Clustal Omega service mounts two directories:
```yaml
volumes:
  - ./clustalo:/app/data   # Input/output data
```

This allows:
- Persistent storage of alignment results
- Access to input/output files from host system

### MMSeqs2 Service (port 8001)
MMSeqs2 service mounts its directory for code and data:
```yaml
volumes:
  - ./mmseqs2:/app
```

This enables:
- Code updates without rebuilding container
- Access to clustering results
- Temporary file management

## Volume Mounting Explained

Docker volumes provide a way to persist data and share files between the host system and containers. In the context of PyEED:

1. **Syntax**: `- ./local/path:/container/path`
   - Left side: Path on your host machine
   - Right side: Path inside the container

2. **Use Cases**:
   - Persisting databases (BLAST)
   - Sharing configuration files
   - Accessing results files
   - Development hot-reloading

3. **Benefits**:
   - Data persists between container restarts
   - Easy data backup
   - Direct file access from host
   - Improved development workflow

## Service URLs

After starting the services, they are available at:
- BLAST: http://localhost:6001
- Clustal Omega: http://localhost:5001
- MMSeqs2: http://localhost:8001

Each service provides a Swagger UI at `/docs` for API documentation.

## Development

To modify the services:
1. Edit code in the respective service directories
2. Changes are reflected immediately due to volume mounting
3. For dependency changes, rebuild the containers:
```bash
docker-compose build
docker-compose up -d
```

## Troubleshooting

1. **Volume Permissions**:
   - Ensure proper read/write permissions on mounted directories
   - Fix with: `chmod -R 755 ./service_directory`

2. **Database Access**:
   - Verify database files are in correct location
   - Check file permissions
   - Ensure proper volume mounting

3. **Port Conflicts**:
   - Change ports in docker-compose.yml if needed
   - Check for other services using same ports

## Complete Example


```yaml
version: "2.2"
services:

  clustalo:
    build:
      context: clustalo/
    container_name: clustalo
    ports:
      - "5001:5001"
    volumes:
      - ./clustalo:/app/data

  blast:
    build:
      context: blast/
    container_name: blast
    ports:
      - "6001:6001"
    volumes:
      - ./blast:/usr/local/bin/data

  mmseqs:
    build:
      context: mmseqs2/
    container_name: mmseqs
    ports:
      - "8001:8001"
    volumes:
      - ./mmseqs2:/app
```

To build the services, run:
```bash
docker-compose build
```

To start the services, run:
```bash
docker-compose up -d
```

