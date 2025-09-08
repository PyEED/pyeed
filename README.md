<div align="center">
<h1 align="center">pyeed</h2>
<h2 align="center" >Python Enzyme Engineering Database</h2>
</div>


## About 📖
pyeed is a library for creating custom Sequence Databases from primary databases such as UniProt. pyeed provides an extensible data model, structuring infomration on sequences, annotatios and their embeddings.

## Installation ⚙️

Install `pyeed` by running
```bash
pip install git+https://github.com/PyEED/pyeed.git
```

## Usage 🚀

### Configure Neo4j Database
pyeed provides a `GraphDatabase` class for interacting with a Neo4j database.
It automatically loads the database credentials from the environment variables `NEO4J_URI`, `NEO4J_USER`, and `NEO4J_PASSWORD`.
```
NEO4J_URI="bolt://localhost:7687"
NEO4J_USER="neo4j"
NEO4J_PASSWORD="12345678"
```

### Ingest Data
pyeed provides a `ingest_uniprot` function for ingesting data from UniProt.
It takes a list of UniProt accessions and ingests the data into the database.

```python
import asyncio
from pyeed import ingest_uniprot, Database

db = Database()
accessions = ["P04182", "Q6QDP7", "P04182"]

asyncio.run(await ingest_uniprot(db, accessions))
```

### Query the Database
Cypher queries can be executed against the connected Neo4j database using the `query` method:

```python
from pyeed import Database

db = Database()
query = "MATCH (p:Protein) RETURN p.accession_id AS AC_ID"

results = db.query(query)
```

### Calculate Sequence Embeddings

...