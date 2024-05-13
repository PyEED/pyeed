import toml
from sdrdm_database import DBConnector

from pyEED.core.proteininfo import ProteinInfo

# Get the protein sequence from NCBI
aldolase = ProteinInfo.from_ncbi("NP_001287541.1")

# Start a blast search with the protein sequence of tem1 as query
print("Starting blast search...")
blast_results = aldolase.pblast(n_hits=10)

# Establish a connection to the database
print("Connecting to database...")
db = DBConnector(**toml.load(open(".github/scripts/env.toml")))

# Insert all blast results into the database
db.insert(*blast_results, verbose=True)
db.connection.table("ProteinInfo")
