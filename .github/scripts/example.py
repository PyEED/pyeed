import toml
from pyEED.core.proteinsequence import ProteinInfo
from sdrdm_database import DBConnector

# Get the protein sequence from NCBI
aldolase = ProteinInfo.from_ncbi("NP_001287541.1")

# Start a blast search with the protein sequence of tem1 as query
blast_results = aldolase.pblast(n_hits=10)

# Establish a connection to the database
db = DBConnector(**toml.load(open(".github/scripts/env.toml")))

# Insert all blast results into the database
db.insert(*blast_results, verbose=True)
db.connection.table("ProteinSequence")

# Lets filter the blast results for a specific organism
target = "Drosophila melanogaster"

# First, join the ProteinSequence table with the ProteinSequence_organism table
prot_seqs = db.connection.table("ProteinSequence")
organisms = db.connection.table("ProteinSequence_organism")
joined = prot_seqs.join(
    organisms,
    prot_seqs.ProteinSequence_id == organisms.ProteinSequence_id,
    rname="organism_{name}",
)

# Next, filter the joined table for the target organism
filtered = joined.filter(joined.organism_name == target)

# Finally, we can get the corresponding ProteinSequence objects
results = db.get("ProteinSequence", filtered)

assert len(results) > 0
