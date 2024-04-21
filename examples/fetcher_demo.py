import asyncio
import json

from pyeed.core import ProteinRecord
from pyeed.fetch.requester import AsyncRequester
from pyeed.fetch.uniprotmapper import UniprotMapper

# Create a ProteinRecord object
seq = ProteinRecord(sequence="MTEITAAMVKELRES")

# Fetch ernty from Uniprot
accession = "P12345"
requester = AsyncRequester(
    ids=[accession],
    url="https://www.ebi.ac.uk/proteins/api/proteins?format=json&accession=",
    batch_size=1,
    rate_limit=5,
    n_concurrent=20,
)

response = asyncio.run(requester.make_request())
response_dict = json.loads(response[0])[0]

# Map the response to a ProteinRecord object
uniprot_mapper = UniprotMapper()
print(type(response_dict))

with open("response_dict.json", "w") as f:
    json.dump(response_dict, f, indent=4)
protein_info = uniprot_mapper.map_uniprot_data(response_dict)
print(protein_info)
