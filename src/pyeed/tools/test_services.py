from pyeed import Pyeed
from pyeed.model import Protein

pyeed = Pyeed(
    uri="bolt://localhost:7688",
    user="neo4j",
    password="12345678",
)

prots = Protein.nodes.all()
