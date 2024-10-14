<div align="center">
<h1 align="center">PyEED

</div>

[![Tests](https://github.com/PyEED/pyeed/actions/workflows/tests.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/tests.yaml)
[![Documentation](https://github.com/PyEED/pyeed/actions/workflows/make_docs.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/make_docs.yaml)

## About 📖
pyeed is a toolkit enabling object-oriented analysis of protein sequences, instead of working with sequences in a file-oriented fashion. This will enable the user to easily access and manipulate sequence information and to perform analyses on the sequence data.  
This library is currently under development and thus the API is subject to change.


## Installation ⚙️

Install `pyeed` by running
```bash
pip install git+https://github.com/PyEED/pyeed.git
```

## Quick start 🚀

### Launch Neo4j database via Docker and mount to a local directory
```bash
docker run --name pyeed-neo4j -p 7474:7474 -p 7687:7687 -v $PWD/data:/data -v $PWD/logs:/logs -v $PWD/import:/var/lib/neo4j/import -v $PWD/plugins:/plugins -e NEO4J_AUTH=neo4j/test -d neo4j
```
