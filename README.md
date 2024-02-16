<div align="center">
<h1 align="center">PyEED

</div>

[![Generate API](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml)
[![Documentation](https://github.com/PyEED/pyeed/actions/workflows/make_docs.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/make_docs.yaml)

## About 📖
pyEED is a toolkit enabling object-oriented analysis of protein sequences, instead of working with sequences in a file-oriented fashion. This will enable the user to easily access and manipulate sequence information and to perform analyses on the sequence data.  
This library is currently under development and thus the API is subject to change. Check out the [Roadmap](#roadmap-%EF%B8%8F) for more information. If you have feature ideas, or notes issues with pyEED, submit an [issue](https://github.com/PyEED/pyeed/issues).


## Installation ⚙️

Install pyEED by running
```bash
pip install git+https://github.com/PyEED/pyeed.git
```
... or update pyEED to get the latest features.
```bash
pip install --upgrade git+https://github.com/PyEED/pyeed.git
```

## Quick start 🚀

In the following example, information of the [aldolase](https://www.ncbi.nlm.nih.gov/protein/NP_001287541.1/) (*Drosophila melanogaster*) is retrieved from the corresponding GenBank entry. Thereafter, a protein blast search ist started and the found sequence information is fetched and stored as `ProteinSequence` objects.

```python
from pyEED.core import ProteinInfo

# Get a protein entry from NCBI by accession id
aldolase = ProteinInfo.from_ncbi("NP_001287541.1")
print(aldolase)

# Start a blast search with the protein sequence of tem1 as query
blast_results = aldolase.pblast(n_hits=10)

# Get the corresponding nucleotide entry
aldolase_cds = aldolase.get_dna()

# print the protein and coding sequence of the 2nd blast result
print(blast_results[1].sequence)

# print the nucleotide sequence of tem1
print(aldolase_cds)
```

## Documentation 📘

Check out the [documentation](https://pyeed.github.io/pyeed/) for in-depth information on how to setup PyEED, 
use the build-in tools, and store sequence data in databases.  
Documentation is in the making 🐛

## Roadmap 🛣️

- [x] `ProteinSequence` data model: Object-oriented representation of a protein sequence database entry.
- [x] `ProteinSequence` query: get protein sequence by accession id from NCBI database
- [x] Blast search: get protein sequnces as a result from a blast search
- [x] Retrieve corresponding coding sequence
- [x] Storing `ProtenSequence` in a SQL database
- [x] Running pairwise alignments
- [x] Network analysis / visualization
- [ ] Create phylogenetic trees
- [ ] Multi-sequence alignments
- [ ] Representative clustering
