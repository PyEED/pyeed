<div align="center">
<h1 align="center">pyEED

</div>

[![Generate API](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml)

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

In the following example, information of the [beta-lactamase](https://www.ncbi.nlm.nih.gov/protein/QGC48744.1) (*Escherichia coli*) is retrieved from the corresponding GenBank entry. Thereafter, a protein blast search ist started and the found sequence information is fetched and stored as `ProteinSequence` objects. Finally, the corresponding coding sequence is retrieved for all sequences.

```python
from pyEED.core.proteinsequence import ProteinSequence
from pyEED.ncbi.utils import get_nucleotide_sequences

# Get a protein entry from NCBI by accession id
tem1 = ProteinSequence.from_ncbi("QGC48744.1")

# Start a blast search with the protein sequence of tem1 as query
blast_results = tem1.pblast(n_hits=10)

# Get the corresponding nucleotide sequences of all blast results
get_nucleotide_sequences(blast_results)

# print the protein and coding sequence of the 2nd blast result
print(blast_results[1].sequence)
print(blast_results[1].nucleotide_seq)
```

A proper documentation is in the making 🐛
## Roadmap 🛣️

- [x] `ProteinSequence` data model: Object-oriented representation of a protein sequence database entry.
- [x] `ProteinSequence` query: get protein sequence by accession id from NCBI database
- [x] Blast search: get protein sequnces as a result from a blast search
- [x] Retrieve corresponding coding sequence
- [x] Storing `ProtenSequence` in a SQL database
- [ ] Running alignments
- [ ] Create phylogenetic trees
- [ ] Network analysis / visualization
