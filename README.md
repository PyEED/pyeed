<div align="center">
<h1 align="center">pyEED

</div>

[![Generate API](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml/badge.svg)](https://github.com/PyEED/pyeed/actions/workflows/generate_api.yaml)


## Installation ⚙️

Install pyEED by running
```bash
pip install git+https://github.com/PyEED/pyeed.git
```
or update pyEED to get the latest features.
```bash
pip install --force-reinstall --no-deps git+https://github.com/PyEED/pyeed.git
```


## 🛣️ Roadmap

- [x] `ProteinSequence` data model: Object-oriented representation of a protein sequence database entry.
- [x] `ProteinSequence` query: get protein sequence by accession id from NCBI database
- [x] Blast search: get protein sequnces as a result from a blast search
- [x] Retrieve corresponding coding sequence
- [ ] Storing `ProtenSequence` in a SQL database
- [ ] Running alignment
- [ ] Network analysis
- [ ] Network visualisation
