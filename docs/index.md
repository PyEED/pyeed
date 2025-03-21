# Python Enzyme Engineering Database

!!! warning "API under construction 🏗️"

    The API is currently under construction and is subject to change.

## 🤔 What is pyeed?

`pyeed` is a Python toolkit for creating Protein and or Nucleotide knowledge graphs for bioinformatic analysis. The knowledge graph is based on the pyeed graph model, structuring Protein and Nucleotide sequences, annotations, and metadata in a Neo4j graph database. pyeed enables seamless data integration of various bioinformatic data sources, such as UniProt and NCBI. Besides the graph model, Python provides a set of tools for sequence analysis, such as sequence alignment or calculation of sequence embeddings.

## 📝 pyeed graph model

The pyeed graph model offers a structure for organizing sequences of proteins and nucleotides and their annotations and metadata. Sequence annotations can describe regions or individual sites of a sequence, such as active sites, binding sites, or domains. Furthermore, the graph model contains information on the sequence's source and Gene Ontology terms of the sequence.

![PyEED Graph Model](./figs/pyeed-model.png)
