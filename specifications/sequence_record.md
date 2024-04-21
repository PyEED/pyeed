---
prefixes:
  EDAM: http://edamontology.org/
---

# Sequence data model

## Macromolecules

### SequenceRecord (EDAM:data_0849)

A molecular sequence and associated annotation data.

- uri
  - Type: string
  - Description: URI of the sequence.
  - Term: EDAM:data_1047
- accession_id
  - Type: string
  - Description: Accession ID of the sequence.
  - Term: EDAM:data_2091
- name
  - Type: string
  - Description: Arbtrary name of the sequence.
  - Term: EDAM:data_2099
- organism
  - Type: Organism
  - Description: The organism from which the sequence was obtained.
  - Term: EDAM:data_2530

### ProteinRecord[_SequenceRecord_] (EDAM:data_2886)

A protein sequence and associated metadata.

- **sequence**
  - Type: string
  - Description: Amino acid sequence of the protein.
  - Term: EDAM:data_2976
- regions
  - Type: Region[]
  - Description: Defines regions within the protein sequence
  - Term: EDAM:data_1255
- sites
  - Type: Site[]
  - Description: Defines sites within the protein sequence
  - Term: EDAM:data_1255
- coding_sequence
  - Type: Region[]
  - Description: Defines the coding sequence of the protein
  - Term: EDAM:topic_3511
- ec_number
  - Type: string
  - Description: An Enzyme Commission (EC) number of an enzyme.
  - Term: EDAM:data_1011
- mol_weight
  - Type: float
  - Description: Calculated molecular weight of the protein
- pdb_uri
  - Type: string
  - Description: Protein Data Bank (PDB) identifier.
  - Term: EDAM:data_1047

### DNARecord[_SequenceRecord_] (EDAM:data_2887)

A nucleic acid sequence and associated metadata. 🧬.

- **sequence**
  - Type: string
  - Description: Nucleotide sequence.
  - Term: EDAM:data_3494
- regions
  - Type: Region[]
  - Description: Defines regions within the nucleotide sequence.
  - Term: EDAM:data_1255
- sites
  - Type: Site[]
  - Description: Defines sites within the nucleotide sequence.
  - Term: EDAM:data_1255
- gc_content
  - Type: float
  - Description: GC content of the sequence.

### AbstractAnnotation

- uri
  - Type: string
  - Description: URI of the annotation.
  - Term: EDAM:data_1047
- accession_id
  - Type: string
  - Description: Accession ID of the annotation.
  - Term: EDAM:data_2091
- name
  - Type: string
  - Description: A name of a sequence feature, e.g. the name of a feature
  - Term: EDAM:data_2099

### Site[_AbstractAnnotation_]

Position(s) constituting a site within a sequence.

- positions
  - Type: integer[]
  - Description: Position of the site(s) within the sequence.
  - Term: EDAM:data_1016

### Region

Regional annotation of a feature within a sequence. The direction of the region is defined by the start and end positions.

- start
  - Type: integer
  - Description: Start position of the site.
- end
  - Type: integer
  - Description: End position of the site.

### Organism

Description of an organism 🦠.

- **taxonomy_id**
  - Type: integer
  - Description: A stable unique identifier for each taxon (for a species, a family, an order, or any other group in the NCBI taxonomy database.
  - Term: EDAM:data_1179
- name
  - Type: string
  - Description: The name of an organism (or group of organisms).
  - Term: EDAM:data_2909
- domain
  - Type: string
  - Description: Domain of the organism
- kingdom
  - Type: string
  - Description: Kingdom of the organism
  - Term: EDAM:data_1044
- phylum
  - Type: string
  - Description: Phylum of the organism
- tax_class
  - Type: string
  - Description: Class of the organism
- order
  - Type: string
  - Description: Order of the organism
- family
  - Type: string
  - Description: The name of a family of organism.
  - Term: EDAM:data_2732
- genus
  - Type: string
  - Description: The name of a genus of organism.
  - Term: EDAM:data_1870
- species
  - Type: string
  - Description: The name of a species (typically a taxonomic group) of organism.
  - Term: EDAM:data_1045

## Sequence Search

### BlastData

- identity
  - Type: float
  - Description: Minimum identity to safe hits.
  - Default: 0.0
- evalue
  - Type: float
  - Description: Expectation value (E) to safe hits.
  - Default: 10.0
- n_hits
  - Type: int
  - Description: Number of hits to return.
  - Default: 100
- substitution_matrix
  - Type: string
  - Description: Substitution matrix to use.
  - Default: "BLOSUM62"
- word_size
  - Type: int
  - Description: Word size of the initial match.
  - Default: 3
  - Inclusivminimum: 2
  - Inclusivemaximum: 7
- gap_open
  - Type: float
  - Description: Gap open penalty.
  - Default: 11.0
- gap_extend
  - Type: float
  - Description: Gap extend penalty.
  - Default: 1.0
- threshold
  - Type: float
  - Description: Minimum score to add a word to the BLAST lookup table.
  - Default: 11
- db_name
  - Type: string
  - Description: Name of the database to search.

## Clusters

### Cluster

- name
  - Type: string
  - Description: Name of the cluster.
- representative
  - Type: Sequence
  - Description: Identifier of the representative sequence of the cluster.
- members
  - Type: Sequence[]
  - Description: Sequences of the cluster.

## Alignments

### Sequence

- sequence_id
  - Type: string
  - Description: Identifier of the sequence in the source database
- sequence
  - Type: string
  - Description: Molecular sequence.

### AlignmentData

- consensus
  - Type: string
  - Description: Consensus sequence of the alignment.
- sequences
  - Type: Sequence[]
  - Description: Sequences of the alignment.
- aligned_sequences
  - Type: Sequence[]
  - Description: Aligned sequences as a result of the alignment.

### PairwiseAlignment(_AlignmentData_)

- score
  - Type: float
  - Description: Alignment score
- identity
  - Type: float
  - Description: Ration of identical residues in the alignment
- similarity
  - Type: float
  - Description: Ration of similar residues in the alignment
- gaps
  - Type: int
  - Description: Number of gaps in the alignment
- mismatches
  - Type: int
  - Description: Number of mismatches in the alignment

### StandardNumbering

- reference_accession_id
  - Type: str
  - Description: Standard numbering of the reference sequence
- numbered_accession_id
  - Type: str
  - Description: Standard numbering of the query sequence
- numbering
  - Type: string[]
  - Description: Standard numbering of the aligned sequence

### ClustalOmegaData(_AlignmentData_)

- version
  - Type: string
  - Description: Version of the Clustal Omega software

## Enumerations

### Ontology

Ontology endponts for different types of sequences.

```python
GO = "https://amigo.geneontology.org/amigo/term/"
SIO = "http://semanticscience.org/resource/"
ECO = "https://www.evidenceontology.org/term/"
```

### Annotation

Ontology terms for different sections of a sequence.

```python
ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
DOMAIN = "http://semanticscience.org/resource/SIO_001379"
FAMILY = "http://semanticscience.org/resource/SIO_001380"
MOTIVE = "http://semanticscience.org/resource/SIO_000131"
CODING_SEQ = "http://semanticscience.org/resource/SIO_001276"
```

### SequenceType

Ontology terms for different types of sequences.

```python
DNA = "http://semanticscience.org/resource/SIO_010018"
PROTEIN = "http://semanticscience.org/resource/SIO_010015"
```
