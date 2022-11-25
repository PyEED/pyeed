# PyEED Data Model

PyEED is a Python-encoded data model of an Enzyme Engineering Database. It supports the scalable and reproducible analysis of sequence and structure data of protein families, and makes data and processes findable, accessible, interoperable, and reusable according to the FAIR data principles.

### ProteinSequence

- __name*__
  - Type: string
  - Description: Systematic name of the protein.
- __amino_acid_sequence*__
  - Type: string
  - Description: The amino acid sequence of the protein sequence object.
- __nr_id__
  - Type: string
  - Description: Identifier for the NCBI NR database
- __uniprot_id__
  - Type: string
  - Description: Identifier for the UniProt database


### Organism

- __ncbi_taxonomy_id*__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism

### Domain

- __name*__
  - Type: string
  - Description: Name of the annotated domain
- __start_position*__
  - Type: integer
  - Description: Position in the sequence where the domain starts
- __end_position*__
  - Type: integer
  - Description: Position in the sequence where the domain ends

### Equivalence

- __reference_position*__
  - Type: integer
  - Description: Equivalent position in the reference sequence
- __sequence_position*__
  - Type: integer
  - Description: Position that is equivalent to the reference sequence position that is also given
  



