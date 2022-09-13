# PyEED Data Model

Insert a general description ...

### ProteinSequence

- __name*__
  - Type: string
  - Description: Systematic name of the protein.
- __amino_acid_seqeunce*__
  - Type: string
  - Description: The amino acid sequence of the protein sequence object.
- __nr_id__
  - Type: string
  - Description: Identifier for the NCBI NR database
- __uniprot_id__
  - Type: string
  - Description: Identifier for the UniProt database
- __pdb_id__
  - Type: string
  - Description: Identifier for the PDB database
  - Multiple: True
- __organism__
  - Type: Organism
  - Description: Corresponding organism 
- __domain__
  - Type: Domain
  - Description: Domain specification
  - Multiple: True
- __reference_sequence__
  - Type: string
  - Description: Identifier of the sequence used as reference
- __equivalence__
  - Type: Equivalence
  - Description: Positions where the given sequence is equivalent to the reference
  - Multiple: True
- __annotation__
  - Type: Annotation
  - Description: Position-wise annotation of the amino acid seqeunce
  - Multiple: True

### Organism

- __ncbi_taxonomy_id*__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism
  - Min_length: 1
  - Max_length: 10

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
  
### Annotation

- __start_position*__
  - Type: integer
  - Description: Start position of the annotation. A single start position without an end corresponds to a single amino acid
- __end_position__
  - Type: integer
  - Description: Optional end position if the annoation contains more than a single amino acid.
- __function*__
  - Type: string
  - Description: Function that is found in the annotated amino acid or sub-sequence