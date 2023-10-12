# PyEED Data Model

PyEED is a Python-encoded data model of an Enzyme Engineering Database. It supports the scalable and reproducible analysis of sequence and structure data of protein families, and makes data and processes findable, accessible, interoperable, and reusable according to the FAIR data principles.

## Objects

### ProteinSequence

- __name__
  - Type: string
  - Description: Systematic name of the protein.
- __amino_acid_sequence__
  - Type: string
  - Description: The amino acid sequence of the protein sequence object.
- nr_id
  - Type: string
  - Description: Identifier for the NCBI NR database
- uniprot_id
  - Type: string
  - Description: Identifier for the UniProt database
- pdb_id
  - Type: string
  - Description: Identifier for the PDB database
- __organism__
  - Description: Corresponding organism
  - Type: Organism
- domains
  - Type: [Domain](#Domain)
  - Description: Domain specification
  - Multiple: True
- reference_sequence
  - Type: string
  - Description: Identifier of the sequence used as reference
- equivalence
  - Type: [Equivalence](#Equivalence)
  - Description: Positions where the given sequence is equivalent to the reference
  - Multiple: True
- annotations
  - Type: [Annotation](#Annotation)
  - Description: Position-wise annotation of the amino acid sequence
  - Multiple: True

### Organism

- name
  - Type: string
  - Description: Name of the organism
- __ncbi_taxonomy_id__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism

### Domain

- __name__
  - Type: string
  - Description: Name of the annotated domain
- __start_position__
  - Type: integer
  - Description: Position in the sequence where the domain starts
- __end_position__
  - Type: integer
  - Description: Position in the sequence where the domain ends

### Equivalence

- __reference_position__
  - Type: integer
  - Description: Equivalent position in the reference sequence
- __sequence_position__
  - Type: integer
  - Description: Position that is equivalent to the reference sequence position that is also given
  
### Annotation

- __start_position__
  - Type: integer
  - Description: Start position of the annotation. A single start position without an end corresponds to a single amino acid
- __end_position__
  - Type: integer
  - Description: Optional end position if the annoation contains more than a single amino acid.
- note
  - Type: string
  - Description: Function that is found in the annotated amino acid or 
- name
  - Type: string
  - Description: Additional note for the annotation
- db_xref
  - Type: string
  - Description: Database cross reference

### DNASequence

- __protein_sequence_id__
  - Type: string
  - Description: Reference to the corresponding protein sequence to which this DNA sequence translates 
