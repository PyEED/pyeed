# PyEED Data Model 𝌭 

The PyEED data model provides a object data structure for protein sequences and their annotations.

## Objects

### ProteinSequence

Description of a protein sequence and its annotations.

<details>
  <summary><i>Inspect attributes</i></summary>

- __name__
  - Type: string
  - Description: Name of the protein
- __sequence__
  - Type: string
  - Description: Amino acid sequence
- __organism__
  - Type: [Organism](#Organism)
  - Description: Corresponding organism
- regions
  - Type: [Annotation](#Annotation)
  - Description: Domains of the protein
  - Multiple: True
- sites
  - Type: [Annotation](#Annotation)
  - Description: Annotations of different sites
  - Multiple: True
- cds
  - Type: 
- nr_id
  - Type: string
  - Description: Identifier for the NCBI NR database
- uniprot_id
  - Type: string
  - Description: Identifier for the UniProt database
- pdb_id
  - Type: string
  - Description: Identifier for the PDB database
- reference_sequence
  - Type: string
  - Description: Identifier of the sequence used as reference
- equivalence
  - Type: [Equivalence](#Equivalence)
  - Description: Positions where the given sequence is equivalent to the reference
  - Multiple: True

</details>

### Organism

Description of an organism.

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the organism
- __taxonomy_id__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism

</details>

### Equivalence

<details>
  <summary><i>Inspect attributes</i></summary>

- __reference_position__
  - Type: integer
  - Description: Equivalent position in the reference sequence
- __sequence_position__
  - Type: integer
  - Description: Position that is equivalent to the reference sequence position that is also given

</details>

### Annotation

Annotation of a protein sequence.

<details>
  <summary><i>Inspect attributes</i></summary>

- __start__
  - Type: integer
  - Description: Start position of the annotation. A single start position without an end corresponds to a single amino acid
- __end__
  - Type: integer
  - Description: Optional end position if the annoation contains more than a single amino acid
- note
  - Type: string
  - Description: Information found in 'note' of an ncbi protein sequence entry
- name
  - Type: string
  - Description: Name of the annotation
- cross_reference
  - Type: string
  - Description: Database cross reference

</details>
  
### DNASequence

<details>
  <summary><i>Inspect attributes</i></summary>

- __sequence__
  - Type: string
  - Description: The DNA sequence
- __organism__
  - Description: Corresponding organism
  - Type: [Organism](#Organism)

</details>
