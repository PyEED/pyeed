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



### DNASequence

- __protein_sequence_id*__
  - Type: string
  - Description: Reference to the corresponding protein sequence to which this DNA sequence translates 
