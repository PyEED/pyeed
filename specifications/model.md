# PyEED Data Model

## Objects

### ProteinSequence

Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence.

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
  - Type: [Region](#Region)
  - Description: Domains of the protein
  - Multiple: True
- sites
  - Type: [Site](#Site)
  - Description: Annotations of different sites
  - Multiple: True
- coding_sequence
  - Type: [NucleotideSequence](#NucleotideSequence)
  - Description: Information about the coding sequence for the protein sequence
- ec_number
  - Type: string
  - Description: Enzyme Commission number
- mol_weight
  - Type: float
  - Description: Calculated molecular weight of the protein
- nr_id
  - Type: string
  - Description: Identifier for the NCBI NR database
- uniprot_id
  - Type: string
  - Description: Identifier for the UniProt database
- pdb_id
  - Type: string
  - Description: Identifier for the PDB database

</details>

### Organism

Description of an organism 🦠.

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the organism
- __taxonomy_id__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism

</details>

### Region

Annotation of a region within a sequence 🗺️.

<details>
  <summary><i>Inspect attributes</i></summary>

- __start__
  - Type: integer
  - Description: Start position of the annotation. A single start position without an end corresponds to a single amino acid
- __end__
  - Type: integer
  - Description: Optional end position if the annotation contains more than a single amino acid
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

### Site

Annotation of a site within a sequence 📍.

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the site
- type
  - Type: string
  - Description: Type of the site
- positions
  - Type: integer
  - Description: Positions of the site
  - Multiple: True
- cross_reference
  - Type: string
  - Description: Database cross reference

</details>


### NucleotideSequence

Description of a nucleotide sequence 🧬.

<details>
  <summary><i>Inspect attributes</i></summary>

- regions
  - Type: [Region](#Region)
  - Description: Defines regions within the nucleotide sequence that code for the protein sequence
  - Multiple: True
- molecule_type
  - Type: string
  - Description: Type of the sequence
- protein_id
  - Type: string
  - Description: Identifier of the corresponding protein sequence
- gene_id
  - Type: string
  - Description: Identifier of the corresponding gene
- sequence
  - Type: string
  - Description: The nucleotide sequence coding for the protein sequence

</details>
