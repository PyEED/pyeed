---
prefixes:
  EDAM: http://edamontology.org/
---

# Sequence data model

## Macromolecules

### SequenceRecord (EDAM:data_0848)

A molecular sequence and associated metadata.

- organism
  - Type: Organism
  - Description: The organism from which the sequence was obtained.

### ProteinRecord[_SequenceRecord_] (EDAM:data_2886)

A protein sequence and associated metadata.

- accession_id
  - Type: string
  - Description: Accession ID of the protein sequence.
  - Term: EDAM:data_2091
- name
  - Type: string
  - Description: Name of the protein sequence.
  - Term: EDAM:data_1009
- **sequence**
  - Type: string
  - Description: Amino acid sequence of the protein.
  - Term: EDAM:data_2976
- families
  - Type: Region[]
  - Description: Family of the protein
  - Term: EDAM:data_1131
- domains
  - Type: Region[]
  - Description: Domains of the protein
  - Term: EDAM:data_1468
- sites
  - Type: Position[]
  - Description: Annotations of different sites
- coding_sequence
  - Type: Region[]
  - Description: Defines the coding sequence of the protein
- ec_number
  - Type: string
  - Description: An Enzyme Commission (EC) number of an enzyme.
  - Term: EDAM:data_1011
- mol_weight
  - Type: float
  - Description: Calculated molecular weight of the protein
- pdb_id
  - Type: string
  - Description: Protein Data Bank (PDB) identifier.
- alphafold_id
  - Type: string
  - Description: AlphaFold identifier.

### DNARecord[_SequenceRecord_] (EDAM:data_0848)

Description of a nucleotide sequence 🧬.

- name
  - Type: string
  - Description: Name of the nucleotide sequence.
- **sequence**
  - Type: string
  - Description: Nucleotide sequence.
- regions
  - Type: Region
  - Description: Defines regions within the nucleotide sequence that code for the protein sequence
- ori
  - Type: Position
  - Description: Origin of replication

## Annotations

### Annotation

- accession_id
  - Type: string
  - Description: Accession ID of the annotation.
  - Term: EDAM:data_3034
- name
  - Type: string
  - Description: A name of a sequence feature, e.g. the name of a feature
- publications
  - Type: DOI[]
  - Description: DOI of the publication(s) referenced in the annotation.

### DOI

- doi
  - Type: string
  - Description: Digital Object Identifier (DOI) of a publication.
  - Term: EDAM:data_0849

### Region[_Annotation_]

Annotation of a region within a sequence 🗺️.

- start
  - Type: integer
  - Description: Start position in the amino acid sequence of the annotated region.
- end
  - Type: integer
  - Description: End position position in the amino acid sequence of the annotated region.

### Position

Annotation of a site within a sequence 📍.

- position
  - Type: integer
  - Description: Position in a sequence.

### Organism

Description of an organism 🦠.

- **taxonomy_id**
  - Type: string
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
