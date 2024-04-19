---
prefixes:
  EDAM: http://edamontology.org/
---

# Sequence data model

## Macromolecules

### SequenceRecord (EDAM:data_0848)

A molecular sequence and associated annotation data.

- uri
  - Type: string
  - Description: URI of the sequence.
  - Term: EDAM:data_0848
- accession_id
  - Type: string
  - Description: Accession ID of the sequence.
  - Term: EDAM:data_2091
- name
  - Type: string
  - Description: Arbtrary name of the sequence.
  - Term: EDAM:data_1009
- organism
  - Type: Organism
  - Description: The organism from which the sequence was obtained.

### ProteinRecord[_SequenceRecord_] (EDAM:data_2886)

A protein sequence and associated metadata.

- **sequence**
  - Type: string
  - Description: Amino acid sequence of the protein.
  - Term: EDAM:data_2976
- regions
  - Type: Region[]
  - Description: Defines regions within the protein sequence
- sites
  - Type: Site[]
  - Description: Defines sites within the protein sequence
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
- pdb_id
  - Type: string
  - Description: Protein Data Bank (PDB) identifier.
  - Term: EDAM:data_1127

### DNARecord[_SequenceRecord_] (EDAM:data_2887)

Description of a nucleotide sequence 🧬.

- **sequence**
  - Type: string
  - Description: Nucleotide sequence.
  - Term: EDAM:data_3494
- regions
  - Type: Region[]
  - Description: Defines regions within the nucleotide sequence.
- sites
  - Type: Site[]
  - Description: Defines sites within the nucleotide sequence that code for the protein sequence.
- gc_content
  - Type: float
  - Description: GC content of the sequence.

### AbstractAnnotation

- uri
  - Type: string
  - Description: URI of the annotation.
- accession_id
  - Type: string
  - Description: Accession ID of the annotation.
  - Term: EDAM:data_3034
- name
  - Type: string
  - Description: A name of a sequence feature, e.g. the name of a feature

### Site[_AbstractAnnotation_]

Position(s) constituting a site within a sequence.

- positions
  - Type: integer[]
  - Description: Position of the site(s) within the sequence.

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

## Enumerations

### AnnotationType

Ontology terms for different sections of a sequence.

```python
ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
DOMAIN = "http://semanticscience.org/resource/SIO_001379"
FAMILY = "http://semanticscience.org/resource/SIO_001380"
MOTIVE = "http://semanticscience.org/resource/SIO_000131"
```

### SequenceType

Ontology terms for different types of sequences.

```python
DNA = "http://semanticscience.org/resource/SIO_010018"
PROTEIN = "http://semanticscience.org/resource/SIO_010015"
```
