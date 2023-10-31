# PyEED Data Model

## Macro Molecules

### ProteinInfo

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
  - Type: [ProteinRegion](#ProteinRegion)
  - Description: Domains of the protein
  - Multiple: True
- sites
  - Type: [Site](#Site)
  - Description: Annotations of different sites
  - Multiple: True
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

### DNAInfo

Description of a nucleotide sequence 🧬

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the nucleotide sequence
- sequence
  - Type: string
  - Description: The nucleotide sequence coding for the protein sequence
- regions
  - Type: [DNARegion](#DNARegion)
  - Description: Defines regions within the nucleotide sequence that code for the protein sequence
  - Multiple: True
- organism
  - Type: [Organism](#Organism)
  - Description: Corresponding organism
- molecule_type
  - Type: string
  - Description: Type of the sequence
- protein_id
  - Type: string
  - Description: Identifier of the corresponding protein sequence
- gene_id
  - Type: string
  - Description: Identifier of the corresponding gene

</details>

## Annotations

### AbstractRegion

Annotation of a region within a sequence 🗺️

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the annotation
- __start__
  - Type: integer
  - Description: Start position of the annotation. A single start position without an end corresponds to a single amino acid
- __end__
  - Type: integer
  - Description: Optional end position if the annotation contains more than a single amino acid
- note
  - Type: string
  - Description: Information found in 'note' of an ncbi protein sequence entry
- cross_reference
  - Type: string
  - Description: Database cross reference

</details>

### DNARegion[_AbstractRegion_]

- type
  - Type: [DNARegionType](#DNARegionType)
  - Description: Type of the region within the nucleotide sequence

### ProteinRegion[_AbstractRegion_]

- type
  - Type: [ProteinRegionType](#ProteinRegionType)
  - Description: Type of the region within the protein sequence

### Site

Annotation of a site within a sequence 📍

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

### Organism

Description of an organism 🦠

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the organism
- __taxonomy_id__
  - Type: string
  - Description: NCBI Taxonomy ID to identify the organism
- domain
  - Type: string
  - Description: Domain of the organism
- kingdom
  - Type: string
  - Description: Kingdom of the organism
- phylum
  - Type: string
  - Description: Phylum of the organism
- class
  - Type: string
  - Description: Class of the organism
- order
  - Type: string
  - Description: Order of the organism
- family
  - Type: string
  - Description: Family of the organism
- genus
  - Type: string
  - Description: Genus of the organism
- species
  - Type: string
  - Description: Species of the organism

</details>

## Enumerations

Differentiation between binding sites and binding site region. (Binding site, Binding region)
Any site annotations for a DNA sequence?

Annotation as solely a "domain"?

### ProteinSiteType

```python
ACTIVE = "active"
BINDING = "binding"
METAL_BINDING = "metal"
POST_TRANS_MODIFICATION = "post-translational modification"
UNANNOTATED = "unannotated"
```

### DNARegionType

```python
CODING_SEQUENCE = "coding sequence"
EXON = "exon"
INTRON = "intron"
PROMOTER = "promoter"
ENHANCER = "enhancer"
UNANNOTATED = "unannotated"
```

### ProteinRegionType

```python
DOMAIN = "domain"
SIGNAL_PEPTIDE = "signal peptide"
TRANSMEMBRANE = "transmembrane"
UNANNOTATED = "unannotated"
```

### TaxonomicRank

```python
DOMAIN = "domain"
KINGDOM = "kingdom"
PHYLUM = "phylum"
CLASS = "class"
ORDER = "order"
FAMILY = "family"
GENUS = "genus"
SPECIES = "species"
```