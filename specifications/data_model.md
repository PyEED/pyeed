# PyEED Data Model

## Macromolecules

### ProteinInfo

Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence.

<details>
  <summary><i>Inspect attributes</i></summary>

- source_id
  - Type: string
  - Description: Identifier of the protein sequence in the source database
- name
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
- coding_sequence_ref
  - Type: [DNARegion](#DNARegion)
  - Description: Defines the coding sequence of the protein
- ec_number
  - Type: string
  - Description: Enzyme Commission number
- mol_weight
  - Type: float
  - Description: Calculated molecular weight of the protein

</details>

### DNAInfo

Description of a nucleotide sequence 🧬

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the nucleotide sequence
- __sequence__
  - Type: string
  - Description: The nucleotide sequence coding for the protein sequence
- organism
  - Type: [Organism](#Organism)
  - Description: Corresponding organism
- regions
  - Type: [DNARegion](#DNARegion)
  - Description: Defines regions within the nucleotide sequence that code for the protein sequence
  - Multiple: True
- source_id
  - Type: string
  - Description: Identifier of the corresponding DNA sequence

</details>

## Annotations

### AbstractRegion

Annotation of a region within a sequence 🗺️

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the annotation
- spans
  - Type: [Span](#Span)
  - Description: Spans of the region. E.g. multiple exons of a gene
  - Multiple: True
- note
  - Type: string
  - Description: Information found in 'note' of an ncbi entry
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

### Span

<details>
  <summary><i>Inspect attributes</i></summary>

- start
  - Type: integer
  - Description: Start position of the span of a region
- end
  - Type: integer
  - Description: End position of the span of a region

</details>

### Site

Annotation of a site within a sequence 📍

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: string
  - Description: Name of the site
- type
  - Type: ProteinSiteType
  - Description: Type of the site
- positions
  - Type: integer
  - Description: Positions of the site
  - Multiple: True
- cross_ref
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
- tax_class
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
GENE = "gene"
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
