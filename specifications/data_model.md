# PyEED Data Model

## Macromolecules

### AbstractSequence

<details>
  <summary><i>Inspect attributes</i></summary>

- source_id
  - Type: string
  - Description: Identifier of the sequence in the source database
- name
  - Type: string
  - Description: Name of the sequence
- __sequence__
  - Type: string
  - Description: Sequence of the molecule
- organism
  - Type: [Organism](#Organism)
  - Description: Corresponding organism
- citation
  - Type: [Citation](#Citation)
  - Description: Publication of the sequence

</details>


### ProteinInfo[_AbstractSequence_]

Description of a protein sequence. Additionally, the `ProteinSequence` contains annotations for sites and regions of the protein sequence alongside information on the organism. Furthermore, the `ProteinSequence` contains information on the coding sequence of the protein sequence, which allows later retrieval of the corresponding nucleotide sequence.

<details>
  <summary><i>Inspect attributes</i></summary>

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
- substrates
  - Type: [Substrate](#Substrate)
  - Description: Promiscuous substrates of the protein
  - Multiple: True

</details>

### Structure

<details>
  <summary><i>Inspect attributes</i></summary>

- pdb_id
  - Type: string
  - Description: PDB ID of the structure
- alphafold_id
  - Type: string
  - Description: AlphaFold ID of the structure
- method
  - Type: string
  - Description: Method used for structure determination
- resolution
  - Type: float
  - Description: Resolution of the structure in angstrom
- chains
  - Type: [Chain](#Chain)
  - Description: Chains of the structure
  - Multiple: True
- ligands
  - Type: [Ligand](#Ligand)
  - Description: Ligands of the structure
  - Multiple: True
- mutations
  - Type: int
  - Description: Mutations of the structure

</details>


### DNAInfo

Description of a nucleotide sequence 🧬

<details>
  <summary><i>Inspect attributes</i></summary>

- regions
  - Type: [DNARegion](#DNARegion)
  - Description: Defines regions within the nucleotide sequence that code for the protein sequence
  - Multiple: True

</details>

## Annotations

### Citation

Information on publication of the entry 📖

<details>
  <summary><i>Inspect attributes</i></summary>

- doi
  - Type: str
  - Description: DOI of the publication
- pubmed_id
  - Type: str
  - Description: PubMed ID of the publication
- medline_id
  - Type: str
  - Description: Medline ID of the publication
- year
  - Type: int
  - Description: Year of publication
- authors
  - Type: Author
  - Description: Authors of the publication
  - Multiple: True

</details>

### Author

<details>
  <summary><i>Inspect attributes</i></summary>

- given_name
  - Type: str
  - Description: Given name of the author
- family_name
  - Type: str
  - Description: Family name of the author

</details>

### Substrate

Promiscuous substrate of an enzyme 🧪

<details>
  <summary><i>Inspect attributes</i></summary>

- name
  - Type: str
  - Description: Name of the substrate
- inchi
  - Type: str
  - Description: InChI code of the substrate
- smiles
  - Type: str
  - Description: SMILES code of the substrate
- chebi_id
  - Type: str
  - Description: ChEBI ID of the substrate

</details>

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

## Alignments

### Alignment

<details>
  <summary><i>Inspect attributes</i></summary>

- reference_seq
  - Type: [AbstractSequence](#AbstractSequence)
  - Description: Protein sequence used as reference
  - Alias: reference
- query_seqs
  - Type: [AbstractSequence](#AbstractSequence)
  - Description: Protein sequence used as query
  - Multiple: True
- method
  - Type: string
  - Description: Method used for the alignment
- consensus
  - Type: string
  - Description: Consensus sequence of the alignment
- score
  - Type: float
  - Description: Alignment score
- standard_numberings
  - Type: [StandardNumbering](#StandardNumbering)
  - Description: Standard numbering of the aligned sequences
  - Multiple: True
- identity
  - Type: float
  - Description: Ration of identical residues in the alignment
- similarity
  - Type: float
  - Description: Ration of similar residues in the alignment
- gaps
  - Type: int
  - Description: Number of gaps in the alignment
- mismatches
  - Type: int
  - Description: Number of mismatches in the alignment


</details>

### StandardNumbering

<details>
  <summary><i>Inspect attributes</i></summary>

- sequence_id
  - Type: string
  - Description:  Identifier of the aligned sequence
- numbering
  - Type: string
  - Description: Standard numbering of the aligned sequence
  - Multiple: True


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
