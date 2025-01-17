---
prefixes:
  sio: http://semanticscience.org/resource/
  edam: http://edamontology.org/
---

# Sequence data model

## Macromolecules

### ProteinRecord (sio:SIO_010043)

A protein sequence and associated metadata.

- **id**
  - Type: string
  - Description: Unique identifier of the sequence.
  - Term: sio:SIO_000729
- name
  - Type: string
  - Description: Arbitrary name of the sequence.
  - Term: sio:SIO_000116
- organism
  - Type: Organism
  - Description: The organism from which the sequence was obtained.
  - Term: sio:SIO_010000
- **sequence**
  - Type: string
  - Description: The letter sequence of the macromolecule.
  - Term: sio:SIO_000030
- embedding
  - Type: float[]
  - Description: 1D embedding vector of the protein sequence.
- seq_length
  - Type: integer
  - Description: Length of the sequence.
  - Term: sio:SIO_000041
- nucleotide_id
  - Type: string
  - Description: Identifier of the nucleotide sequence.
- locus_tag
  - Type: string
  - Description: Locus tag of the protein within the nucleotide sequence.
- sites
  - Type: Site[]
  - Description: Defines sites within the nucleotide sequence.
- regions
  - Type: Region[]
  - Description: Defines regions within the nucleotide sequence.
- structure_ids
  - Type: string[]
  - Description: Identifiers of the structures of the protein.
  - Term: sio:SIO_000729
- ec_number
  - Type: string
  - Description: An Enzyme Commission (EC) number of an enzyme.
  - Term: edam:data_1011
- mol_weight
  - Type: float
  - Description: Calculated molecular weight of the protein based on the sequence.
  - Term: edam:data_1505
- annotations
  - Type: Annotation[]
  - Description: Annotations of the protein sequence.
- go_terms
  - Type: string[]
  - Description: Gene Ontology terms associated with the protein.

### DNARecord (sio:SIO_010008)

A nucleic acid sequence and associated metadata 🧬

- **id**
  - Type: string
  - Description: Unique identifier of the sequence.
  - Term: sio:SIO_000729
- name
  - Type: string
  - Description: Arbitrary name of the sequence.
  - Term: sio:SIO_000116
- organism
  - Type: Organism
  - Description: The organism from which the sequence was obtained.
  - Term: sio:SIO_010000
- **sequence**
  - Type: string
  - Description: The letter sequence of the macromolecule.
  - Term: sio:SIO_000030
- seq_length
  - Type: integer
  - Description: Length of the sequence.
  - Term: sio:SIO_000041
- embedding
  - Type: float[]
  - Description: 1D embedding vector of the protein sequence.
- sites
  - Type: Site[]
  - Description: Defines sites within the nucleotide sequence.
- regions
  - Type: Region[]
  - Description: Defines regions within the nucleotide sequence.
- gc_content
  - Type: float
  - Description: GC content of the sequence.
- annotations
  - Type: Annotation[]
  - Description: Annotations of the DNA sequence.
- go_terms
  - Type: string[]
  - Description: Gene Ontology terms associated with the DNA.

### Site (sio:sio:010049)

Position(s) constituting a site within a sequence.

- name
  - Type: string
  - Description: Name of the site.
- **annotation**
  - Type: Annotation
  - Description: Annotation of the site.
- **positions**
  - Type: integer[]
  - Description: Position of the site(s) within the sequence.
  - Term: sio:SIO_000056

### Region ( _AbstractAnnotation_ ) (sio:SIO_000370)

Regional annotation of a feature within a sequence.

- **id**
  - Type: string
  - Description: Unique identifier of the site.
- name
  - Type: string
  - Description: Name of the site.
- **annotation**
  - Type: Annotation
  - Description: Annotation of the site.
- **start**
  - Type: integer
  - Description: Start position of the site.
  - Term: sio:SIO_000943
- **end**
  - Type: integer
  - Description: End position of the site.
  - Term: sio:SIO_000953

### Organism

Description of an organism 🦠.

- **taxonomy_id**
  - Type: integer
  - Description: A stable unique identifier for each taxon for a species, a family, an order, or any other group in the NCBI taxonomy database.
  - Term: edam:data_1179
- name
  - Type: string
  - Description: The name of an organism (or group of organisms).
  - Term: edam:data_2909
- domain
  - Type: string
  - Description: Domain of the organism
- kingdom
  - Type: string
  - Description: Kingdom of the organism
  - Term: edam:data_1044
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
  - Term: edam:data_2732
- genus
  - Type: string
  - Description: The name of a genus of organism.
  - Term: edam:data_1870
- species
  - Type: string
  - Description: The name of a species (typically a taxonomic group) of organism.
  - Term: edam:data_1045

## Sequence Search

### BlastData

- identity
  - Type: float
  - Description: Minimum identity to safe hits.
  - Default: 0.0
- evalue
  - Type: float
  - Description: Expectation value (E) to safe hits.
  - Default: 10.0
- n_hits
  - Type: integer
  - Description: Number of hits to return.
  - Default: 100
- substitution_matrix
  - Type: string
  - Description: Substitution matrix to use.
  - Default: "BLOSUM62"
- word_size
  - Type: integer
  - Description: Word size of the initial match.
  - Default: 3
  - Inclusivminimum: 2
  - Inclusivemaximum: 7
- gap_open
  - Type: float
  - Description: Gap open penalty.
  - Default: 11.0
- gap_extend
  - Type: float
  - Description: Gap extend penalty.
  - Default: 1.0
- threshold
  - Type: float
  - Description: Minimum score to add a word to the BLAST lookup table.
  - Default: 11
- db_name
  - Type: string
  - Description: Name of the database to search.

## Alignments

### Sequence

- sequence_id
  - Type: string
  - Description: Identifier of the sequence in the source database
- sequence
  - Type: string
  - Description: Molecular sequence.

### AlignmentResult

- consensus
  - Type: string
  - Description: Consensus sequence of the alignment.
- sequences
  - Type: Sequence[]
  - Description: Sequences of the alignment.
- aligned_sequences
  - Type: Sequence[]
  - Description: Aligned sequences as a result of the alignment.
- standard_numbering
  - Type: StandardNumbering
  - Description: Standard numbering of the aligned sequences.

### PairwiseAlignmentResult (_AlignmentResult_)

- score
  - Type: float
  - Description: Alignment score
- identity
  - Type: float
  - Description: Ratio of identical residues in the alignment
- similarity
  - Type: float
  - Description: Ratio of similar residues in the alignment
- gaps
  - Type: integer
  - Description: Number of gaps in the alignment
- mismatches
  - Type: integer
  - Description: Number of mismatches in the alignment

### StandardNumbering

- reference_id
  - Type: string
  - Description: Standard numbering of the reference sequence
- numberd_sequences
  - Type: NumberedSequence[]
  - Description: Numbered sequence of the aligned sequence

### NumberedSequence

- numbered_id
  - Type: string
  - Description: Identifier of the numbered sequence
- numbering
  - Type: string[]
  - Description: Standard numbering of the aligned sequence

## Enumerations

### Annotation

Ontology terms for different sections of a sequence.

```python
ACTIVE_SITE = "http://semanticscience.org/resource/SIO_010041"
BINDING_SITE = "http://semanticscience.org/resource/SIO_010040"
ALLOSTERIC_SITE = "http://semanticscience.org/resource/SIO_010050"
DOMAIN = "http://semanticscience.org/resource/SIO_001379"
FAMILY = "http://semanticscience.org/resource/SIO_001380"
MOTIVE = "http://semanticscience.org/resource/SIO_000131"
CODING_SEQ = "http://semanticscience.org/resource/SIO_001276"
ALPHAHELIX = "http://semanticscience.org/resource/SIO_010468"
BETASTRAND = "http://semanticscience.org/resource/SIO_010469"
DNA = "http://semanticscience.org/resource/SIO_010018"
PROTEIN = "http://semanticscience.org/resource/SIO_010015"
```