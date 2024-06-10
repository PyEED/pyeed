---
hide:
    - navigation
---

# Sequence data model

This page provides comprehensive information about the structure and components of the data model, including detailed descriptions of the types and their properties, information on enumerations, and an overview of the ontologies used and their associated prefixes. Below, you will find a graph that visually represents the overall structure of the data model.

??? quote "Graph"
    ``` mermaid
    flowchart TB
        sequencerecord(SequenceRecord)
        proteinrecord(ProteinRecord)
        dnarecord(DNARecord)
        abstractannotation(AbstractAnnotation)
        site(Site)
        region(Region)
        regionset(RegionSet)
        organism(Organism)
        blastdata(BlastData)
        sequence(Sequence)
        alignmentresult(AlignmentResult)
        pairwisealignmentresult(PairwiseAlignmentResult)
        standardnumbering(StandardNumbering)
        numberedsequence(NumberedSequence)
        clustalomegaresult(ClustalOmegaResult)
        ontology(Ontology)
        annotation(Annotation)
        sequencetype(SequenceType)
        sequencerecord(SequenceRecord) --> organism(Organism)
        sequencerecord(SequenceRecord) --> site(Site)
        sequencerecord(SequenceRecord) --> region(Region)
        sequencerecord(SequenceRecord) --> regionset(RegionSet)
        proteinrecord(ProteinRecord) --> region(Region)
        regionset(RegionSet) --> region(Region)
        alignmentresult(AlignmentResult) --> sequence(Sequence)
        alignmentresult(AlignmentResult) --> sequence(Sequence)
        alignmentresult(AlignmentResult) --> standardnumbering(StandardNumbering)
        standardnumbering(StandardNumbering) --> numberedsequence(NumberedSequence)

        click sequencerecord "#sequencerecord" "Go to SequenceRecord"
        click proteinrecord "#proteinrecord" "Go to ProteinRecord"
        click dnarecord "#dnarecord" "Go to DNARecord"
        click abstractannotation "#abstractannotation" "Go to AbstractAnnotation"
        click site "#site" "Go to Site"
        click region "#region" "Go to Region"
        click regionset "#regionset" "Go to RegionSet"
        click organism "#organism" "Go to Organism"
        click blastdata "#blastdata" "Go to BlastData"
        click sequence "#sequence" "Go to Sequence"
        click alignmentresult "#alignmentresult" "Go to AlignmentResult"
        click pairwisealignmentresult "#pairwisealignmentresult" "Go to PairwiseAlignmentResult"
        click standardnumbering "#standardnumbering" "Go to StandardNumbering"
        click numberedsequence "#numberedsequence" "Go to NumberedSequence"
        click clustalomegaresult "#clustalomegaresult" "Go to ClustalOmegaResult"
        click ontology "#ontology" "Go to Ontology"
        click annotation "#annotation" "Go to Annotation"
        click sequencetype "#sequencetype" "Go to SequenceType"
    ```


## Ontologies
- [sio](http://semanticscience.org/resource/)
- [edam](http://edamontology.org/)


## Types


### SequenceRecord
A molecular sequence and associated annotation data.

__sequence__* `string`

- The letter sequence of the macromolecule.


__id__ `string`

- Unique identifier of the sequence.


__name__ `string`

- Arbitrary name of the sequence.


__organism__ [`Organism`](#organism)

- The organism from which the sequence was obtained.


__seq_length__ `integer`

- Length of the sequence.


__sites__ [`list[Site]`](#site)

- Defines sites within the nucleotide sequence.


__regions__ [`list[Region]`](#region)

- Defines regions within the nucleotide sequence.


__region_sets__ [`list[RegionSet]`](#regionset)

- Multiple regions forming a higher order structure or feature of a sequence.


------

### ProteinRecord
[SequenceRecord]A protein sequence and associated metadata.

__structure_id__ `string`

- Protein Data Bank (PDB) identifier.


__coding_sequence__ [`list[Region]`](#region)

- Defines the coding sequence of the protein


__ec_number__ `string`

- An Enzyme Commission (EC) number of an enzyme.


__mol_weight__ `float`

- Calculated molecular weight of the protein based on the sequence.


------

### DNARecord
[SequenceRecord]A nucleic acid sequence and associated metadata 🧬

__gc_content__ `float`

- GC content of the sequence.


------

### AbstractAnnotation


__url__ `string`

- URI of the annotation.


__accession_id__ `string`

- Accession ID of the annotation.


__name__ `string`

- A name of a sequence feature, e.g. the name of a feature


------

### Site
[AbstractAnnotation]Position(s) constituting a site within a sequence.

__positions__ `list[integer]`

- Position of the site(s) within the sequence.


------

### Region
[AbstractAnnotation]Regional annotation of a feature within a sequence. The direction of the region is defined by the start and end positions.

__start__ `integer`

- Start position of the site.


__end__ `integer`

- End position of the site.


------

### RegionSet
A set of regions forming a higher order structure. For example, a set of exons in a gene, or a set of secondary structures forming a super-secondary structure.

__regions__ [`list[Region]`](#region)

- Regions of the cluster.


------

### Organism
Description of an organism 🦠.

__taxonomy_id__* `integer`

- A stable unique identifier for each taxon (for a species, a family, an order, or any other group in the NCBI taxonomy database.


__name__ `string`

- The name of an organism (or group of organisms).


__domain__ `string`

- Domain of the organism


__kingdom__ `string`

- Kingdom of the organism


__phylum__ `string`

- Phylum of the organism


__tax_class__ `string`

- Class of the organism


__order__ `string`

- Order of the organism


__family__ `string`

- The name of a family of organism.


__genus__ `string`

- The name of a genus of organism.


__species__ `string`

- The name of a species (typically a taxonomic group) of organism.


------

### BlastData


__identity__ `float`

- Minimum identity to safe hits.

- `Default`: 0.0

__evalue__ `float`

- Expectation value (E) to safe hits.

- `Default`: 10.0

__n_hits__ `integer`

- Number of hits to return.

- `Default`: 100

__substitution_matrix__ `string`

- Substitution matrix to use.

- `Default`: ""blosum62""

__word_size__ `integer`

- Word size of the initial match.

- `Default`: 3- `Inclusivminimum`: 2- `Inclusivemaximum`: 7

__gap_open__ `float`

- Gap open penalty.

- `Default`: 11.0

__gap_extend__ `float`

- Gap extend penalty.

- `Default`: 1.0

__threshold__ `float`

- Minimum score to add a word to the BLAST lookup table.

- `Default`: 11

__db_name__ `string`

- Name of the database to search.


------

### Sequence


__sequence_id__ `string`

- Identifier of the sequence in the source database


__sequence__ `string`

- Molecular sequence.


------

### AlignmentResult


__consensus__ `string`

- Consensus sequence of the alignment.


__sequences__ [`list[Sequence]`](#sequence)

- Sequences of the alignment.


__aligned_sequences__ [`list[Sequence]`](#sequence)

- Aligned sequences as a result of the alignment.


__standard_numbering__ [`StandardNumbering`](#standardnumbering)

- Standard numbering of the aligned sequences.


------

### PairwiseAlignmentResult
[AlignmentResult]

__score__ `float`

- Alignment score


__identity__ `float`

- Ratio of identical residues in the alignment


__similarity__ `float`

- Ratio of similar residues in the alignment


__gaps__ `integer`

- Number of gaps in the alignment


__mismatches__ `integer`

- Number of mismatches in the alignment


------

### StandardNumbering


__reference_id__ `string`

- Standard numbering of the reference sequence


__numberd_sequences__ [`list[NumberedSequence]`](#numberedsequence)

- Numbered sequence of the aligned sequence


------

### NumberedSequence


__numbered_id__ `string`

- Identifier of the numbered sequence


__numbering__ `list[string]`

- Standard numbering of the aligned sequence


------

### ClustalOmegaResult
[AlignmentResult]

__version__ `string`

- Version of the Clustal Omega software


## Enumerations

### Ontology

| Alias | Value |
|-------|-------|
| `ECO` | https://www.evidenceontology.org/term/ |
| `GO` | https://amigo.geneontology.org/amigo/term/ |
| `SIO` | http://semanticscience.org/resource/ |

### Annotation

| Alias | Value |
|-------|-------|
| `ACTIVE_SITE` | http://semanticscience.org/resource/SIO_010041 |
| `ALLOSTERIC_SITE` | http://semanticscience.org/resource/SIO_010050 |
| `ALPHAHELIX` | http://semanticscience.org/resource/SIO_010468 |
| `BETASTRAND` | http://semanticscience.org/resource/SIO_010469 |
| `BINDING_SITE` | http://semanticscience.org/resource/SIO_010040 |
| `CODING_SEQ` | http://semanticscience.org/resource/SIO_001276 |
| `DOMAIN` | http://semanticscience.org/resource/SIO_001379 |
| `FAMILY` | http://semanticscience.org/resource/SIO_001380 |
| `MOTIVE` | http://semanticscience.org/resource/SIO_000131 |

### SequenceType

| Alias | Value |
|-------|-------|
| `DNA` | http://semanticscience.org/resource/SIO_010018 |
| `PROTEIN` | http://semanticscience.org/resource/SIO_010015 |