```mermaid
classDiagram
    SequenceRecord <-- ProteinRecord
    SequenceRecord <-- DNARecord
    AbstractAnnotation <-- Site
    AlignmentData <-- PairwiseAlignment
    AlignmentData <-- ClustalOmegaData
    SequenceRecord *-- Organism
    ProteinRecord *-- Site
    ProteinRecord *-- Region
    DNARecord *-- Site
    DNARecord *-- Region
    Cluster *-- Sequence
    AlignmentData *-- Sequence
    
    class SequenceRecord {
        +string uri
        +string accession_id
        +string name
        +Organism organism
    }
    
    class ProteinRecord {
        +string sequence*
        +Region[0..*] regions
        +Site[0..*] sites
        +Region[0..*] coding_sequence
        +string ec_number
        +float mol_weight
        +string pdb_uri
    }
    
    class DNARecord {
        +string sequence*
        +Region[0..*] regions
        +Site[0..*] sites
        +float gc_content
    }
    
    class AbstractAnnotation {
        +string uri
        +string accession_id
        +string name
    }
    
    class Site {
        +integer[0..*] positions
    }
    
    class Region {
        +integer start
        +integer end
    }
    
    class Organism {
        +string taxonomy_id*
        +string name
        +string domain
        +string kingdom
        +string phylum
        +string tax_class
        +string order
        +string family
        +string genus
        +string species
    }
    
    class BlastData {
        +float identity
        +float evalue
        +int n_hits
        +string substitution_matrix
        +int word_size
        +float gap_open
        +float gap_extend
        +float threshold
        +string db_name
    }
    
    class Cluster {
        +string name
        +Sequence representative
        +Sequence[0..*] members
    }
    
    class Sequence {
        +string sequence_id
        +string sequence
    }
    
    class AlignmentData {
        +string consensus
        +Sequence[0..*] sequences
        +Sequence[0..*] aligned_sequences
    }
    
    class PairwiseAlignment {
        +float score
        +float identity
        +float similarity
        +int gaps
        +int mismatches
    }
    
    class StandardNumbering {
        +str reference_accession_id
        +str numbered_accession_id
        +string[0..*] numbering
    }
    
    class ClustalOmegaData {
        +string version
    }
    
    class AnnotationType {
        << Enumeration >>
        +ACTIVE_SITE
        +BINDING_SITE
        +DOMAIN
        +FAMILY
        +MOTIVE
        +CODING_SEQ
    }
    
    class SequenceType {
        << Enumeration >>
        +DNA
        +PROTEIN
    }
    
```