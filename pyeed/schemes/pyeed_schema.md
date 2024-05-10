```mermaid
classDiagram
    SequenceRecord <-- ProteinRecord
    SequenceRecord <-- DNARecord
    AbstractAnnotation <-- Site
    AlignmentResult <-- PairwiseAlignmentResult
    AlignmentResult <-- ClustalOmegaResult
    SequenceRecord *-- Organism
    ProteinRecord *-- Site
    ProteinRecord *-- Region
    DNARecord *-- Site
    DNARecord *-- Region
    Cluster *-- Sequence
    AlignmentResult *-- Sequence
    
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
        +integer taxonomy_id*
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
    
    class AlignmentResult {
        +string consensus
        +Sequence[0..*] sequences
        +Sequence[0..*] aligned_sequences
    }
    
    class PairwiseAlignmentResult {
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
    
    class ClustalOmegaResult {
        +string version
    }
    
    class Ontology {
        << Enumeration >>
        +GO
        +SIO
        +ECO
    }
    
    class Annotation {
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