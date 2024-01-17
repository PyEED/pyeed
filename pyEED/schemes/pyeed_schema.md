```mermaid
classDiagram
    AbstractSequence <-- ProteinInfo
    AbstractRegion <-- DNARegion
    AbstractRegion <-- ProteinRegion
    AbstractSequence *-- Citation
    AbstractSequence *-- Organism
    ProteinInfo *-- Substrate
    ProteinInfo *-- DNARegion
    ProteinInfo *-- ProteinRegion
    ProteinInfo *-- Site
    DNAInfo *-- DNARegion
    Citation *-- Author
    AbstractRegion *-- Span
    DNARegion *-- DNARegionType
    ProteinRegion *-- ProteinRegionType
    Site *-- ProteinSiteType
    Alignment *-- AbstractSequence
    Alignment *-- StandardNumbering
    StandardNumbering *-- AbstractSequence
    
    class AbstractSequence {
        +string source_id
        +string name
        +string sequence*
        +Organism organism
        +Citation citation
    }
    
    class ProteinInfo {
        +ProteinRegion[0..*] regions
        +Site[0..*] sites
        +DNARegion coding_sequence_ref
        +string ec_number
        +float mol_weight
        +Substrate[0..*] substrates
    }
    
    class DNAInfo {
        +DNARegion[0..*] regions
    }
    
    class Citation {
        +str doi
        +str pubmed_id
        +str medline_id
        +int year
        +Author[0..*] authors
    }
    
    class Author {
        +str given_name
        +str family_name
    }
    
    class Substrate {
        +str name
        +str inchi
        +str smiles
        +str chebi_id
    }
    
    class AbstractRegion {
        +string name
        +Span[0..*] spans
        +string note
        +string cross_reference
    }
    
    class DNARegion {
        +DNARegionType type
    }
    
    class ProteinRegion {
        +ProteinRegionType type
    }
    
    class Span {
        +integer start
        +integer end
    }
    
    class Site {
        +string name
        +ProteinSiteType type
        +integer[0..*] positions
        +string cross_ref
    }
    
    class Organism {
        +string name
        +string taxonomy_id*
        +string domain
        +string kingdom
        +string phylum
        +string tax_class
        +string order
        +string family
        +string genus
        +string species
    }
    
    class Alignment {
        +AbstractSequence reference_seq
        +AbstractSequence[0..*] query_seqs
        +string method
        +string consensus
        +float score
        +StandardNumbering[0..*] standard_numberings
        +float identity
        +float similarity
        +int gaps
        +int mismatches
    }
    
    class StandardNumbering {
        +AbstractSequence sequence_id
        +string[0..*] numbering
    }
    
    class ProteinSiteType {
        << Enumeration >>
        +ACTIVE
        +BINDING
        +METAL_BINDING
        +POST_TRANS_MODIFICATION
        +UNANNOTATED
    }
    
    class DNARegionType {
        << Enumeration >>
        +CODING_SEQUENCE
        +EXON
        +INTRON
        +GENE
        +PROMOTER
        +ENHANCER
        +UNANNOTATED
    }
    
    class ProteinRegionType {
        << Enumeration >>
        +DOMAIN
        +SIGNAL_PEPTIDE
        +TRANSMEMBRANE
        +UNANNOTATED
    }
    
```