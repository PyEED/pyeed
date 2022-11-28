```mermaid
classDiagram
    
    class ProteinSequence {
        +string name*
        +string amino_acid_sequence*
        +string nr_id
        +string uniprot_id
    }
    
    class Organism {
        +string ncbi_taxonomy_id*
    }
    
```