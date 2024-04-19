# Alignment

## Objects

### AlignmentResult

- input_sequences
  - Type: Sequence[]
  - Description: Sequences of the alignment
- method
  - Type: string
  - Description: Applied alignment method
- consensus
  - Type: string
  - Description: Consensus sequence of the alignment
- aligned_sequences
  - Type: Sequence[]
  - Description: Aligned sequences of the alignment
- standard_numberings
  - Type: StandardNumbering
  - Description: Standard numbering of the aligned sequences
  - Multiple: True

### PairwiseAlignment[_AlignmentData_]

- score
  - Type: float
  - Description: Alignment score
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

### Sequence

- sequence_id
  - Type: string
  - Description: Identifier of the sequence in the source database
- sequence
  - Type: string
  - Description: Sequence of the alignment. Gaps are represented by '-'

### StandardNumbering

- reference_id
  - Type: str
  - Description: Standard numbering of the reference sequence
- numbered_id
  - Type: str
  - Description: Standard numbering of the query sequence
- numbering
  - Type: string
  - Description: Standard numbering of the aligned sequence
  - Multiple: True
