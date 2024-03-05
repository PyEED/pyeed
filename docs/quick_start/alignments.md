# Aligning Sequences

## 🪢 The `Alignment` object

The `Alignment` is the central object of working with sequence alignments. 
Both alignment objects contain the following attributes:  

- `input_sequences`, containing the sequences of the alignment
- `method`, denoting the applied alignment method
- `consensus`, containing the consensus sequence of the alignment
- `aligned_sequences`, containing the aligned sequences as a result of the alignment
- `standard_numberings`, containing the standard numbering of the aligned sequences to a reference sequence

Besides the `Alignment` object, PyEED also provides a `PairwiseAlignment` object, containing the alignment score, identity, similarity, gaps, and mismatches of a pairwise alignment.


Before running the alignment, an `Alignment` needs to be created. This can be done by passing a list of `ProteinInfo` or `DNAInfo` objects to the constructor. The alignment can then be run by calling the `align()` method, passing the alignment method as an argument. The method returns the alignment object, containing the aligned sequences.

``` py
from pyeed.core import ProteinInfo, Alignment

# Get two ProteinInfo objects
tem1 = ProteinInfo.from_ncbi("QGC48744.1")
tem109 = ProteinInfo.from_ncbi("AAT46413.1")

# Create an Alignment
alignment = Alignment(input_sequences=[tem1, tem109])
```

Alternatively, the `from_sequneces()` class method can be used to create an alignment from a list of `ProteinInfo` or `DNAInfo` objects.

``` py
from pyeed.core import ProteinInfo, Alignment

# Get two ProteinInfo objects
tem1 = ProteinInfo.from_ncbi("QGC48744.1")
tem109 = ProteinInfo.from_ncbi("AAT46413.1")
list_of_sequences = [tem1, tem109]

# Create an Alignment
alignment = Alignment.from_sequences(list_of_sequences)
```

### 🔗 Pairwise Alignments

=== "Local Alignment"

    ``` py
    from pyeed.core import ProteinInfo, Alignment
    from pyeed.aligners import PairwiseAligner

    # Get two ProteinInfo objects
    tem1 = ProteinInfo.from_ncbi("QGC48744.1")
    tem109 = ProteinInfo.from_ncbi("AAT46413.1")

    # Create and run alignment
    alignment = PairwiseAlignment([tem1, tem109], aligner=PairwiseAligner, mode="local")
    ```

=== "Global Alignment"
    ``` py
    from pyeed.core import ProteinInfo, PairwiseAlignment
    from pyeed.aligners import PairwiseAligner

    # Get two ProteinInfo objects
    tem1 = ProteinInfo.from_ncbi("QGC48744.1")
    tem109 = ProteinInfo.from_ncbi("AAT46413.1")

    # Create and run alignment
    alignment = PairwiseAlignment([tem1, tem109], aligner=PairwiseAligner, mode="global")
    ```

=== "Multi-Pairwise Alignment"
    If more than two sequences are subjected to a `PairwiseAligner`, all unique sequence pairs are aligned. 
    Instead of an `Alignment` object, a `List[PairwiseAlignment]` object is returned, containing all individual pairwise alignments.

    ``` py
    from pyeed.core import ProteinInfo, Alignment
    from pyeed.aligners import PairwiseAligner

    # Get sequences
    ncbi_accessions = ["QGC48744.1", "AAT46413.1", "AAT46414.1", "AAT46415.1"]
    sequences = ProteinInfo.from_ncbi(ncbi_accessions)

    # Create and run alignment
    alignment = Alignment.from_sequences(sequences, aligner=PairwiseAligner, mode="global)

    ```

### ⛓️ Multiple Sequence Alignments

Most sequence alignment tools are implemented as command line tools, which need to be set up. PyEED implements common alignment tools such as `ClustalOmega` wrapped as Docker containers. As a result, no manual installation of tools is required, since the tools are automatically installed when needed.

=== "ClustalOmega"

    ``` py
    from pyeed.core import ProteinInfo, Alignment
    from pyeed.aligners import ClustalOmega

    # Get sequences
    ncbi_accessions = ["QGC48744.1", "AAT46413.1", "AAT46414.1", "AAT46415.1"]
    sequences = ProteinInfo.from_ncbi(ncbi_accessions)

    # Create and run alignment
    alignment = Alignment.from_sequences(sequences, aligner=ClustalOmega)
    ```

=== "MUSCLE"

    ``` py
    # Not yet implemented
    ```

## 💯 Standard Numbering

Standard numbering is a way to express the alignment of a sequence to a reference sequence in the form of numbers.
The Standard Numbering is created based on aligned sequences in which aligned positions are at the same position within the alignment. If a position is not aligned, and a gap is introduced and denoted with `-`.  
