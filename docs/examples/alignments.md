# Aligning sequences

## Pairwise sequence alignment

PyEED allows to perform pairwise sequence alignments. Therefore, an `PairwiseAligner` is created to perform the alignment.

``` py
from pyEED.core import ProteinInfo
from pyEED.alignment import PairwiseAligner

# Create two ProteinInfo objects
```




## Multiple sequence alignments

Most sequence alignment tools are implemented as command line tools, which need to be set up. PyEED provides a wrapper for common alignment tools, circumventing the need to setup tools manually. In this process, PyEED builds the alignment tools as a DOCKER container, ensuring the tools run on each system.

### Clustal Omega

To run align sequences with Clustal Omega, a `CLustalOmega` object is ...