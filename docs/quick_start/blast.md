# Using BLAST

## Using NCBI BLAST
NCBI offers a web interface for blasting. With PyEED this can be programmatically accessed. A BLAST search can be initiated by calling the `ncbi_blast()` method on a `ProteinRecord` object. The method returns the found sequences as a list of `ProteinRecord` objects. As additional parameters, the `ncbi_blast()` method accepts the following arguments:

- n_hits (int): The number of hits to return.
- e_value (float): The e-value threshold for the search.
- db (str): The database to search in. The default is `swissprot`.
- matrix (str): The matrix to use for the search. The default is `BLOSUM62`.
- identity (float): The minimum identity percentage for the search. The default is `0.0`.

``` py
from pyeed.core import ProteinRecord

# Create a ProteinInfo object
protein = ProteinRecord.get_id("UCS38941.1")

# Perform a BLAST search
blast_results = protein.ncbi_blast()
```
!!! info "NCBI BLAST performance"

    Due to server-side limitations of NCBI, BLAST searches might be slowed down or even be blocked, if multiple searches are performed in a short time.
