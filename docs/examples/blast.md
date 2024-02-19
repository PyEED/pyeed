# Using BLAST

## Using NCBI BLAST

NCBI serves offer a web interface for blasting. With PyEED this can be programmatically accessed. A BLAST search can be initiated by calling the `ncbi_blast()` method on a `ProteinInfo` object. The method returns the found sequences as a list of `ProteinInfo` objects.

``` py
from pyEED.core import ProteinInfo

# Create a ProteinInfo object
protein = ProteinInfo.from_db("UCS38941.1")

# Perform a BLAST search
blast_results = protein.ncbi_blast()
```
!!! info "NCBI BLAST performance"

    Due to server-side limitations of NCBI, BLAST searches might be slowed down or even be blocked, if multiple searches are performed in a short time.




## Using BLAST with a local database

Building a local BLAST database is a good way to speed up BLAST searches. PyEED allows to perform BLAST searches on local databases. The `local_blast()` method can be called on a `ProteinInfo` object. The method returns the found sequences as a list of `ProteinInfo` objects.

``` py