# Using BLAST

## Using NCBI BLAST
NCBI offers a web interface for blasting. With PyEED this can be programmatically accessed. A BLAST search can be initiated by calling the `ncbi_blast()` method on a `ProteinInfo` object. The method returns the found sequences as a list of `ProteinInfo` objects.

``` py
from pyEED.core import ProteinInfo

# Create a ProteinInfo object
protein = ProteinInfo.get_id("UCS38941.1")

# Perform a BLAST search
blast_results = protein.ncbi_blastp()
```
!!! info "NCBI BLAST performance"

    Due to server-side limitations of NCBI, BLAST searches might be slowed down or even be blocked, if multiple searches are performed in a short time.




## Using BLAST with a local database

Building a local BLAST database is a good way to speed up BLAST searches. PyEED allows BLAST searches against local databases. The `blastp()` method can be called on a `ProteinInfo` object. The method returns the found sequences as a list of `ProteinInfo` objects.

``` py

    blast_results = protein.blastp(
        db_path="/PATH/TO/LOCAL/BLAST/DB",
        n_hits=200,
        e_value=0.001,
        word_size=3,
    )
```