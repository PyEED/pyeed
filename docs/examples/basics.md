# Basics

## 🧬 Create a sequence
PyEED treats Protein sequences the same as DNA sequences. The central object of PyEED is the `ProteinInfo` and `DNAInfo` object, which both are an `AbstractSequence`.  
A sequence object can be created by passing a sequence string to the constructor.  

=== "Protein"

    ``` py
    from pyEED.core import ProteinInfo

    protein = ProteinInfo(sequence="MTEITAAMVKELREDKAVQLLREKGLGK")
    ```

=== "DNA"

    ``` py
    from pyEED.core import DNAInfo

    dna = DNAInfo(sequence="ATGCGTACGTCGATCGATCGATCGATCGATCGATCGATCGATCGTAGTC")
    ```


## 🔎 Search for a sequence

Besides adding sequence information manually, PyEED also allows to search for sequences in the NCBI and UniProt databases. Therefore, the `from_db()` method can be used. In addition to the sequence itself, the method also returns the sequence's annotations and maps them to the corresponding attributes of the sequence object.

=== "Protein"

    ``` py
    protein = ProteinInfo.from_db("UCS38941.1")
    ```

=== "DNA"

    ``` py
    dna = DNAInfo.from_db("NC_000913.3")
    ```

Alternatively, the sequence can be initiated from a sequence string, triggering a BLAST search in the NCBI database. If the sequence is found, the sequence object is filled with the corresponding information.

=== "Protein"

    ``` py
    # Not yet implemented
    ```

=== "DNA"

    ``` py
    # Not yet implemented
    ```

## ⬇️ Save a sequence

### To file

The sequence can be stored in a `FASTA`, `JSON`, `YAML`, or `XML`file format. Therefore, the respective method can be used.

=== "FASTA"

    ``` py
    protein.to_fasta("protein.fasta")
    ```

=== "JSON"

    ``` py
    protein.to_json("protein.json")
    ```

=== "YAML"

    ``` py
    protein.to_yaml("protein.yaml")
    ```

=== "XML"

    ``` py
    protein.to_xml("protein.xml")
    ```

### To database
Alternatively, sequence data can be stored in a `PostgreSQL` database. Therefore, the `to_db()` method can be used.

```py
# Feature is currently implemented
```
