# Basics

## 🧬 Create a sequence
PyEED treats Protein sequences the same as DNA sequences. The central object of PyEED is the `ProteinInfo` and `DNAInfo` object, which both are an `AbstractSequence`.  
A sequence object can be created by passing a sequence string to the constructor.  

=== "Protein"

    ``` py
    from pyeed.core import ProteinRecord

    protein = ProteinRecord(
        name="My Protein",
        sequence="MTEITAAMVKELREDKAVQLLREKGLGK"
    )
    ```

=== "DNA"

    ``` py
    from pyeed.core import DNARecord

    dna = DNARecord(sequence="ATGCGTACGTCGATCGATCGATCGATCGATCGATCGATCGATCGTAGTC")
    ```


## 🔎 Search for a sequence

Besides adding sequence information manually, PyEED also allows searching for sequences in the NCBI and UniProt databases. Therefore, the `get_id()` method can be used. In addition to the sequence itself, the method also returns the sequence's annotations and maps them to the corresponding attributes of the sequence object.

=== "Protein"

    ``` py
    protein = ProteinRecord.get_id("UCS38941.1")
    ```

=== "DNA"

    ``` py
    dna_record = DNARecord.get_id('AF188200.1')
    ```

## ⬇️ Save a sequence

### To file

The sequence can be stored in a `FASTA`, `JSON`, `YAML`, or `XML` file format. Therefore, the respective method can be used. The methods below are written for protein, but can be used the same for `dna_record`
The file path is passed as an argument to the method.

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

=== "FASTA"

    ``` py
    protein.to_fasta("protein.fasta")
    ```

### To database
Alternatively, sequence data can be stored in a graph database. Therefore, the `to_db()` method can be used.

```py
# Feature is currently implemented
```
