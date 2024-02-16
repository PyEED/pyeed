# Python Enzyme Engineering Database

!!! warning "API under construction 🏗️"

    The API is currently under construction and is subject to change. Check out the current [development status](##Installation)

## What is PyEED?

PyEED is a Python toolkit, that allows easy creation, annotation, and analysis of custom sequence data. All functionalities are based on a data model, which integrates all information on a given nucleotide or protein sequence in a single object. This allows the bundling of all information on a given sequence, making it available in all creation, annotation, and analysis steps. The entire system is generic and is also capable of modeling different scenarios.

## PyEED data structure

The object structure of PyEED is based on a [data model](https://github.com/PyEED/pyeed/blob/main/specifications/data_model.md)(1), describing the relation between all attributes of a sequence. These attributes include the sequence, the organism, and annotations of the sequence. Furthermore, the information is marked with annotations, marking the origin of the information. 
{ .annotate }

1.   PyEED uses the [sdRDM framework](https://github.com/JR-1991/software-driven-rdm) to define the architecture of its data as a Markdown document. The hierarchical structure defined in the Markdown document is used to generate Python classes, mirroring the structure of the data model. PyEED can thus be used to read and write data from SQL databases and apply its tools to the data.


### Relational data base

Besides defining the data structure, the data model also defines tables and their relations in a relational database. This allows the native transfer of data that is currently worked with in e.g. a JupyterNotebook and a relational database. Thus, the database structures are inferred from the data model in the same way as the Python classes.

<details>
  <summary>Why do we use a database?</summary>
  The objects we are creating in our JupyterNotebook are only stored for the runtime of the program. After we stop our JupyterNotebook the data would be lost. A database is a persistent data storage and lets us structure our data in such a way, that we need less space on our memory. Therefore we create tables based on the datamodel and use relations in relational database systems to represent those. We are modeling normal references as 1:n and data structures such as sets or lists with n:m by using linking tables.
</details>
