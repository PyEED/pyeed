# Python Enzyme Engineering Database

!!! warning "API under construction 🏗️"

    The API is currently under construction and is subject to change.

## 🤔 What is PyEED?

`pyeed` is a Python toolkit, that allows easy creation, annotation, and analysis of sequence data. All functionalities are based on a data model, which integrates all information on a given nucleotide or protein sequence in a single object. This allows the bundling of all information on a given sequence, making it available in all creation, annotation, and analysis steps. The entire system is generic and applies to various research scenarios.  
`pyeed` is designed to enable object-oriented programming for bioinformatics. 

## 📝 Data Structure

The data structure of `pyeed` is based on a [data model](https://github.com/PyEED/pyeed/blob/main/specifications/data_model.md)(1), describing the relation between all attributes of a sequence. These attributes include the sequence, the organism, and annotations of the sequence such as sites and regions within the sequence. Furthermore, the information is marked with annotations, marking the origin of the information. 

## 🧰 Tools

`pyeed` implements common tools for clustering, aligning, and visualizing sequences. CLI tools such as `Clustal Omega` are implemented as a Docker Service, allowing easy installation and usage of these tools.