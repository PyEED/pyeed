from Bio import SeqIO

from pyeed.core.annotation import Annotation

from ..core.organism import Organism


def SeqIO_to_pyeed(entry: SeqIO):
    """Handel SeqIO entry and return pyeed object."""

    if entry.annotations["molecule_type"] == "protein":
        return _seqio_to_protein_sequence(entry)

    elif entry.annotations["molecule_type"] == "DNA":
        raise NotImplementedError("DNA is not implemented yet.")

    else:
        raise ValueError(
            f"{entry.id} of type {entry.annotations['molecule_type']} is not 'protein' or 'DNA'."
        )


def _seqio_to_protein_sequence(cls, entry: SeqIO):
    id = entry.id
    name = entry.description
    amino_acid_sequence = str(entry.seq)

    if "db_source" in entry.annotations:
        if "pdb" in entry.annotations["db_source"]:
            pdb_id = entry.id
        else:
            pdb_id = None

        if "REFSEQ" in entry.annotations["db_source"]:
            reference_sequence = entry.annotations["db_source"].split(" ")[-1]
        else:
            reference_sequence = None

    annotations = []
    for feature in entry.features:
        if feature.type == "source":
            organism = Organism(
                name=feature.qualifiers["organism"][0],
                ncbi_taxonomy_id=feature.qualifiers["db_xref"][0],
                id=feature.qualifiers["db_xref"][0],
            )

        if feature.type == "Region":
            annotations.append(
                Annotation(
                    name=feature.qualifiers["region_name"][0],
                    start_position=int(feature.location.start),
                    end_position=int(feature.location.end),
                    db_xref=feature.qualifiers["db_xref"][0],
                    note=feature.qualifiers["note"][0],
                )
            )

        if feature.type == "Site":
            annotations.append(
                Annotation(
                    name=feature.qualifiers["site_type"][0],
                    start_position=int(feature.location.start),
                    end_position=int(feature.location.end),
                    db_xref=feature.qualifiers["db_xref"][0],
                )
            )

    return cls(
        id=id,
        name=name,
        amino_acid_sequence=amino_acid_sequence,
        reference_sequence=reference_sequence,
        organism=organism,
        annotations=annotations,
        pdb_id=pdb_id,
    )
