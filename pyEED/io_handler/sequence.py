from Bio import SeqIO

from pyeed.core.region import Region
from pyeed.core.site import Site

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
    """Handel SeqIO entry and return `ProteinSequence`"""

    if "db_source" in entry.annotations:
        if "pdb" in entry.annotations["db_source"]:
            pdb_id = entry.id
        else:
            pdb_id = None

        if "REFSEQ" in entry.annotations["db_source"]:
            reference_sequence = entry.annotations["db_source"].split(" ")[-1]
        else:
            reference_sequence = None

    sites = []
    regions = []
    for feature in entry.features:
        # TODO: assert that only one protein is in the file
        if feature.type == "Protein":
            print("yes")
            if "product" in feature.qualifiers:
                protein_name = feature.qualifiers["product"][0]

            if "calculated_mol_wt" in feature.qualifiers:
                mol_weight = feature.qualifiers["calculated_mol_wt"][0]
            else:
                mol_weight = None

            if "EC_number" in feature.qualifiers:
                ec_number = feature.qualifiers["EC_number"][0]
                print(ec_number)
            else:
                ec_number = None
        else:
            protein_name = entry.description
            mol_weight = None
            ec_number = None

        if feature.type == "source":
            organism = Organism(
                name=feature.qualifiers["organism"][0],
                taxonomy_id=feature.qualifiers["db_xref"][0].split(":")[1],
                id=feature.qualifiers["db_xref"][0],
            )

        if feature.type == "Region":
            regions.append(
                Region(
                    name=feature.qualifiers["region_name"][0],
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    cross_reference=feature.qualifiers["db_xref"][0],
                    note=feature.qualifiers["note"][0],
                )
            )

        if feature.type == "Site":
            site_type = feature.qualifiers["site_type"][0]

            if "note" in feature.qualifiers:
                name = feature.qualifiers["note"][0]
            else:
                name = site_type

            sites.append(
                Site(
                    name=name,
                    positions=[loc for loc in feature.location],
                    cross_reference=feature.qualifiers["db_xref"][0],
                    type=site_type,
                )
            )

    return cls(
        id=entry.id,
        name=protein_name,
        sequence=str(entry.seq),
        reference_sequence=reference_sequence,
        ec_number=ec_number,
        mol_weight=mol_weight,
        organism=organism,
        sites=sites,
        regions=regions,
        pdb_id=pdb_id,
    )
