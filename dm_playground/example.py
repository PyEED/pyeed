from pyeed_ontology.core import ProteinRecord, Region

protein = ProteinRecord(
    name="my_protein",
    sequence="MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG",
)

# Add annotations to the protein
atp_hydrolysis_activity = "http://purl.obolibrary.org/obo/GO_0016887"
protein.annotations_.append(atp_hydrolysis_activity)

# Add a region, annotated as a domain to the protein
domain1 = protein.add_to_domains(
    start=1, end=34, name="domain1", accession_id="PF00001"
)

# Add a region, annotated as a family to the protein
family1 = protein.add_to_families(
    name="family1", accession_id="pfam00001", start=1, end=34
)

print("annotations of the protein: ", protein.annotations_)
print("annotations of the domain", protein.domains[0].annotations_)
print("annotations of the family", protein.families[0].annotations_)
print(protein.accession_id.__annotations__)




# with open("dm_playground/protein.json", "w") as f:
#     f.write(protein.json())
