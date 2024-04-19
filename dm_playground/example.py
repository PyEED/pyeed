from pyeed_ontology.core import ProteinRecord, Region, AnnotationType

protein = ProteinRecord(
    name="my_protein",
    sequence="MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG",
)

protein.add_to_sites(
    positions=[5, 7, 11],
    annotation=AnnotationType.ACTIVE_SITE.value,
    name="my_active_site",
)

with open("protein.json", "w") as f:
    f.write(protein.model_dump_json())
