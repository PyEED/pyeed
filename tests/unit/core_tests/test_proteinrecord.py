# Generated by CodiumAI
import pytest

from pyeed.core.organism import Organism
from pyeed.core.proteinrecord import ProteinRecord
from pyeed.core.region import Region
from pyeed.core.site import Site


class TestProteinRecord:

    # In Docs
    def test_instantiation_from_ncbi(self):
        proteinRecord = ProteinRecord.get_id("UCS38941.1")

    # create a ProteinInfo object with all required fields
    def test_create_protein_info_with_required_fields(self):
        protein_record = ProteinRecord(
            name="MyProtein",
            id="ABC123",
            sequence="MVKLQWERTY",
            organism=Organism(
                name="Homo sapiens",
                taxonomy_id="12345",
                domain="Eukaryota",
                kingdom="Animalia",
                phylum="Chordata",
                tax_class="Mammalia",
                order="Primates",
                family="Hominidae",
                genus="Homo",
                species="sapiens",
            ),
        )

        assert protein_record.name == "MyProtein"
        assert protein_record.id == "ABC123"
        assert protein_record.sequence == "MVKLQWERTY"
        assert protein_record.organism.name == "Homo sapiens"
        assert protein_record.organism.taxonomy_id == 12345
        assert protein_record.organism.domain == "Eukaryota"
        assert protein_record.organism.kingdom == "Animalia"
        assert protein_record.organism.phylum == "Chordata"
        assert protein_record.organism.tax_class == "Mammalia"
        assert protein_record.organism.order == "Primates"
        assert protein_record.organism.family == "Hominidae"
        assert protein_record.organism.genus == "Homo"
        assert protein_record.organism.species == "sapiens"
    
    
    # add a region to ProteinInfo object
    def test_add_region_to_protein_info(self):
        protein_record = ProteinRecord(
            name="MyProtein",
            source_id="ABC123",
            sequence="MVKLQWERTY",
            organism=Organism(
                name="Homo sapiens",
                taxonomy_id="12345",
                domain="Eukaryota",
                kingdom="Animalia",
                phylum="Chordata",
                tax_class="Mammalia",
                order="Primates",
                family="Hominidae",
                genus="Homo",
                species="sapiens",
            ),
        )

        region = protein_record.add_to_regions(
            start=1, 
            end=10,
            id="SpecialRegion",
        )

        assert len(protein_record.regions) == 1
        assert protein_record.regions[0] == region
        assert protein_record.regions[0].id == "SpecialRegion"
        assert protein_record.regions[0].start == 1
        assert protein_record.regions[0].end == 10

    # add a site to ProteinInfo object
    def test_add_site_to_protein_info(self):
        protein_record = ProteinRecord(
            name="MyProtein",
            source_id="ABC123",
            sequence="MVKLQWERTY",
            organism=Organism(
                name="Homo sapiens",
                taxonomy_id="12345",
                domain="Eukaryota",
                kingdom="Animalia",
                phylum="Chordata",
                tax_class="Mammalia",
                order="Primates",
                family="Hominidae",
                genus="Homo",
                species="sapiens",
            ),
        )

        site = protein_record.add_to_sites(
            name="Site1",
            uri="TestUri",
            accession_id='1234',
            id='idTest',
            positions=[1, 2, 3],
        )

        assert len(protein_record.sites) == 1
        assert protein_record.sites[0] == site
        assert protein_record.sites[0].name == "Site1"
        assert protein_record.sites[0].accession_id == '1234'
        assert protein_record.sites[0].positions == [1, 2, 3]
