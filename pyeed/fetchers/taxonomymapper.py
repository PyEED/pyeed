import logging
import xml.etree.ElementTree as ET

from pyeed.core.organism import Organism
from pyeed.fetchers.requester import AsyncRequester


LOGGER = logging.getLogger(__name__)


class NCBITaxonomyMapper:

    def __init__(self):
        pass

    def map(self, xml_data: str) -> Organism:
        """
        Maps XML data to an Organism object.

        Parameters:
            xml_data (str): The XML data to be mapped.

        Returns:
            Organism: The mapped Organism object.
        """

        root = ET.fromstring(xml_data)
        taxon_info = root.find("taxon")

        organism = Organism(
            name=taxon_info.get("ScientificName"),
            taxonomy_id=taxon_info.get("taxId"),
            domain=self.get_scientific_name_by_rank(root, "superkingdom"),
            kingdom=self.get_scientific_name_by_rank(root, "kingdom"),
            phylum=self.get_scientific_name_by_rank(root, "phylum"),
            tax_class=self.get_scientific_name_by_rank(root, "class"),
            order=self.get_scientific_name_by_rank(root, "order"),
            family=self.get_scientific_name_by_rank(root, "family"),
            genus=self.get_scientific_name_by_rank(root, "genus"),
            species=taxon_info.get("ScientificName"),
        )

        return organism

    def get_scientific_name_by_rank(self, root: ET, rank):
        """
        Extracts the scientific name of a taxon based on its rank.

        Parameters:
            root (ElementTree): The root element of the XML data.
            rank (str): The rank of the taxon.

        Returns:
            str or None: The scientific name of the taxon if found, None otherwise.
        """

        for taxon in root.iter("taxon"):
            if taxon.get("rank") == rank:
                return taxon.get("scientificName")

        return None
