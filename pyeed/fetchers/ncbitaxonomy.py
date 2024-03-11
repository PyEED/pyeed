from typing import Generator, List
from Bio import Entrez
from tqdm import tqdm

from pyeed.core.organism import Organism
from pyeed.fetchers.abstractfetcher import AbstractFetcher
from pyeed.fetchers.abstractfetcher import LOGGER


class NCBITaxonomyFetcher(AbstractFetcher):
    """
    The `NCBITaxonomyFetcher` class is a subclass of the `AbstractFetcher` class and is used
    to fetch taxonomy data from the NCBI database.
    It provides methods to fetch taxonomy data for a single ID or multiple IDs in chunks.
    The fetched data is then mapped to an instance of the `Organism` class.

    Example Usage:
        ```py
        # Create an instance of NCBITaxonomyParser
        fetcher = NCBITaxonomyFetcher(foreign_id=[9606, 10090], email="example@gmail.com", api_key="API_KEY")

        # Fetch taxonomy data for multiple IDs
        results = fetcher.fetch(Organism)
        ```
    """

    def __init__(
        self, foreign_id: int | List[int], email: str = None, api_key: str = None
    ):
        super().__init__(foreign_id)

        self.api_key: str = api_key
        if email is None:
            self.email: str = self.get_substitute_email()
        self.taxonomy_dicts: List[dict] = None

    def fetch(self, cls: "Organism"):
        """
        Fetches taxonomy data from NCBI and returns a list of instances of the 'Organism' class.
        """
        tax_dicts = self.get()
        return self.map(tax_dicts, cls)

    def get(self) -> List[dict]:
        """
        Fetches taxonomy data from NCBI and returns a list of dictionaries of the results.
        """

        if isinstance(self.foreign_id, list):
            self.taxonomy_dicts = self.get_multiple_ids()
            return self.taxonomy_dicts
        else:
            self.taxonomy_dicts = self.get_single_id()
            return self.taxonomy_dicts

    def make_request(self, request_string: str) -> Generator:
        """
        Makes a request to the NCBI taxonomy database and returns the results.
        """

        Entrez.email = self.email
        Entrez.api_key = self.api_key

        with Entrez.efetch(
            db="taxonomy",
            id=request_string,
            retmode="xml",
            api_key=self.api_key,
        ) as handle:
            return Entrez.read(handle)

    def get_single_id(self):
        """Gets a single texonomy entry from NCBI."""
        return self.make_request(self.foreign_id)

    def get_multiple_ids(self):
        """
        Gets multiple taxonomy entries from NCBI. Requests are made in chunks, making the
        process more reliable.
        """

        request_chunks = tqdm(
            self.make_chunks(self.foreign_id, 100),
            desc=f"f⬇️ Fetching {len(self.foreign_id)} taxonomy entries...",
        )

        results = []
        for chunk in request_chunks:
            request_string = ",".join(map(str, chunk))

            results.extend(self.make_request(request_string))

        return results

    @staticmethod
    def map(taxonomy_dicts: List[dict], cls: "Organism") -> List["Organism"]:

        organisms = []
        for taxonomy_dict in taxonomy_dicts:
            tax_id = taxonomy_dict.get("TaxId")
            organism = cls(taxonomy_id=tax_id)

            organism.name = taxonomy_dict.get("ScientificName")
            organism.species = taxonomy_dict.get("ScientificName")

            lineage = taxonomy_dict.get("LineageEx")

            if not lineage:
                LOGGER.debug(f"No lineage found for {tax_id}: {taxonomy_dict}")
                return organism

            for tax_rank in lineage:
                if tax_rank.get("Rank") == "superkingdom":
                    organism.domain = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "phylum":
                    organism.phylum = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "class":
                    organism.tax_class = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "order":
                    organism.order = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "family":
                    organism.family = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "genus":
                    organism.genus = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "species":
                    organism.species = tax_rank.get("ScientificName")
                elif tax_rank.get("Rank") == "kingdom":
                    organism.kingdom = tax_rank.get("ScientificName")
                else:
                    continue

            organisms.append(organism)

        return organisms


if __name__ == "__main__":
    from pyeed.core import Organism

    single_tax_id = "9606"
    multiple_tax_ids = [9606, "10090", 10116]

    mul = NCBITaxonomyFetcher(multiple_tax_ids).fetch(Organism)

    other = NCBITaxonomyFetcher(single_tax_id).fetch(Organism)
