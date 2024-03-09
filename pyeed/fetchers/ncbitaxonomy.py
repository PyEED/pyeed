from typing import Generator, List
from Bio import Entrez
from tqdm import tqdm

from pyeed.core.organism import Organism
from pyeed.fetchers.abstractfetcher import AbstractFetcher
from pyeed.fetchers.abstractfetcher import LOGGER


class NCBITaxonomyParser(AbstractFetcher):

    def __init__(self, foreign_id: List[str], email: str = None, api_key: str = None):
        super().__init__(foreign_id)
        self.api_key = api_key
        if email is None:
            self.email = self.get_substitute_email()

    def get(self):
        """
        Fetches taxonomy data from NCBI and returns a list of dictionaries of the results.
        """

        if isinstance(self.foreign_id, list):
            return list(self.get_multiple_ids())
        else:
            return list(self.get_single_id())

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
            request_string = ",".join(chunk)

            results.extend(self.make_request(request_string))

        return results

    def map(self, cls: "Organism"):

        tax_id = self.source.get("TaxId")
        organism = cls(taxonomy_id=tax_id)

        organism.name = self.source.get("ScientificName")
        organism.species = self.source.get("ScientificName")

        lineage = self.source.get("LineageEx")

        if not lineage:
            LOGGER.debug(f"No lineage found for {tax_id}: {self.source}")
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

        return organism


if __name__ == "__main__":
    single_tax_id = "9606"
    multiple_tax_ids = ["9606", "10090", "10116"]

    mul_res = NCBITaxonomyParser(multiple_tax_ids).get()

    single_res = NCBITaxonomyParser(single_tax_id).get()
