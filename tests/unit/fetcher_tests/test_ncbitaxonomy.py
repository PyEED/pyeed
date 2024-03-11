import pytest
import json

from pyeed.fetchers.ncbitaxonomy import NCBITaxonomyFetcher
from pyeed.core import Organism

TAXONOMY_ID = 9606


@pytest.fixture
def mock_ncbi_get(mocker):
    with open("tests/fixtures/data/ecoli_taxonomy_response.json") as f:
        mock_response = json.load(f)
    mocker.patch.object(NCBITaxonomyFetcher, "get", return_value=mock_response)
    return mock_response


class TestNCBITaxonomyParser:

    def test_map(self, mock_ncbi_get):
        tax_fetcher = NCBITaxonomyFetcher(TAXONOMY_ID)
        tax_dict = tax_fetcher.get()
        organisms = tax_fetcher.map(tax_dict, Organism)

        assert len(organisms) == 1
        assert organisms[0].taxonomy_id == "9606"
        assert organisms[0].name == "Homo sapiens"
        assert organisms[0].species == "Homo sapiens"
