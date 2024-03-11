import pytest
import json

from pyeed.fetchers.ncbitaxonomy import NCBITaxonomyParser
from pyeed.core import Organism

TAXONOMY_ID = 9606


@pytest.fixture
def mock_ncbi_get(mocker):
    with open("tests/fixtures/data/ecoli_taxonomy_response.json") as f:
        mock_response = json.load(f)
    mocker.patch.object(NCBITaxonomyParser, "get", return_value=mock_response)
    return mock_response


class TestNCBITaxonomyParser:

    def test_make_chunks_single_input(self):
        single_list = [1]
        tax_fetcher = NCBITaxonomyParser(single_list)

        chunks = tax_fetcher.make_chunks(single_list, 100)

        assert len(chunks) == 1
        assert len(chunks[0]) == 1
        assert chunks[0][0] == 1

    def test_make_chunks(self):
        long_list = list(range(987))
        tax_fetcher = NCBITaxonomyParser(long_list)

        chunks = tax_fetcher.make_chunks(long_list, 100)

        assert len(chunks) == 10
        assert len(chunks[0]) == 100
        assert len(chunks[-1]) == 87
        assert chunks[0][0] == 0
        assert chunks[-1][-1] == 986
        assert chunks[0][-1] == 99

    def test_map(self, mock_ncbi_get):
        tax_fetcher = NCBITaxonomyParser(TAXONOMY_ID)
        tax_fetcher.taxonomy_dicts = tax_fetcher.get()
        organisms = tax_fetcher.map(Organism)

        assert len(organisms) == 1
        assert organisms[0].taxonomy_id == "9606"
        assert organisms[0].name == "Homo sapiens"
        assert organisms[0].species == "Homo sapiens"
