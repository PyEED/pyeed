from os import name
import unittest
from pyEED.core.organism import Organism
from pyEED.core.proteininfo import ProteinInfo


class TestProteinInfo(unittest.TestCase):
    def setUp(self):
        organism = Organism(
            name="test_organism",
            taxonomy_id=12345,
            domain="test_domain",
            kingdom="test_kingdom",
            phylum="test_phylum",
            tax_class="test_tax_class",
            order="test_order",
            family="test_family",
            genus="test_genus",
            species="test_species",
        )

        self.proteininfo = ProteinInfo(
            name="test_protein",
            source_id="VD12345.34",
            sequence="qweasdyxc",
            ec_number="1.2.3.13",
            mol_weight=123.45,
            organism=organism,
        )

    def test_content(self):
        self.assertEqual(self.proteininfo.name, "test_protein")
        self.assertEqual(self.proteininfo.source_id, "VD12345.34")
        self.assertEqual(self.proteininfo.sequence, "qweasdyxc")
        self.assertEqual(self.proteininfo.ec_number, "1.2.3.13")
        self.assertEqual(self.proteininfo.mol_weight, 123.45)


if __name__ == "__main__":
    unittest.main()
