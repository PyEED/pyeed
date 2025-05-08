import json
from collections import defaultdict
from typing import Any, List, Optional

import requests
from bs4 import BeautifulSoup, Tag
from httpx import Response
from loguru import logger
from SPARQLWrapper import JSON, SPARQLWrapper

from pyeed.adapter.primary_db_adapter import PrimaryDBMapper
from pyeed.model import (
    Annotation,
    GOAnnotation,
    Molecule,
    Organism,
    Protein,
    Reaction,
    Region,
    Site,
)


class UniprotToPyeed(PrimaryDBMapper):
    def add_to_db(self, response: Response) -> None:
        # Organism information
        records = self.parse_response(response)

        for record in records:
            taxonomy_id = record["organism"]["taxonomy"]
            organism = Organism.get_or_save(
                taxonomy_id=taxonomy_id,
                name=record["organism"]["names"][0]["value"],
            )

            try:
                ec_number = record["protein"]["recommendedName"]["ecNumber"][0]["value"]
            except KeyError:
                ec_number = None

            # Protein name
            protein_data = record.get("protein", {})
            name = protein_data.get("recommendedName", {}).get("fullName", {}).get(
                "value"
            ) or protein_data.get("submittedName", [{}])[0].get("fullName", {}).get(
                "value"
            )

            try:
                protein = Protein.get_or_save(
                    accession_id=record["accession"],
                    sequence=record["sequence"]["sequence"],
                    mol_weight=float(record["sequence"]["mass"]),
                    ec_number=ec_number,
                    name=name,
                    seq_length=len(record["sequence"]["sequence"]),
                )
            except KeyError as e:
                logger.warning(
                    f"Error during mapping of {record['accession']} to graph model: {e}"
                )
                return

            protein.organism.connect(organism)
            self.add_reaction(record, protein)

        self.add_sites(record, protein)
        self.add_regions(record, protein)
        self.add_catalytic_activity(record, protein)
        self.add_go(record, protein)

    def add_sites(self, record: dict[str, Any], protein: Protein) -> None:
        data_dict: dict[str, list[int]] = defaultdict(list)

        for feature in record.get("features", []):
            if feature["type"] == "BINDING":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    data_dict[feature["ligand"]["name"] + "$binding"].append(position)
            elif feature["type"] == "ACT_SITE":
                for position in range(int(feature["begin"]), int(feature["end"]) + 1):
                    data_dict[feature["category"] + "$site"].append(position)

        for entry, positions in data_dict.items():
            if entry.split("$")[1] == "binding":
                annotation = Annotation.BINDING_SITE.value
            elif entry.split("$")[1] == "site":
                annotation = Annotation.ACTIVE_SITE.value

            site = Site(
                name=entry.split("$")[0],
                annotation=annotation,
            )
            site.save()

            protein.site.connect(site, {"positions": positions})

    def add_regions(self, record: dict[str, Any], protein: Protein) -> None:
        data_list: list[tuple[str, tuple[int, int]]] = []

        for feature in record.get("features", []):
            if feature["type"] == "HELIX":
                data_list.append(
                    (
                        feature["category"] + "$helix",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "STRAND":
                data_list.append(
                    (
                        feature["category"] + "$strand",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "TURN":
                data_list.append(
                    (
                        feature["category"] + "$turn",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "SIGNAL":
                data_list.append(
                    (
                        feature["category"] + "$signal",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )
            elif feature["type"] == "PROPEP":
                data_list.append(
                    (
                        feature["category"] + "$propep",
                        (int(feature["begin"]), int(feature["end"])),
                    )
                )

        for name, positions in data_list:
            if name.split("$")[1] == "helix":
                annotation = Annotation.ALPHAHELIX.value
            elif name.split("$")[1] == "strand":
                annotation = Annotation.BETASTRAND.value
            elif name.split("$")[1] == "turn":
                annotation = Annotation.TURN.value
            elif name.split("$")[1] == "signal":
                annotation = Annotation.SIGNAL.value
            elif name.split("$")[1] == "propep":
                annotation = Annotation.PROPEP.value

            region = Region(
                name=name,
                annotation=annotation,
            )
            region.save()

            protein.region.connect(region, {"start": positions[0], "end": positions[1]})

    def add_catalytic_activity(self, record: dict[str, Any], protein: Protein) -> None:
        try:
            for reference in record["comments"]:
                if reference["type"] == "CATALYTIC_ACTIVITY":
                    catalytic_annotation = Reaction.get_or_save(
                        rhea_id=str(reference["id"]) if reference.get("id") else None,
                        # Optionally, you can add name=reference["reaction"]["name"] if Reaction supports it
                    )
                    # If protein has a reaction relationship, connect it
                    if hasattr(protein, "reaction"):
                        protein.reaction.connect(catalytic_annotation)

        except Exception as e:
            logger.error(
                f"Error saving catalytic activity for {protein.accession_id}: {e}"
            )

    def get_substrates_and_products_from_rhea(
        self, rhea_id: str
    ) -> dict[str, List[str]]:
        """Fetch substrates and products from Rhea by parsing the side URI (_L = substrate, _R = product).

        Args:
            rhea_id (str or int): The Rhea reaction ID (e.g., 49528)

        Returns:
            dict: {
                'substrates': [list of chebi URIs],
                'products': [list of chebi URIs]
            }
        """
        rhea_id = rhea_id.strip().replace("RHEA:", "")
        rhea_id_str = str(rhea_id).strip()
        sparql = SPARQLWrapper("https://sparql.rhea-db.org/sparql")
        sparql.setQuery(f"""
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

        SELECT DISTINCT ?participant ?compound ?chebi ?side
        WHERE {{
        rh:{rhea_id_str} rh:side ?side .
        ?side rh:contains ?participant .
        ?participant rh:compound ?compound .
        OPTIONAL {{ ?compound rh:chebi ?chebi . }}
        OPTIONAL {{ ?compound rh:underlyingChebi ?chebi . }}
        OPTIONAL {{
            ?compound rdfs:seeAlso ?chebi .
            FILTER STRSTARTS(STR(?chebi), "http://purl.obolibrary.org/obo/CHEBI_")
        }}
        }}
        """)
        sparql.setReturnFormat(JSON)
        sparql.addCustomHttpHeader("User-Agent", "MyPythonClient/1.0")

        results_raw = sparql.query().convert()
        if not isinstance(results_raw, dict):
            raise TypeError("Expected dict from SPARQL query")

        results: dict[str, Any] = results_raw

        substrates = set()
        products = set()

        for r in results["results"]["bindings"]:
            chebi_uri = r.get("chebi", {}).get("value")
            if not chebi_uri:
                logger.info(f"No ChEBI URI found for compound {r['compound']['value']}")

            side_uri = r["side"]["value"]
            if side_uri.endswith("_L"):
                substrates.add(chebi_uri)
            elif side_uri.endswith("_R"):
                products.add(chebi_uri)

        return {"substrates": sorted(substrates), "products": sorted(products)}

    def get_smiles_from_chebi_web(self, chebi_url: str) -> Optional[str]:
        """
        Extract SMILES from the official ChEBI page using HTML scraping.
        """
        chebi_id = chebi_url.split("_")[-1]
        url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{chebi_id}"

        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")

        # Look for table rows that contain the SMILES label
        for table in soup.find_all("table", class_="chebiTableContent"):
            if not isinstance(table, Tag):
                continue
            for row in table.find_all("tr"):
                if not isinstance(row, Tag):
                    continue
                headers = row.find_all("td", class_="chebiDataHeader")
                if (
                    headers
                    and isinstance(headers[0], Tag)
                    and "SMILES" in headers[0].text
                ):
                    data_cells = row.find_all("td")
                    if data_cells:
                        return f"{data_cells[-1].text.strip()}"
        return None

    def add_reaction(self, record: dict[str, Any], protein: Protein) -> None:
        for reference in record.get("comments", []):  # Safe retrieval with .get()
            if reference.get("type") == "CATALYTIC_ACTIVITY":
                rhea_id = None  # Default value

                for db_ref in reference.get("reaction", {}).get("dbReferences", []):
                    if db_ref.get("id", "").startswith("RHEA:"):
                        rhea_id = db_ref["id"]
                        break  # Stop after finding the first match

                catalytic_annotation = Reaction.get_or_save(
                    rhea_id=rhea_id,
                )
                if rhea_id is not None:
                    self.add_molecule(rhea_id, catalytic_annotation)
                    protein.reaction.connect(catalytic_annotation)

    def add_molecule(self, rhea_id: str, reaction: Reaction) -> None:
        chebi = self.get_substrates_and_products_from_rhea(rhea_id)

        substrate_ids = chebi["substrates"]
        product_ids = chebi["products"]

        for i in substrate_ids:
            smiles = self.get_smiles_from_chebi_web(i)

            chebi_id = i.split("_")[-1]
            chebi_id = f"CHEBI:{chebi_id}"
            substrate = Molecule.get_or_save(
                chebi_id=chebi_id,
                smiles=smiles,
            )
            reaction.substrate.connect(substrate)

        for i in product_ids:
            smiles = self.get_smiles_from_chebi_web(i)

            chebi_id = i.split("_")[-1]
            chebi_id = f"CHEBI:{chebi_id}"
            product = Molecule.get_or_save(
                chebi_id=chebi_id,
                smiles=smiles,
            )
            reaction.product.connect(product)

    def add_go(self, record: dict[str, Any], protein: Protein) -> None:
        for reference in record["dbReferences"]:
            if reference["type"] == "GO":
                go_annotation = GOAnnotation.get_or_save(
                    go_id=reference["id"],
                    term=reference["properties"]["term"],
                )

                protein.go_annotation.connect(go_annotation)

    def parse_response(self, response: Response) -> Any:
        """Parse UniProt JSON response."""
        if not response.content:
            return []

        try:
            return response.json()
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse UniProt response: {e}")
            return []
