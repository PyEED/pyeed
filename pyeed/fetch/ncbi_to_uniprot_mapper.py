import json
import logging
import os
import sys
from typing import List

import httpx
from crc64iso import crc64iso
from pysam import FastaFile

logger = logging.getLogger(__name__)


class NCBIToUniprotMapper:
    def __init__(self, ids: List[str], file: str):
        self.ids = ids
        self.file = file
        self.uniparc_url = "https://www.ebi.ac.uk/proteins/api/uniparc?offset=0&size=100&sequencechecksum="
        self.ncbi_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    def download_fasta(self, refseq_id: str) -> None:
        """
        Downloads a FASTA file for a given RefSeq ID using httpx and saves it locally.

        Args:
            refseq_id str: NCBI ID
        """

        params = {
            "db": "protein",
            "id": refseq_id,
            "rettype": "fasta",
            "retmode": "text",
        }

        try:
            response = httpx.get(self.ncbi_url, params=params, timeout=10.0)

            if response.status_code == 200:
                filename = f"{refseq_id}.fasta"
                with open(filename, "w") as f:
                    f.write(response.text)
                print(f"✅ Downloaded: {filename}")
            else:
                print(
                    f"❌ Failed to download {refseq_id} (Status: {response.status_code})"
                )

        except httpx.HTTPError as e:
            print(f"❌ HTTP error occurred while downloading {refseq_id}: {e}")

    def get_checksum(self, refseq_id: str) -> str:
        """Fetches and calculates the checksum for a given RefSeq ID.

        Args:
            refseq_id str: NCBI ID

        Returns:
            str: checksum ID
        """

        self.download_fasta(refseq_id)
        fa = FastaFile(f"{refseq_id}.fasta")
        seq = fa.fetch(fa.references[0])
        return f"{crc64iso.crc64(seq)}"

    def checksum_list(self, refseq_ids: List[str]) -> List[str]:
        """Creates a list of checksum IDs and deletes the FASTA files after processing.

        Args:
            refseq_ids str: NCBI IDs

        Returns:
            List[str]: cheksum IDs
        """

        checksums = []
        for refseq_id in refseq_ids:
            checksums.append(self.get_checksum(refseq_id))
            fasta_file_path = f"{refseq_id}.fasta"
            fai_file_path = f"{refseq_id}.fasta.fai"

            if os.path.exists(fasta_file_path):
                os.remove(fasta_file_path)  # Delete the fasta file

            if os.path.exists(fai_file_path):
                os.remove(fai_file_path)
        return checksums

    def execute_request(self) -> None:
        """Fetches the uniparc and uniprot ids for the given refseq ids and saves them in a json file."""

        checksum_list = self.checksum_list(self.ids)

        id_mapping_uniprot = {}
        id_mapping_uniparc = {}
        counter = 0

        for checksum in checksum_list:
            url = f"{self.uniparc_url}{checksum}"

            # perform request and get response as JSON
            with httpx.Client() as client:
                response = client.get(url, headers={"Accept": "application/json"})

            # check if the request was successful
            if response.status_code != 200:
                print(f"Request failed with status code {response.status_code}")
                response.raise_for_status()  # Raise exception for any non-200 response
                sys.exit()

            # Check if the response body is empty
            if not response.content.strip():  # Check if the body is empty
                print("The response body is empty.")
                sys.exit()

            # extracts the uniprot and the uniparc id from the repsonse and saves them in a dictionary
            response_body = response.json()
            for item in response_body:
                uniparc_id = item.get("accession", None)
                for ref in item.get("dbReference", []):
                    if ref.get("type") == "UniProtKB/TrEMBL":
                        uniprot_id = ref.get("id", None)
                        id_mapping_uniparc[self.ids[counter]] = uniparc_id
                        id_mapping_uniprot[self.ids[counter]] = uniprot_id
            counter += 1

        with open(f"{self.file}_uniprot.json", "w") as f:
            json.dump(id_mapping_uniprot, f)

        with open(f"{self.file}_uniparc.json", "w") as f:
            json.dump(id_mapping_uniparc, f)
