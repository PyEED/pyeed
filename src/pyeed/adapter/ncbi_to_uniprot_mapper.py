import json
import logging
import os
import sys
from typing import List, Optional 

import httpx
from crc64iso import crc64iso
from pysam import FastaFile

logger = logging.getLogger(__name__)


class NCBIToUniprotMapper:
    def __init__(self, ids: List[str], file_name: str):
        self.ids = ids
        self.file_name = file_name
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

    # def get_checksum(self, refseq_id: str) -> str:
    #     """Fetches and calculates the checksum for a given RefSeq ID.

    #     Args:
    #         refseq_id str: NCBI ID

    #     Returns:
    #         str: checksum ID
    #     """

    #     self.download_fasta(refseq_id)
    #     fa = FastaFile(f"{refseq_id}.fasta")
    #     seq = fa.fetch(fa.references[0])
    #     return f"{crc64iso.crc64(seq)}"

    # def checksum_list(self, refseq_ids: List[str]) -> List[str]:
    #     """Creates a list of checksum IDs and deletes the FASTA files after processing.

    #     Args:
    #         refseq_ids str: NCBI IDs

    #     Returns:
    #         List[str]: cheksum IDs
    #     """

    #     checksums = []
    #     for refseq_id in refseq_ids:
    #         checksums.append(self.get_checksum(refseq_id))
    #         fasta_file_path = f"{refseq_id}.fasta"
    #         fai_file_path = f"{refseq_id}.fasta.fai"

    #         if os.path.exists(fasta_file_path):
    #             os.remove(fasta_file_path)  # Delete the fasta file

    #         if os.path.exists(fai_file_path):
    #             os.remove(fai_file_path)
    #     return checksums

    # def execute_request(self) -> None:
    #     """Fetches the uniparc and uniprot ids for the given refseq ids and saves them in a json file."""

    #     checksum_list = self.checksum_list(self.ids)

    #     id_mapping_uniprot = {}
    #     id_mapping_uniparc = {}
    #     counter = 0

    #     for checksum in checksum_list:
    #         url = f"{self.uniparc_url}{checksum}"

    #         # perform request and get response as JSON
    #         with httpx.Client() as client:
    #             response = client.get(url, headers={"Accept": "application/json"})

    #         # check if the request was successful
    #         if response.status_code != 200:
    #             print(f"Request failed with status code {response.status_code}")
    #             response.raise_for_status()  # Raise exception for any non-200 response
    #             sys.exit()

    #         # Check if the response body is empty
    #         if not response.content.strip():  # Check if the body is empty
    #             print("The response body is empty.")
    #             sys.exit()

    #         # extracts the uniprot and the uniparc id from the repsonse and saves them in a dictionary
    #         response_body = response.json()
    #         for item in response_body:
    #             uniparc_id = item.get("accession", None)
    #             uniprot_ids = []
    #             for ref in item.get("dbReference", []):
    #                 if (
    #                     ref.get("type") == "UniProtKB/TrEMBL"
    #                     or ref.get("type") == "UniProtKB/Swiss-Prot"
    #                 ) and ref.get("active") == "Y":
    #                     uniprot_ids.append(ref.get("id"))
    #                 id_mapping_uniparc[self.ids[counter]] = uniparc_id
    #                 id_mapping_uniprot[self.ids[counter]] = uniprot_ids
    #         counter += 1

    #     with open(f"{self.file_name}_uniprot.json", "w") as f:
    #         json.dump(id_mapping_uniprot, f)

    #     with open(f"{self.file_name}_uniparc.json", "w") as f:
    #         json.dump(id_mapping_uniparc, f)

    def get_checksum(self, refseq_id: str) -> Optional[str]:
        """Fetches and calculates the checksum for a given RefSeq ID.

        Args:
            refseq_id (str): NCBI ID

        Returns:
            Optional[str]: checksum ID or None if failed
        """

        try:
            self.download_fasta(refseq_id)
            fasta_path = f"{refseq_id}.fasta"
            if not os.path.exists(fasta_path):
                raise FileNotFoundError(f"{fasta_path} not found.")

            fa = FastaFile(fasta_path)
            seq = fa.fetch(fa.references[0])
            fa.close()
            return f"{crc64iso.crc64(seq)}"

        except Exception as e:
            print(f"❌ Failed to process {refseq_id}: {e}")
            with open("missing_fasta_ids.txt", "a") as log_file:
                log_file.write(f"{refseq_id}\n")
            return None

    def checksum_list(self, refseq_ids: List[str]) -> List[str]:
        """Creates a list of checksum IDs and deletes the FASTA files after processing.

        Args:
            refseq_ids (List[str]): NCBI IDs

        Returns:
            List[str]: checksum IDs
        """

        checksums = []
        for refseq_id in refseq_ids:
            checksum = self.get_checksum(refseq_id)
            if checksum:
                checksums.append(checksum)
            else:
                print(f"⚠️ Skipping ID with missing or invalid FASTA: {refseq_id}")

            # Clean up files regardless of success
            for ext in ["fasta", "fasta.fai"]:
                file_path = f"{refseq_id}.{ext}"
                if os.path.exists(file_path):
                    try:
                        os.remove(file_path)
                    except Exception as e:
                        print(f"⚠️ Could not delete {file_path}: {e}")

        return checksums
    
    def execute_request(self) -> None:
        """Fetches the UniParc and UniProt IDs for the given RefSeq IDs and saves them in JSON files."""

        checksum_list = self.checksum_list(self.ids)

        id_mapping_uniprot = {}
        id_mapping_uniparc = {}

        for idx, checksum in enumerate(checksum_list):
            refseq_id = self.ids[idx] if idx < len(self.ids) else f"index_{idx}"
            url = f"{self.uniparc_url}{checksum}"

            try:
                with httpx.Client(timeout=10.0) as client:
                    response = client.get(url, headers={"Accept": "application/json"})

                if response.status_code != 200:
                    print(f"❌ Request failed for {refseq_id} (Checksum: {checksum}) - Status: {response.status_code}")
                    continue

                if not response.content.strip():
                    print(f"⚠️ Empty response for {refseq_id} (Checksum: {checksum})")
                    continue

                try:
                    response_body = response.json()
                except json.JSONDecodeError:
                    print(f"❌ Invalid JSON for {refseq_id} (Checksum: {checksum})")
                    continue

                for item in response_body:
                    uniparc_id = item.get("accession")
                    uniprot_ids = [
                        ref.get("id")
                        for ref in item.get("dbReference", [])
                        if ref.get("type") in {"UniProtKB/TrEMBL", "UniProtKB/Swiss-Prot"}
                        and ref.get("active") == "Y"
                    ]
                    id_mapping_uniparc[refseq_id] = uniparc_id
                    id_mapping_uniprot[refseq_id] = uniprot_ids

            except httpx.RequestError as e:
                print(f"🚨 Request error for {refseq_id} (Checksum: {checksum}): {e}")
                continue
            except Exception as e:
                print(f"⚠️ Unexpected error for {refseq_id} (Checksum: {checksum}): {e}")
                continue

        with open(f"{self.file_name}_uniprot.json", "w") as f:
            json.dump(id_mapping_uniprot, f, indent=2)

        with open(f"{self.file_name}_uniparc.json", "w") as f:
            json.dump(id_mapping_uniparc, f, indent=2)

        print("✅ Mapping complete.")