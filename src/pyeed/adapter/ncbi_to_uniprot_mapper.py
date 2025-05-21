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
    def __init__(self, ids: List[str], filename: str):
        self.ids = ids
        self.filename = filename
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
        """Fetches the UniParc and UniProt IDs for the given RefSeq IDs and saves them in JSON files regularly."""

        checksum_list = self.checksum_list(self.ids)

        id_mapping_uniprot = {}
        id_mapping_uniparc = {}

        for idx, checksum in enumerate(checksum_list):
            if isinstance(self.ids, list):
                refseq_id = self.ids[idx]
            elif isinstance(self.ids, dict):
                refseq_id = list(self.ids.values())[idx]
            else:
                raise TypeError(f"Unsupported type for self.ids: {type(self.ids)}")

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
                    
                    with open(f"{self.filename}_uniprot.json", "w") as f:
                        json.dump(id_mapping_uniprot, f, indent=2)
                    with open(f"{self.filename}_uniparc.json", "w") as f:
                        json.dump(id_mapping_uniparc, f, indent=2)

                    print(f"💾 Saved mapping for {refseq_id} at index {idx + 1}/{len(checksum_list)}")

            except httpx.RequestError as e:
                print(f"🚨 Request error for {refseq_id} (Checksum: {checksum}): {e}")
                continue
            except Exception as e:
                print(f"⚠️ Unexpected error for {refseq_id} (Checksum: {checksum}): {e}")
                continue

        print("✅ Mapping complete.")