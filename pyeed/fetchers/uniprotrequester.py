import json
import requests
from typing import Union, List
from requests.exceptions import HTTPError


class UniprotRequester:
    def __init__(
        self, foreign_id: Union[str, List[str]], offset: int = 0, size: int = -1
    ):
        self._offset = offset
        self._size = size
        self.foreign_id = foreign_id
        self.url = f"https://www.ebi.ac.uk/proteins/api/proteins?offset={self._offset}&size={self._size}&accession="

    def _check_is_multiple(self) -> bool:
        """Checks if the foreign_id attribute is a list."""
        if not isinstance(self.foreign_id, list):
            return False

        if len(self.foreign_id) == 1:
            self.foreign_id = self.foreign_id[0]
            return False

        return True

    def _construct_request_string(self, request_list: list) -> str:
        return ",".join(map(str, request_list))

    def make_request(self) -> List[dict]:
        """Makes a request to the Uniprot database and returns the json response."""
        if self._check_is_multiple():
            sequence_results = []
            print(f"⬇️ Fetching {len(self.foreign_id)} entries from Uniprot...")
            for chunk in self.make_chunks(self.foreign_id, 1000):
                sequence_results.extend(
                    self.fetch(self._construct_request_string(chunk))
                )
            return sequence_results
        else:
            return self.fetch(self.foreign_id)

    def fetch(self, ids: str) -> str:
        """Makes a request to Uniprot and gets entries for the given ID(s)."""

        request_url = self.url + ids

        try:
            request = requests.get(request_url, headers={"Accept": "application/json"})
            request.raise_for_status()
        except HTTPError as http_err:
            raise HTTPError(f"ID(s) {ids} not found in Uniprot: {http_err}")

        if not request.ok:
            request.raise_for_status()

        return json.loads(request.text)

    @staticmethod
    def make_chunks(input_list: list, chunk_size: int = 100) -> List[list]:
        """
        Splits a list into chunks of a given size.
        """
        if input_list is None:
            raise ValueError("input_list cannot be None.")

        if not isinstance(input_list, list):
            raise TypeError("input_list must be a list")

        return [
            input_list[i : i + chunk_size]
            for i in range(0, len(input_list), chunk_size)
        ]
