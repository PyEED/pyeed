from typing import Generator, List
from Bio import Entrez
from tqdm import tqdm


class NCBIRequester:

    def __init__(
        self,
        foreign_id: int | List[int],
        email: str,
        api_key: str,
        retmode: str,
        db: str,
        chunk_size: int = 100,
    ) -> None:
        self.foreign_id = foreign_id
        self.email = email
        self.api_key = api_key
        self.retmode = retmode
        self.db = db
        self.chunk_size = chunk_size
        self._is_multiple: bool = isinstance(foreign_id, list)

    def _construct_request_string(self) -> str:
        """
        Constructs a request string from the foreign_id attribute.
        """

        if self._is_multiple:
            return ",".join(map(str, self.foreign_id))
        else:
            return str(self.foreign_id)

    def make_request(self) -> Generator:
        """
        Makes a request to the NCBI taxonomy database and returns the response.
        """

        Entrez.email = self.email
        Entrez.api_key = self.api_key

        if self._is_multiple:
            request_chunks = tqdm(
                self.make_chunks(self.foreign_id, self.chunk_size),
                desc=f"f⬇️ Fetching {len(self.foreign_id)} {self.db} entries for NCBI...",
            )

            results = []
            for chunk in request_chunks:
                results.extend(self.fetch(chunk))

            return results

        else:
            return self.fetch(self._construct_request_string())

    def fetch(self, request_string: str):
        """
        Fetches data from NCBI using the Entrez.efetch method.
        """
        with Entrez.efetch(
            db=self.db,
            id=request_string,
            retmode=self.retmode,
        ) as handle:
            return Entrez.read(handle)

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
