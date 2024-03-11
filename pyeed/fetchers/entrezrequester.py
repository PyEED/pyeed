from typing import List
from Bio import Entrez, SeqIO
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
        rettype: str = None,
    ) -> None:
        self.foreign_id = foreign_id
        self.email = email
        self.api_key = api_key
        self.rettype = rettype
        self.retmode = retmode
        self.db = db
        self.chunk_size = chunk_size
        self._is_multiple: bool = isinstance(foreign_id, list)

    @staticmethod
    def _construct_request_string(request_list: list) -> str:
        """
        Constructs a request string from the foreign_id attribute.
        """
        return ",".join(map(str, request_list))

    def make_request(self) -> List:
        """
        Makes a request to the NCBI taxonomy database and returns the response.
        """

        if self._is_multiple:
            sequence_results = []
            request_chunks = tqdm(
                self.make_chunks(self.foreign_id, self.chunk_size),
                desc=f"f⬇️ Fetching {len(self.foreign_id)} {self.db} entries for NCBI...",
            )

            for chunk in request_chunks:
                sequence_results.extend(
                    self.fetch(self._construct_request_string(chunk))
                )

            return sequence_results

        else:
            return self.fetch(self.foreign_id)

    def fetch(self, request_string: str) -> List[dict]:
        """
        Fetches data from NCBI using the Entrez.efetch method.
        """

        Entrez.email = self.email
        Entrez.api_key = self.api_key
        with Entrez.efetch(
            db=self.db,
            id=request_string,
            retmode=self.retmode,
            rettype=self.rettype,
        ) as handle:
            if self.rettype == "genbank":
                results = []
                for record in SeqIO.parse(handle, "genbank"):
                    results.append(record)
                return results

            else:
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


if __name__ == "__main__":
    # Example usage
    email = "d@des.de"
    # ids = [9606, 10090, 10116]
    # ids = 9606
    # rentrez = NCBIRequester(
    #     foreign_id=ids,
    #     email=email,
    #     api_key=None,
    #     retmode="xml",
    #     db="taxonomy",
    # )

    protein_ids = ["NP_001191", "UCS38941.1", "NP_001191", "UCS38941.1"]

    rentrez = NCBIRequester(
        foreign_id=protein_ids,
        email=email,
        api_key=None,
        retmode="text",
        db="protein",
        rettype="genbank",
    )

    res = rentrez.make_request()

    print(res)
