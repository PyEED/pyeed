import asyncio

import nest_asyncio
from loguru import logger

from pyeed.adapter.ncbi_dna_mapper import NCBIDNAToPyeed
from pyeed.adapter.ncbi_protein_mapper import NCBIProteinToPyeed
from pyeed.adapter.primary_db_adapter import PrimaryDBAdapter
from pyeed.adapter.uniprot_mapper import UniprotToPyeed
from pyeed.dbconnect import DatabaseConnector
from pyeed.dbsort import DBPattern
from pyeed.embedding import (
    free_memory,
    get_batch_embeddings,
    load_model_and_tokenizer,
    update_protein_embeddings_in_db,
)


class Pyeed:
    def __init__(
        self,
        uri: str,
        user: str | None = None,
        password: str | None = None,
    ):
        self.db = DatabaseConnector(uri, user, password)

    def fetch_from_primary_db(self, ids: list[str], db="UNIPROT"):
        """
        Fetches sequences and corresponding annotations from primary sequence databases
        and adds them to local database.
        """
        nest_asyncio.apply()

        if isinstance(ids, str):
            ids = [ids]

        # Remove accessions that are already in the database
        query = """
        MATCH (p:Protein)
        RETURN collect(p.accession_id) as accessions
        """
        accessions = self.db.execute_read(query)[0]["accessions"]
        ids = [id for id in ids if id not in accessions]

        if db == DBPattern.UNIPROT.name:
            params_template = {
                "format": "json",
            }

            # set up UniProt adapter
            adapter = PrimaryDBAdapter(
                ids=ids,
                ids_attr_name="accession",
                url="https://www.ebi.ac.uk/proteins/api/proteins",
                rate_limit=10,
                n_concurrent=5,
                batch_size=5,
                data_mapper=UniprotToPyeed(),
                progress=None,
                task_id=None,
                request_params=params_template,
            )

        elif db == DBPattern.NCBI.name:
            params_template = {
                "retmode": "text",
                "rettype": "genbank",
                "db": "protein",
            }

            # set up NCBI adapter
            adapter = PrimaryDBAdapter(
                ids=ids,
                ids_attr_name="id",
                url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                rate_limit=2,
                n_concurrent=5,
                batch_size=10,
                data_mapper=NCBIProteinToPyeed(),
                progress=None,
                task_id=None,
                request_params=params_template,
            )

        asyncio.run(adapter.make_request())

    def calculate_sequence_embeddings(self, batch_size=16):
        """
        Calculates embeddings for all sequences in the database that do not have embeddings, processing in batches.
        """

        # Load the model, tokenizer, and device
        model, tokenizer, device = load_model_and_tokenizer()

        # Cypher query to retrieve proteins without embeddings and with valid sequences
        query = """
        MATCH (p:Protein)
        WHERE p.embedding IS NULL AND p.sequence IS NOT NULL
        RETURN p.accession_id AS accession, p.sequence AS sequence
        """

        # Execute the query and retrieve the results
        results = self.db.execute_read(query)
        data = [(result["accession"], result["sequence"]) for result in results]
        if not data:
            logger.info("No sequences to process.")
            return

        accessions, sequences = zip(*data)
        total_sequences = len(sequences)
        logger.debug(f"Calculating embeddings for {total_sequences} sequences.")

        # Process and save embeddings batch by batch
        for batch_start in range(0, total_sequences, batch_size):
            batch_end = min(batch_start + batch_size, total_sequences)
            batch_sequences = sequences[batch_start:batch_end]
            batch_accessions = accessions[batch_start:batch_end]
            logger.debug(
                f"Processing batch {batch_start // batch_size + 1}/"
                f"{(total_sequences + batch_size - 1) // batch_size + 1}"
            )

            # Get embeddings for the current batch
            embeddings_batch = get_batch_embeddings(
                batch_sequences, model, tokenizer, device
            )

            # Update the database for the current batch
            update_protein_embeddings_in_db(self.db, batch_accessions, embeddings_batch)

        # Free memory after processing all batches
        del model, tokenizer
        free_memory()

    def getProtein(self, accession: str):
        """
        Fetches a protein from the database by accession ID.
        """
        query = f"""
        MATCH (p:Protein {{accession_id: '{accession}'}})
        RETURN p
        """
        return self.db.execute_read(query)

    def getProteins(self):
        """
        Fetches all proteins from the database.
        """
        query = """
        MATCH (n:Protein) WHERE (n.accession_id) IS NOT NULL 
        RETURN DISTINCT  n.accession_id AS accession_id
        """
        return self.db.execute_read(query)

    def fetch_nucleotide_from_db(self, ids: list[str]):
        """
        Fetches a nucleotide sequence from the remote db and adds it to the local db.
        """

        if isinstance(ids, str):
            ids = [ids]

        # Remove accessions that are already in the database
        query = """
        MATCH (p:DNA)
        RETURN collect(p.accession_id) as accessions
        """
        accessions = self.db.execute_read(query)[0]["accessions"]
        ids = [id for id in ids if id not in accessions]

        params_template = {
            "retmode": "text",
            "rettype": "genbank",
            "db": "nuccore",
        }

        # set up NCBI adapter
        adapter = PrimaryDBAdapter(
            ids=ids,
            ids_attr_name="id",
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            rate_limit=2,
            n_concurrent=5,
            batch_size=10,
            data_mapper=NCBIDNAToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        asyncio.run(adapter.make_request())

    def getDNAs(self):
        """
        Fetches all DNA sequences from the database.
        """
        query = """
        MATCH (n:DNA) WHERE (n.accession_id) IS NOT NULL 
        RETURN DISTINCT  n.accession_id AS accession_id
        """
        return self.db.execute_read(query)

    def getDNA(self, accession: str):
        """
        Fetches a DNA sequence from the database by accession ID.
        """
        query = f"""
        MATCH (p:DNA {{accession_id: '{accession}'}})
        RETURN p
        """
        return self.db.execute_read(query)

    def fetchRemoteCodingSequences(self):
        """
        Fetches all of the coding sequences from the remote database and adds them to the local database.
        The coding sequences are saved in the protein records.
        """

        # Get all proteins and a list of coding sequences ids
        query = """
        MATCH (p:Protein) 
        WHERE p.nucleotide_id IS NOT NULL 
        RETURN p.nucleotide_id AS nucleotide_id
        """

        nucleotide_ids = self.db.execute_read(query)
        nucleotide_ids = [record["nucleotide_id"] for record in nucleotide_ids]

        logger.info(f"Fetching {len(nucleotide_ids)} coding sequences.")
        logger.info(f"Fetching coding sequences: {nucleotide_ids}")

        # Fetch the coding sequences
        self.fetch_nucleotide_from_db(nucleotide_ids)

        # we need to update the protein records with the coding sequences
        # the connection between the protein and the coding sequence is the nucleotide_id
        # the relationship is called "HAS_CODING_SEQUENCE"
        query = """
        MATCH (p:Protein) 
        WHERE p.nucleotide_id IS NOT NULL
        MATCH (n:DNA {accession_id: p.nucleotide_id})
        MERGE (p)-[:HAS_CODING_SEQUENCE]->(n)
        """
        self.db.execute_write(query)


if __name__ == "__main__":
    import os

    from dotenv import load_dotenv

    load_dotenv()

    eed = Pyeed(
        uri=os.getenv("NEO4J_URI"),
        user=os.getenv("NEO4J_USER"),
        password=os.getenv("NEO4J_PASSWORD"),
    )
    eed.db._wipe_database()

    dna_ids = ["NM_001300612.1"]

    eed.fetch_nucleotide_from_db(dna_ids)

    eed.db.close()
