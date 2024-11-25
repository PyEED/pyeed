import asyncio
from typing import Literal

import nest_asyncio
from loguru import logger

from pyeed.adapter.ncbi_dna_mapper import NCBIDNAToPyeed
from pyeed.adapter.ncbi_protein_mapper import NCBIProteinToPyeed
from pyeed.adapter.primary_db_adapter import PrimaryDBAdapter
from pyeed.adapter.uniprot_mapper import UniprotToPyeed
from pyeed.analysis.standard_numbering import StandardNumberingTool
from pyeed.dbchat import DBChat
from pyeed.dbconnect import DatabaseConnector
from pyeed.embedding import (
    free_memory,
    get_batch_embeddings,
    load_model_and_tokenizer,
    update_protein_embeddings_in_db,
)
from pyeed.model import StandardNumbering


class Pyeed:
    """
    Main class to interact with the pyeed graph database.
    """

    def __init__(
        self,
        uri: str,
        user: str | None = None,
        password: str | None = None,
    ):
        self.db = DatabaseConnector(uri, user, password)

    def chat(self, question: str, openai_key: str, retry: bool = False) -> list[dict]:
        """Query the database using natural language via OpenAI's GPT-4 model.

        Args:
            question (str): Question to ask the database.
            openai_key (str): OpenAI API key.
            retry (bool, optional): Whether to retry once if the query if it fails.
                Defaults to False.

        Returns:
            list[dict]: List of responses from the database.
        """
        chat = DBChat(self.db)
        return chat.run(question=question, openai_key=openai_key, retry=retry)

    def fetch_from_primary_db(
        self,
        ids: list[str],
        db: Literal["uniprot", "ncbi_protein", "ncbi_nucleotide"],
    ):
        """
        Fetches sequences and corresponding annotations from primary sequence databases
        and adds them to local database.

        Args:
            ids (list[str]): List of sequence IDs to fetch from the primary database.
            db (str): Name of the primary database to fetch from. Options are "uniprot",
                "ncbi_protein", and "ncbi_nucleotide".
        """
        dbs = ("uniprot", "ncbi_protein", "ncbi_nucleotide")

        nest_asyncio.apply()

        if isinstance(ids, str):
            ids = [ids]

        # Remove accessions that are already in the database
        if db.lower() == "ncbi_nucleotide":
            query = """
            MATCH (p:DNA)
            RETURN collect(p.accession_id) as accessions
            """
        else:
            query = """
            MATCH (p:Protein)
            RETURN collect(p.accession_id) as accessions
            """

        accessions = self.db.execute_read(query)[0]["accessions"]
        ids = [id for id in ids if id not in accessions]

        logger.info(f"Fetching {len(ids)} sequences from {db}.")
        if db.lower() == "uniprot":
            self.fetch_uniprot(ids)

        elif db.lower() == "ncbi_protein":
            self.fetch_ncbi_protein(ids)

        elif db.lower() == "ncbi_protein":
            self.fetch_ncbi_nucleotide(ids)

        else:
            raise ValueError(
                f"Invalid database name '{db}'. Options are {', '.join(dbs)}."
            )

    def fetch_uniprot(self, ids: list[str]):
        """
        Fetches protein sequences from UniProt and adds them to the local database.
        """
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

        asyncio.run(adapter.make_request())
        nest_asyncio.apply()

    def fetch_ncbi_protein(self, ids: list[str]):
        """
        Fetches protein sequences from NCBI and adds them to the local database.

        Args:
            ids (list[str]): List of protein IDs to fetch from NCBI.
        """

        params_template = {
            "retmode": "text",
            "rettype": "genbank",
            "db": "protein",
        }

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
        nest_asyncio.apply()

    def fetch_ncbi_nucleotide(self, ids: list[str]):
        """
        Fetches nucleotide sequences from NCBI and adds them to the local database.

        Args:
            ids (list[str]): List of nucleotide IDs to fetch from NCBI.
        """

        params_template = {
            "retmode": "text",
            "rettype": "genbank",
            "db": "nuccore",
        }

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
        nest_asyncio.apply()

    def calculate_sequence_embeddings(
        self,
        batch_size: int = 16,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
    ):
        """
        Calculates embeddings for all sequences in the database that do not have embeddings, processing in batches.

        Args:
            batch_size (int): Number of sequences to process in each batch.
            model_name (str): Name of the model to use for calculating embeddings.
                Defaults to "facebook/esm2_t33_650M_UR50D".
                Available model can be found at https://huggingface.co/facebook/esm2_t6_8M_UR50D.
        """

        # Load the model, tokenizer, and device
        model, tokenizer, device = load_model_and_tokenizer(model_name)

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
            update_protein_embeddings_in_db(
                self.db, list(batch_accessions), embeddings_batch
            )

        # Free memory after processing all batches
        del model, tokenizer
        free_memory()

    def get_proteins(self, accession_ids: list[str]):
        """
        Fetches a protein from the database by accession ID.
        """

        if isinstance(accession_ids, str):
            accession_ids = [accession_ids]

        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $accession_ids
        RETURN p
        """
        return self.db.execute_read(query, {"accession_ids": accession_ids})

    def get_dnas(self, accession_ids: list[str]):
        """
        Fetches a DNA sequence from the database by accession ID.

        Args:
            accession_ids (list[str]): List of DNA sequence accession IDs to fetch.
        """

        if isinstance(accession_ids, str):
            accession_ids = [accession_ids]

        query = """
        MATCH (p:DNA)
        WHERE p.accession_id IN $accession_ids
        RETURN p
        """
        return self.db.execute_read(query, {"accession_ids": accession_ids})

    def fetch_dna_entries_for_proteins(self):
        """
        Fetches DNA sequences for proteins that have a nucleotide id, set in the database.
        The fetching is done from NCBI nucleotide database.
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
        self.fetch_ncbi_nucleotide(nucleotide_ids)

        # we need to update the protein records with the coding sequences
        # the connection between protein and DNA is ENCODES (fom DNA to protein)
        # but this connection could already exist, so we need to check if it exists, and onyl add it if it does not
        # the start and end positions of nucleotide sequence are stored in protein record
        # the protein record has the attribute nucleotide_id, which is the id of the coding sequence
        query = """
        MATCH (p:Protein)
        WHERE p.nucleotide_id IS NOT NULL
        RETURN p
        """
        proteins = self.db.execute_read(query)

        for protein in proteins:
            protein = protein["p"]
            # check wether the connection already exists
            # for that we take a look at the DNA node, with the nucleotide_id
            # and check if there is a connection to the protein node
            query = f"""
            MATCH (p:Protein {{accession_id: '{protein["accession_id"]}'}})
            MATCH (d:DNA {{accession_id: '{protein["nucleotide_id"]}'}})
            RETURN EXISTS((d)-[:ENCODES]->(p)) AS exists
            """
            result = self.db.execute_read(query)
            try:
                exists = result[0]["exists"]
            except IndexError:
                logger.debug(f"No connection between {protein['accession_id']} and {protein['nucleotide_id']} found.")
                continue

            if exists:
                logger.info(f"Connection between {protein['accession_id']} and {protein['nucleotide_id']} already exists.")
                continue


            query = f"""
            MATCH (p:Protein {{accession_id: '{protein["accession_id"]}'}})
            MATCH (d:DNA {{accession_id: '{protein["nucleotide_id"]}'}})
            MERGE (d)-[r:ENCODES]->(p)
            SET r.start = {protein["start"]}, r.end = {protein["end"]}
            """
            self.db.execute_write(query)
