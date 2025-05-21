import asyncio
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Literal

import nest_asyncio
import torch
from loguru import logger
import os 

from pyeed.adapter.ncbi_dna_mapper import NCBIDNAToPyeed
from pyeed.adapter.ncbi_protein_mapper import NCBIProteinToPyeed
from pyeed.adapter.ncbi_to_uniprot_mapper import NCBIToUniprotMapper
from pyeed.adapter.primary_db_adapter import PrimaryDBAdapter
from pyeed.adapter.uniprot_mapper import UniprotToPyeed
from pyeed.dbchat import DBChat
from pyeed.dbconnect import DatabaseConnector
from pyeed.embedding import (
    load_model_and_tokenizer,
    process_batches_on_gpu,
)


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

    def chat(
        self, question: str, openai_key: str, retry: bool = False
    ) -> list[dict[str, Any]]:
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
    ) -> None:
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
        # count how many sequences are already in the database
        logger.info(f"Found {len(accessions)} sequences in the database.")

        logger.info(f"Fetching {len(ids)} sequences from {db}.")
        if db.lower() == "uniprot":
            self.fetch_uniprot(ids)

        elif db.lower() == "ncbi_protein":
            self.fetch_ncbi_protein(ids)

        elif db.lower() == "ncbi_nucleotide":
            self.fetch_ncbi_nucleotide(ids)

        else:
            raise ValueError(
                f"Invalid database name '{db}'. Options are {', '.join(dbs)}."
            )

    def fetch_uniprot(self, ids: list[str]) -> None:
        """
        Fetches protein sequences from UniProt and adds them to the local database.
        """
        params_template = {
            "format": "json",
        }

        # set up UniProt adapter
        adapter = PrimaryDBAdapter(
            ids=ids,
            id_param_name="accession",
            url="https://www.ebi.ac.uk/proteins/api/proteins",
            rate_limit=10,
            max_concurrent=5,
            batch_size=5,
            data_mapper=UniprotToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        # Fix: call nest_asyncio.apply() first, then run the adapter's coroutine
        nest_asyncio.apply()
        asyncio.get_event_loop().run_until_complete(adapter.execute_requests())

    def fetch_ncbi_protein(self, ids: list[str]) -> None:
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
            id_param_name="id",
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            rate_limit=2,
            max_concurrent=5,
            batch_size=10,
            data_mapper=NCBIProteinToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        # Fix: use run_until_complete instead of asyncio.run
        nest_asyncio.apply()
        asyncio.get_event_loop().run_until_complete(adapter.execute_requests())

    def fetch_ncbi_nucleotide(self, ids: list[str]) -> None:
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
            id_param_name="id",
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            rate_limit=2,
            max_concurrent=5,
            batch_size=10,
            data_mapper=NCBIDNAToPyeed(),
            progress=None,
            task_id=None,
            request_params=params_template,
        )

        # Fix: apply nest_asyncio and then run the coroutine with the event loop
        nest_asyncio.apply()
        asyncio.get_event_loop().run_until_complete(adapter.execute_requests())

    def database_id_mapper(self, ids: list[str], file: str) -> None:
        """
        Maps IDs from one database to another using the UniProt ID mapping service

        Args:
            ids (list[str]): List of IDs to map.
        """

        mapper = NCBIToUniprotMapper(ids, file)
        mapper.execute_request()

        nest_asyncio.apply()

    def calculate_sequence_embeddings(
        self,
        batch_size: int = 16,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        num_gpus: int = 1,  # Number of GPUs to use
    ) -> None:
        """
        Calculates embeddings for all sequences in the database that do not have embeddings,
        distributing the workload across available GPUs.

        Args:
            batch_size (int): Number of sequences to process in each batch.
            model_name (str): Model used for calculating embeddings.
            num_gpus (int, optional): Number of GPUs to use. If None, use all available GPUs.
        """

        # Get the available GPUs
        os.environ.pop("CUDA_VISIBLE_DEVICES", None)
        
        print("CUDA available:", torch.cuda.is_available())
        print("Device count:", torch.cuda.device_count())
        print("Torch CUDA version:", torch.version.cuda)
        print("Torch version:", torch.__version__)
        
        available_gpus = torch.cuda.device_count()
        if num_gpus is None or num_gpus > available_gpus:
            num_gpus = available_gpus

        if num_gpus == 0:
            logger.warning("No GPU available! Running on CPU.")

        # Load separate models for each GPU
        devices = (
            [torch.device(f"cuda:{i}") for i in range(num_gpus)]
            if num_gpus > 0
            else [torch.device("cpu")]
        )

        models_and_tokenizers = [
            load_model_and_tokenizer(model_name, device) for device in devices
        ]

        # Retrieve sequences without embeddings
        query = """
        MATCH (p:Protein)
        WHERE p.embedding IS NULL AND p.sequence IS NOT NULL
        RETURN p.accession_id AS accession, p.sequence AS sequence
        """
        results = self.db.execute_read(query)
        data = [(result["accession"], result["sequence"]) for result in results]

        if not data:
            logger.info("No sequences to process.")
            return

        accessions, sequences = zip(*data)
        total_sequences = len(sequences)
        logger.debug(f"Total sequences to process: {total_sequences}")

        # Split the data into num_gpus chunks
        gpu_batches = [
            list(zip(accessions[i::num_gpus], sequences[i::num_gpus]))
            for i in range(num_gpus)
        ]

        start_time = time.time()

        # Process batches in parallel across GPUs
        with ThreadPoolExecutor(max_workers=num_gpus) as executor:
            futures = []
            for i, gpu_data in enumerate(gpu_batches):
                if not gpu_data:
                    continue  # Skip empty GPU batches

                model, tokenizer, device = models_and_tokenizers[i]
                futures.append(
                    executor.submit(
                        process_batches_on_gpu,
                        gpu_data,
                        batch_size,
                        model,
                        tokenizer,
                        self.db,
                        device,
                    )
                )

            for future in futures:
                future.result()  # Wait for all threads to complete

        end_time = time.time()
        logger.info(
            f"Total embedding calculation time: {end_time - start_time:.2f} seconds"
        )

        # Cleanup
        for model, _, _ in models_and_tokenizers:
            del model

    def get_proteins(self, accession_ids: list[str]) -> list[dict[str, Any]]:
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

    def get_dnas(self, accession_ids: list[str]) -> list[dict[str, Any]]:
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

    def fetch_dna_entries_for_proteins(self, ids: list[str] | None = None) -> None:
        """
        Fetches DNA sequences for proteins that have a nucleotide id, set in the database.
        The fetching is done from NCBI nucleotide database in batches.

        Args:
            ids (list[str], optional): List of protein IDs to fetch DNA sequences for.
                Defaults to None.
        """
        BATCH_SIZE = 100

        # Get all proteins and a list of coding sequences ids
        if ids is None:
            query = """
            MATCH (p:Protein) 
            WHERE p.nucleotide_id IS NOT NULL
            RETURN p.nucleotide_id AS nucleotide_id
            """
            response = self.db.execute_read(query)
        else:
            query = """
            MATCH (p:Protein) 
            WHERE p.nucleotide_id IS NOT NULL AND p.accession_id IN $ids
            RETURN p.nucleotide_id AS nucleotide_id
            """
            response = self.db.execute_read(query, {"ids": ids})

        nucleotide_ids = [str(record["nucleotide_id"]) for record in response]

        logger.info(f"Found {len(nucleotide_ids)} coding sequences.")

        # Process nucleotide IDs in batches to check which ones are already in DB
        all_existing_sequences = set()
        for i in range(0, len(nucleotide_ids), BATCH_SIZE):
            batch_ids = nucleotide_ids[i : i + BATCH_SIZE]
            try:
                query = """
                MATCH (n:DNA)
                WHERE n.accession_id IN $nucleotide_ids
                RETURN n.accession_id AS accession_id
                """
                coding_sequences_resp = self.db.execute_read(
                    query, {"nucleotide_ids": batch_ids}
                )
                batch_existing = {
                    str(record["accession_id"]) for record in coding_sequences_resp
                }
                all_existing_sequences.update(batch_existing)
            except Exception as e:
                logger.error(
                    f"Error checking existing sequences for batch {i}: {str(e)}"
                )
                continue

        # Filter out existing sequences
        nucleotide_ids = [
            id for id in nucleotide_ids if id not in all_existing_sequences
        ]

        logger.info(f"Fetching {len(nucleotide_ids)} new coding sequences.")

        # Fetch coding sequences in batches
        for i in range(0, len(nucleotide_ids), BATCH_SIZE):
            try:
                batch_ids = nucleotide_ids[i : i + BATCH_SIZE]
                self.fetch_ncbi_nucleotide(batch_ids)
                logger.info(f"Successfully fetched batch {i // BATCH_SIZE + 1}")
            except Exception as e:
                logger.error(f"Error fetching batch {i // BATCH_SIZE + 1}: {str(e)}")
                continue

        # Process protein-DNA relationships in batches
        if ids is None:
            query = """
            MATCH (p:Protein)
            WHERE p.nucleotide_id IS NOT NULL
            RETURN p
            """
            proteins = self.db.execute_read(query)
        else:
            query = """
            MATCH (p:Protein)
            WHERE p.nucleotide_id IS NOT NULL AND p.accession_id IN $ids
            RETURN p
            """
            proteins = self.db.execute_read(query, {"ids": ids})

        for i in range(0, len(proteins), BATCH_SIZE):
            try:
                batch_proteins = proteins[i : i + BATCH_SIZE]

                # Build batch query for checking existing relationships
                batch_check_query = """
                UNWIND $proteins AS protein
                MATCH (p:Protein {accession_id: protein.p.accession_id})
                MATCH (d:DNA {accession_id: protein.p.nucleotide_id})
                RETURN 
                    protein.p.accession_id AS protein_id,
                    protein.p.nucleotide_id AS dna_id,
                    EXISTS((d)-[:ENCODES]->(p)) AS exists,
                    protein.p.nucleotide_start AS start,
                    protein.p.nucleotide_end AS end
                """

                results = self.db.execute_read(
                    batch_check_query, {"proteins": batch_proteins}
                )

                # Filter relationships that need to be created
                new_relationships = []
                for result in results:
                    if not result["exists"]:
                        new_relationships.append(
                            {
                                "protein_id": result["protein_id"],
                                "dna_id": result["dna_id"],
                                "start": result["start"],
                                "end": result["end"],
                            }
                        )
                    else:
                        logger.info(
                            f"Connection between {result['protein_id']} and {result['dna_id']} already exists."
                        )

                if new_relationships:
                    # Create new relationships in batch
                    batch_create_query = """
                    UNWIND $relationships AS rel
                    MATCH (p:Protein {accession_id: rel.protein_id})
                    MATCH (d:DNA {accession_id: rel.dna_id})
                    MERGE (d)-[r:ENCODES]->(p)
                    SET r.start = rel.start, r.end = rel.end
                    """
                    self.db.execute_write(
                        batch_create_query, {"relationships": new_relationships}
                    )
                    logger.info(
                        f"Successfully processed relationship batch {i // BATCH_SIZE + 1}"
                    )
            except Exception as e:
                logger.error(
                    f"Error processing relationship batch {i // BATCH_SIZE + 1}: {str(e)}"
                )
                continue

    def create_coding_sequences_regions(self) -> None:
        """
        Creates coding sequences regions for all proteins in the database.

        It finds the nucleotide start and end positions and create a Region object for the corresponding DNA sequence.
        Create the region object with the right annotation. And then connect it to the DNA sequence.
        """

        # in case of multiple DNA entires per Protein we need to create a Region for each DNA entry
        # some of the DNA entries might even not have start and end vlaues on the ENCODES endge, in this case please take the entire sequence length
        """
        This Cypher query creates coding sequence regions for DNA sequences that don't already have them.
        
        The query:
        1. Finds all Protein-DNA pairs connected by an ENCODES relationship
        2. Filters for cases where the DNA doesn't already have a coding sequence Region for that protein
        3. Creates a new Region node with 'coding sequence' annotation and the protein's ID
        4. Creates a HAS_REGION relationship from the DNA to the new Region
        5. Sets the start position to either:
           - The start value from the ENCODES relationship if it exists
           - 0 (beginning of sequence) if no start value is specified
        6. Sets the end position to either:
           - The end value from the ENCODES relationship if it exists
           - The full DNA sequence length minus 1 if no end value is specified
        """
        query = """
        MATCH (d:DNA)-[rel_encode:ENCODES]->(p:Protein)
        WHERE NOT EXISTS((d)-[:HAS_REGION]->(:Region {annotation: 'coding sequence', sequence_id: p.accession_id}))
        CREATE (r:Region {annotation: 'coding sequence', sequence_id: p.accession_id})
        CREATE (d)-[rel:HAS_REGION {
            start: CASE 
                WHEN rel_encode.start IS NOT NULL THEN rel_encode.start - 1
                ELSE 0 
            END, 
            end: CASE 
                WHEN rel_encode.end IS NOT NULL THEN rel_encode.end - 1
                ELSE size(d.sequence) - 1 
            END
        }]->(r)
        """
        self.db.execute_write(query)

        # for dna where ther is no protein encoded and no Region with coding sequence annotation, create a Region with the entire sequence length
        # make the start at 0 and the end at the sequence length minus 1
        query = """
        MATCH (d:DNA)
        WHERE NOT EXISTS((d)-[:HAS_REGION]->(:Region {annotation: 'coding sequence', sequence_id: d.accession_id}))
        CREATE (r:Region {annotation: 'coding sequence', sequence_id: d.accession_id})
        CREATE (d)-[:HAS_REGION {start: 0, end: size(d.sequence) - 1}]->(r)
        """
        self.db.execute_write(query)

        # Log the number of regions created
        count_query = """
        MATCH (d:DNA)-[:HAS_REGION]->(r:Region {annotation: 'coding sequence'})
        RETURN count(r) as region_count
        """
        result = self.db.execute_read(count_query)
        logger.info(f"Created {result[0]['region_count']} coding sequence regions")
