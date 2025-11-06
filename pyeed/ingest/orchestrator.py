"""Fluent API for orchestrating data ingestion and enrichment."""

from __future__ import annotations

from typing import Self

from ..db.neo4j import Database
from ..embedding.sinks import EmbeddingSink
from .embed import embed_proteins
from .functions import (
    enrich_molecules,
    enrich_reactions,
    fetch_proteins_by_ids,
    fetch_proteins_by_interpro,
)
from .progress import create_progress

__all__ = ["Ingester"]


class Ingester:
    """Fluent API for ingesting and enriching protein data.

    Provides a chainable interface for configuring and executing
    data ingestion pipelines.

    Example:
        >>> db = Database()
        >>> await Ingester(db).from_interpro("IPR002133").run()

        >>> # With enrichment
        >>> await (Ingester(db)
        ...        .from_ids(["P12345", "Q9Y6K9"])
        ...        .with_reactions()
        ...        .with_molecules()
        ...        .run())

        >>> # Full pipeline
        >>> await (Ingester(db)
        ...        .from_interpro("IPR002133")
        ...        .with_full_enrichment()
        ...        .quiet()
        ...        .run())
    """

    def __init__(
        self,
        db: Database,
    ) -> None:
        """Initialize the ingester.

        Args:
            db: Database instance to use for ingestion.
        """
        self.db = db
        self._source: tuple[str, list[str] | str] | None = None
        self._enrich_reactions = False
        self._enrich_molecules = False
        self._embed_proteins = False
        self._show_progress = True
        self._skip_schema_sync = False
        self._chunk_size = 30
        self._page_size = 30
        self._batch_size = 200
        # Embedding-specific configuration
        self._embedding_model = "facebook/esm2_t33_650M_UR50D"
        self._pooling_methods: list | None = None  # None means [mean_pooling] default
        self._embedding_sinks: list[EmbeddingSink] = []
        self._embedding_batch_size = 16
        self._embedding_chunk_size = 1000
        self._embedding_model_dtype = "float16"
        self._embedding_return_dtype = "float16"

    def from_ids(self, accessions: list[str]) -> Self:
        """Fetch proteins by UniProt accession IDs.

        Args:
            accessions: List of UniProt accession IDs.

        Returns:
            Self for method chaining.
        """
        self._source = ("ids", accessions)
        return self

    def from_interpro(self, interpro_id: str) -> Self:
        """Fetch proteins by InterPro ID.

        Args:
            interpro_id: InterPro ID (e.g., "IPR002133").

        Returns:
            Self for method chaining.
        """
        self._source = ("interpro", interpro_id)
        return self

    def with_reactions(self) -> Self:
        """Include Rhea reaction enrichment.

        Returns:
            Self for method chaining.
        """
        self._enrich_reactions = True
        return self

    def with_molecules(self) -> Self:
        """Include ChEBI molecule enrichment.

        Returns:
            Self for method chaining.
        """
        self._enrich_molecules = True
        return self

    def with_full_enrichment(self) -> Self:
        """Include all enrichments (reactions + molecules).

        Returns:
            Self for method chaining.
        """
        self._enrich_reactions = True
        self._enrich_molecules = True
        return self

    def quiet(self) -> Self:
        """Disable progress bars.

        Returns:
            Self for method chaining.
        """
        self._show_progress = False
        return self

    def skip_schema_sync(self) -> Self:
        """Skip automatic schema synchronization.

        Returns:
            Self for method chaining.
        """
        self._skip_schema_sync = True
        return self

    def with_batch_size(self, batch_size: int) -> Self:
        """Set batch size for database operations.

        Args:
            batch_size: Number of items per batch.

        Returns:
            Self for method chaining.
        """
        self._batch_size = batch_size
        return self

    def with_chunk_size(self, chunk_size: int) -> Self:
        """Set chunk size for API requests.

        Args:
            chunk_size: Number of items per API request.

        Returns:
            Self for method chaining.
        """
        self._chunk_size = chunk_size
        return self

    def with_page_size(self, page_size: int) -> Self:
        """Set page size for UniProt API pagination.

        Args:
            page_size: Results per page.

        Returns:
            Self for method chaining.
        """
        self._page_size = page_size
        return self

    def with_embeddings(
        self,
        sinks: list[EmbeddingSink] | None = None,
        pooling_methods: list | None = None,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        batch_size: int = 16,
        chunk_size: int = 1000,
        model_dtype: str = "float16",
        return_dtype: str = "float16",
    ) -> Self:
        """Include protein embedding computation.

        Args:
            sinks: List of EmbeddingSink instances
            pooling_methods: List of pooling functions (default: [mean_pooling])
            model_name: ESM2 model identifier
            batch_size: Sequences per GPU batch
            chunk_size: Sequences per embedding chunk
            model_dtype: Model computation dtype (float32/float16/bfloat16)
            return_dtype: Returned embedding dtype (float32/float16/bfloat16)

        Returns:
            Self for method chaining
        """
        self._embed_proteins = True
        self._embedding_sinks = sinks or []
        if pooling_methods:
            self._pooling_methods = pooling_methods
        self._embedding_model = model_name
        self._embedding_batch_size = batch_size
        self._embedding_chunk_size = chunk_size
        self._embedding_model_dtype = model_dtype
        self._embedding_return_dtype = return_dtype
        return self

    async def run(self) -> dict[str, int]:
        """Execute the configured ingestion pipeline.

        Returns:
            Dictionary with counts: {
                "proteins": int,
                "reactions": int,
                "molecules": int,
                "embeddings": int  (if embeddings are enabled)
            }

        Raises:
            ValueError: If no data source was configured.

        Example:
            >>> db = Database()
            >>> ingester = Ingester(db).from_interpro("IPR002133")
            >>> stats = await ingester.run()
            >>> stats["proteins"] > 0
            True
        """
        if self._source is None:
            raise ValueError("No data source configured. Use from_ids() or from_interpro().")

        stats: dict[str, int] = {"proteins": 0, "reactions": 0, "molecules": 0, "embeddings": 0}

        source_type, source_value = self._source

        # Fetch proteins
        if source_type == "ids":
            proteins = await fetch_proteins_by_ids(
                self.db,
                source_value,  # type: ignore
                chunk_size=self._chunk_size,
                page_size=self._page_size,
                batch_size=self._batch_size,
                show_progress=self._show_progress,
                skip_schema_sync=self._skip_schema_sync,
            )
        else:  # interpro
            proteins = await fetch_proteins_by_interpro(
                self.db,
                source_value,  # type: ignore
                chunk_size=self._chunk_size,
                page_size=self._page_size,
                batch_size=self._batch_size,
                show_progress=self._show_progress,
                skip_schema_sync=self._skip_schema_sync,
            )

        stats["proteins"] = len(proteins)

        # Create shared progress instance for enrichment/embedding steps if needed
        needs_progress = self._enrich_reactions or self._enrich_molecules or self._embed_proteins
        if self._show_progress and needs_progress:
            progress = create_progress()
            with progress:
                # Enrich reactions
                if self._enrich_reactions:
                    stats["reactions"] = await enrich_reactions(
                        self.db,
                        batch_size=self._batch_size,
                        show_progress=self._show_progress,
                        progress=progress,
                    )

                # Enrich molecules
                if self._enrich_molecules:
                    stats["molecules"] = await enrich_molecules(
                        self.db,
                        batch_size=self._batch_size,
                        show_progress=self._show_progress,
                        progress=progress,
                    )

                # Compute embeddings
                if self._embed_proteins:
                    sinks = self._embedding_sinks
                    stats["embeddings"] = await embed_proteins(
                        self.db,
                        model_name=self._embedding_model,
                        pooling_methods=self._pooling_methods,
                        sinks=sinks,
                        batch_size=self._embedding_batch_size,
                        chunk_size=self._embedding_chunk_size,
                        model_dtype=self._embedding_model_dtype,
                        return_dtype=self._embedding_return_dtype,
                        show_progress=self._show_progress,
                        progress=progress,
                    )
        else:
            # No shared progress needed
            if self._enrich_reactions:
                stats["reactions"] = await enrich_reactions(
                    self.db,
                    batch_size=self._batch_size,
                    show_progress=self._show_progress,
                )

            if self._enrich_molecules:
                stats["molecules"] = await enrich_molecules(
                    self.db,
                    batch_size=self._batch_size,
                    show_progress=self._show_progress,
                )

            if self._embed_proteins:
                sinks = self._embedding_sinks
                stats["embeddings"] = await embed_proteins(
                    self.db,
                    model_name=self._embedding_model,
                    pooling_methods=self._pooling_methods,
                    sinks=sinks,
                    batch_size=self._embedding_batch_size,
                    chunk_size=self._embedding_chunk_size,
                    model_dtype=self._embedding_model_dtype,
                    return_dtype=self._embedding_return_dtype,
                    show_progress=self._show_progress,
                )

        return stats
