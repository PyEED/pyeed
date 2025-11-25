"""Methods to add individual nodes to the database."""

from __future__ import annotations

import asyncio
import concurrent.futures

from pyeed.db.neo4j import GraphDB
from pyeed.db.queries import (
    _upsert_nodes_with_session,
    create_reaction_molecule_relationships,
)
from pyeed.ingest.model.molecule import Molecule
from pyeed.ingest.model.reaction import Reaction


def _run_async(coro):
    """Run async coroutine synchronously, handling both regular Python and Jupyter notebooks.

    Works with existing event loops (notebooks) or creates a new one (regular Python).
    """
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        # No running loop, create a new one
        return asyncio.run(coro)
    else:
        # Running loop exists (e.g., in Jupyter), create a new thread with its own event loop
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(asyncio.run, coro)
            return future.result()


def add_reactions(
    db: GraphDB,
    reaction_data: list[tuple[Reaction, list[Molecule], list[Molecule]]],
    tx_size: int = 5000,
) -> None:
    """Add reactions with substrates and products to the database.

    Synchronous wrapper that works in both regular Python and Jupyter notebooks.

    Args:
        db: GraphDB instance
        reaction_data: List of tuples (reaction, substrates, products)
            where substrates and products are lists of Molecule objects
        tx_size: Transaction batch size

    Example:
        >>> reaction = Reaction(id="RHEA:10000", description="Some reaction")
        >>> substrates = [Molecule(id="CHEBI:15377"), Molecule(id="CHEBI:15378")]
        >>> products = [Molecule(id="CHEBI:15379")]
        >>> add_reactions(db, [(reaction, substrates, products)])
    """
    if not reaction_data:
        return

    # Collect all nodes to upsert
    all_reactions: list[Reaction] = []
    all_molecules: list[Molecule] = []
    substrate_map: dict[str, list[str]] = {}
    product_map: dict[str, list[str]] = {}

    for reaction, substrates, products in reaction_data:
        all_reactions.append(reaction)

        # Collect substrates
        if substrates:
            all_molecules.extend(substrates)
            substrate_map[reaction.id] = [mol.id for mol in substrates]

        # Collect products
        if products:
            all_molecules.extend(products)
            product_map[reaction.id] = [mol.id for mol in products]

    async def _async_add() -> None:
        # Upsert all nodes and create relationships in a single session
        # This ensures nodes are visible when creating relationships
        async with db.async_driver.session() as session:
            # Step 1: Upsert all nodes (reactions and molecules)
            all_nodes = all_reactions + all_molecules
            if all_nodes:
                await _upsert_nodes_with_session(session, all_nodes, tx_size=tx_size)

            # Step 2: Create relationships
            if substrate_map or product_map:
                await create_reaction_molecule_relationships(
                    session=session,
                    substrate_map=substrate_map,
                    product_map=product_map,
                    tx_size=tx_size,
                )

    _run_async(_async_add())
