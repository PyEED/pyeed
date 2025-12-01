from typing import Any

import fastapi

app = fastapi.FastAPI()


@app.get("/")
def read_root() -> dict[str, str]:
    return {"message": "Hello, World!"}


@app.get("/id-similarity-search")
def id_similarity_search(query: str, n_hits: int = 10) -> list[dict[str, Any]]:
    pass


@app.get("/sequence-similarity-search")
def sequence_similarity_search(query: str, n_hits: int = 10) -> list[dict[str, Any]]:
    pass
