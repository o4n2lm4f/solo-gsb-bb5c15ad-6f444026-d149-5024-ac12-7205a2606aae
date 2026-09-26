"""DNA restriction-enzyme digest simulator API."""

from __future__ import annotations

from fastapi import FastAPI

from .digest import build_fragments, enzyme_cuts
from .schemas import (
    DigestRequest,
    DigestResponse,
    EnzymeResult,
    FragmentOut,
)

app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description=(
        "Simulates restriction-enzyme digestion of linear or circular DNA. "
        "Recognition sequences and their reverse complements are searched; "
        "circular templates match sites spanning the origin."
    ),
)


@app.get("/health", tags=["meta"])
def health() -> dict:
    return {"status": "ok"}


@app.post("/digest", response_model=DigestResponse, tags=["digest"])
def digest(req: DigestRequest) -> DigestResponse:
    circular = req.topology == "circular"

    enzyme_results: list[EnzymeResult] = []
    all_cuts: set[int] = set()
    for enz in req.enzymes:
        cuts = enzyme_cuts(req.sequence, enz.recognition, enz.cut_offset, circular)
        enzyme_results.append(
            EnzymeResult(
                name=enz.name,
                recognition=enz.recognition,
                cut_offset=enz.cut_offset,
                cuts=cuts,
            )
        )
        all_cuts.update(cuts)

    cuts_sorted = sorted(all_cuts)
    fragments = build_fragments(req.sequence, cuts_sorted, circular)

    return DigestResponse(
        topology=req.topology,
        sequence_length=len(req.sequence),
        cut_positions=cuts_sorted,
        fragment_count=len(fragments),
        fragments=[FragmentOut(**vars(f)) for f in fragments],
        enzymes=enzyme_results,
    )
