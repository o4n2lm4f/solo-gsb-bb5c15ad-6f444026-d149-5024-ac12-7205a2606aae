"""Request/response models with validation (clear 422 reasons)."""

from __future__ import annotations

from typing import List, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

_VALID_BASES = frozenset("ACGT")


def _check_dna(value: str, what: str) -> str:
    seq = value.strip().upper()
    if not seq:
        raise ValueError(f"{what} must not be empty")
    bad = sorted(set(seq) - _VALID_BASES)
    if bad:
        raise ValueError(
            f"{what} contains invalid base(s): {', '.join(bad)} "
            "(only A, C, G, T are allowed)"
        )
    return seq


class Enzyme(BaseModel):
    name: str = Field(..., examples=["EcoRI"])
    recognition: str = Field(..., examples=["GAATTC"])
    cut_offset: int = Field(
        ...,
        description="Cut offset from the left end of the forward recognition "
        "sequence; both ends are valid: 0 <= cut_offset <= len(recognition)",
        examples=[1],
    )

    @field_validator("name")
    @classmethod
    def _name_not_empty(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("enzyme name must not be empty")
        return v

    @field_validator("recognition")
    @classmethod
    def _recognition_is_dna(cls, v: str) -> str:
        return _check_dna(v, "recognition sequence")

    @model_validator(mode="after")
    def _offset_in_range(self) -> "Enzyme":
        k = len(self.recognition)
        if not 0 <= self.cut_offset <= k:
            raise ValueError(
                f"cut_offset {self.cut_offset} is out of range for enzyme "
                f"'{self.name}': must be within [0, {k}] "
                "(recognition length, both ends inclusive)"
            )
        return self


class DigestRequest(BaseModel):
    sequence: str = Field(..., examples=["AAAAGAATTCTTTT"])
    topology: Literal["linear", "circular"] = "linear"
    enzymes: List[Enzyme] = Field(default_factory=list)

    @field_validator("sequence")
    @classmethod
    def _sequence_is_dna(cls, v: str) -> str:
        return _check_dna(v, "sequence")


class FragmentOut(BaseModel):
    index: int
    sequence: str
    length: int
    start: int
    end: int


class EnzymeResult(BaseModel):
    name: str
    recognition: str
    cut_offset: int
    cuts: List[int]


class DigestResponse(BaseModel):
    topology: str
    sequence_length: int
    cut_positions: List[int]
    fragment_count: int
    fragments: List[FragmentOut]
    enzymes: List[EnzymeResult]
