"""Core restriction-digest logic (no web framework dependencies).

Conventions
-----------
* Cut offset is counted from the left end of the (forward) recognition
  sequence; both ends are valid, i.e. offset 0 cuts immediately before the
  first base of the recognition site and offset == len(recognition) cuts
  immediately after the last base.
* A cut at position ``p`` on a sequence of length ``n`` sits *between*
  bases ``p - 1`` and ``p`` (0-based), so valid linear cut positions are
  ``0 <= p <= n``.
* Both the recognition sequence and its reverse complement are searched.
  A reverse-complement hit at start ``i`` cuts at ``i + (k - offset)``
  where ``k`` is the recognition length.  Palindromic recognition
  sequences are searched only once so a single site is not cut twice.
* Circular templates also match sites spanning the origin; cut positions
  are reduced modulo ``n`` and de-duplicated.  Out-of-range cut positions
  on linear templates are rejected (dropped) instead of wrapped.
"""

from __future__ import annotations

from dataclasses import dataclass

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Fragment:
    index: int
    sequence: str
    length: int
    start: int  # inclusive, coordinate on the original (0-based) sequence
    end: int    # exclusive; for circular fragments that wrap the origin end <= start


def _find_starts(haystack: str, needle: str, limit: int) -> list[int]:
    """All (overlapping) start positions of ``needle`` with start < ``limit``."""
    starts = []
    i = haystack.find(needle)
    while i != -1 and i < limit:
        starts.append(i)
        i = haystack.find(needle, i + 1)
    return starts


def enzyme_cuts(sequence: str, recognition: str, offset: int, circular: bool) -> list[int]:
    """Sorted, de-duplicated cut positions produced by one enzyme."""
    n = len(sequence)
    k = len(recognition)
    if circular:
        # Extend the template so sites spanning the origin can match;
        # only matches starting within the first n bases are considered.
        repeats = (n + k - 1) // n + 1
        haystack = (sequence * repeats)[: n + k - 1]
        limit = n
    else:
        haystack = sequence
        limit = max(0, n - k + 1)

    # (pattern, offset relative to the pattern's left end on the forward strand)
    targets = [(recognition, offset)]
    rc = reverse_complement(recognition)
    if rc != recognition:
        targets.append((rc, k - offset))

    cuts: set[int] = set()
    for pattern, off in targets:
        for i in _find_starts(haystack, pattern, limit):
            cut = i + off
            if circular:
                cuts.add(cut % n)
            elif 0 <= cut <= n:
                cuts.add(cut)
            # linear templates: out-of-bounds cut positions are rejected
    return sorted(cuts)


def build_fragments(sequence: str, cuts: list[int], circular: bool) -> list[Fragment]:
    """Slice ``sequence`` at sorted cut positions, in original strand order."""
    n = len(sequence)
    if not cuts:
        return [Fragment(index=0, sequence=sequence, length=n, start=0, end=n)]

    if not circular:
        interior = [c for c in cuts if 0 < c < n]
        bounds = [0] + interior + [n]
        return [
            Fragment(index=i, sequence=sequence[s:e], length=e - s, start=s, end=e)
            for i, (s, e) in enumerate(zip(bounds, bounds[1:]))
        ]

    fragments = []
    for i, s in enumerate(cuts):
        e = cuts[(i + 1) % len(cuts)]
        if e > s:
            frag_seq = sequence[s:e]
        else:  # wraps the origin (also the single-cut case, where e == s)
            frag_seq = sequence[s:] + sequence[:e]
        fragments.append(
            Fragment(index=i, sequence=frag_seq, length=len(frag_seq), start=s, end=e)
        )
    return fragments
