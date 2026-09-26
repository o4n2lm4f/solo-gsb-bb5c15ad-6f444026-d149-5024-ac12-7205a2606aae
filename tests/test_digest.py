"""End-to-end tests via FastAPI TestClient."""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def digest(payload):
    return client.post("/digest", json=payload)


# ---------- happy paths ----------

def test_linear_forward_cut():
    # EcoRI-like: GAATTC, offset 1 -> cut between position 4 and 5
    r = digest({
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}],
    })
    assert r.status_code == 200
    body = r.json()
    assert body["cut_positions"] == [5]
    assert [f["sequence"] for f in body["fragments"]] == ["AAAAG", "AATTCTTTT"]
    assert [(f["start"], f["end"], f["length"]) for f in body["fragments"]] == [
        (0, 5, 5),
        (5, 14, 9),
    ]


def test_linear_reverse_complement_cut():
    # Non-palindromic recognition AACG (revcomp CGTT), offset 1.
    # Reverse hit "CGTT" starts at index 3 -> cut at 3 + (4 - 1) = 6.
    r = digest({
        "sequence": "TTTCGTTAAA",
        "topology": "linear",
        "enzymes": [{"name": "X", "recognition": "AACG", "cut_offset": 1}],
    })
    assert r.status_code == 200
    body = r.json()
    assert body["cut_positions"] == [6]
    assert [f["sequence"] for f in body["fragments"]] == ["TTTCGT", "TAAA"]


def test_circular_wraparound_site_and_modulo():
    # "GTA" occurs at index 2 and again across the origin (index 6..8).
    r = digest({
        "sequence": "ACGTACGT",
        "topology": "circular",
        "enzymes": [{"name": "Y", "recognition": "GTA", "cut_offset": 0}],
    })
    assert r.status_code == 200
    body = r.json()
    assert body["cut_positions"] == [2, 6]
    frags = body["fragments"]
    assert [f["sequence"] for f in frags] == ["GTAC", "GTAC"]
    assert (frags[0]["start"], frags[0]["end"]) == (2, 6)
    # second fragment wraps the origin: end <= start
    assert (frags[1]["start"], frags[1]["end"]) == (6, 2)


def test_circular_palindrome_site_not_doubled():
    # ACGT is its own reverse complement; each site must cut only once.
    r = digest({
        "sequence": "ACGTACGT",
        "topology": "circular",
        "enzymes": [{"name": "Z", "recognition": "ACGT", "cut_offset": 2}],
    })
    assert r.status_code == 200
    assert r.json()["cut_positions"] == [2, 6]


def test_circular_single_cut_linearizes():
    r = digest({
        "sequence": "TTGAATTCAA",
        "topology": "circular",
        "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}],
    })
    assert r.status_code == 200
    body = r.json()
    assert body["cut_positions"] == [3]
    assert body["fragment_count"] == 1
    frag = body["fragments"][0]
    assert frag["sequence"] == "AATTCAATTG"  # rotated at the cut, same length
    assert frag["length"] == 10
    assert (frag["start"], frag["end"]) == (3, 3)


def test_no_cuts_returns_full_sequence():
    for topo in ("linear", "circular"):
        r = digest({
            "sequence": "ACGTACGT",
            "topology": topo,
            "enzymes": [{"name": "None", "recognition": "TTTT", "cut_offset": 2}],
        })
        assert r.status_code == 200
        body = r.json()
        assert body["cut_positions"] == []
        assert body["fragment_count"] == 1
        assert body["fragments"][0]["sequence"] == "ACGTACGT"
        assert body["fragments"][0]["length"] == 8


def test_multiple_enzymes_combined_cuts():
    r = digest({
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [
            {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1},
            {"name": "T-cutter", "recognition": "TTTT", "cut_offset": 2},
        ],
    })
    assert r.status_code == 200
    body = r.json()
    # "TTTT" also matches as its own reverse complement "AAAA" at index 0,
    # giving a reverse-hit cut at 0 + (4 - 2) = 2.
    assert body["cut_positions"] == [2, 5, 12]
    assert [f["sequence"] for f in body["fragments"]] == [
        "AA", "AAG", "AATTCTT", "TT",
    ]


def test_lowercase_sequence_is_normalized():
    r = digest({
        "sequence": "aaaagaattctttt",
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "recognition": "gaattc", "cut_offset": 1}],
    })
    assert r.status_code == 200
    assert r.json()["fragments"][0]["sequence"] == "AAAAG"


# ---------- validation errors (422 with reasons) ----------

def test_invalid_base_rejected():
    r = digest({
        "sequence": "ACGTXACGT",
        "topology": "linear",
        "enzymes": [],
    })
    assert r.status_code == 422
    assert "invalid base" in str(r.json()).lower()


def test_empty_sequence_rejected():
    r = digest({"sequence": "", "topology": "linear", "enzymes": []})
    assert r.status_code == 422
    assert "empty" in str(r.json()).lower()


def test_offset_out_of_range_rejected():
    r = digest({
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 7}],
    })
    assert r.status_code == 422
    detail = str(r.json())
    assert "out of range" in detail and "[0, 6]" in detail


def test_invalid_recognition_base_rejected():
    r = digest({
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "Bad", "recognition": "GAATN", "cut_offset": 1}],
    })
    assert r.status_code == 422
    assert "invalid base" in str(r.json()).lower()


def test_invalid_topology_rejected():
    r = digest({
        "sequence": "ACGT",
        "topology": "supercoiled",
        "enzymes": [],
    })
    assert r.status_code == 422


def test_health():
    assert client.get("/health").json() == {"status": "ok"}
