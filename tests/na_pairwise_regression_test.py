"""
Tests written by Claude that may catch some problems
"""

from collections import Counter
import os
import tempfile
import urllib.error
import urllib.request

import pytest

from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure
from fr3d.classifiers.hydrogen_bonds import check_hydrogen_bond
from fr3d.classifiers.hydrogen_bonds import load_ideal_basepair_hydrogen_bonds
from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure
from fr3d.classifiers.NA_pairwise_interactions import focus_basepair_cutoffs
from fr3d.classifiers.NA_pairwise_interactions import load_structure
from fr3d.classifiers.class_limits_2024 import nt_nt_cutoffs
import fr3d.search.FR3D as fr3d_runner
from fr3d.search.query_processing import readQueryFromJSON


def _normalize_direct_family(label):
    """Normalize direct (non-near) 3-character interaction labels."""
    if len(label) != 3 or label[0] not in ("c", "t"):
        return None
    return label[0] + label[1:].upper()


def _count_direct_basepair_families(interaction_to_pairs):
    counts = Counter()
    for interaction, pair_list in interaction_to_pairs.items():
        family = _normalize_direct_family(interaction)
        if family:
            counts[family] += len(pair_list)
    return counts


def _fetch_online_1s72_counts():
    url = "https://rna.bgsu.edu/rna3dhub/pdb/1s72/interactions/fr3d/basepairs/tsv"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            raw = response.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError) as exc:
        pytest.skip("Could not fetch online 1S72 FR3D TSV: %s" % exc)

    counts = Counter()
    for line in raw.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        family = _normalize_direct_family(fields[1].strip())
        if family:
            counts[family] += 1

    assert len(counts) >= 10
    return counts


def _build_component(sequence, atoms):
    return Component(
        atoms=[Atom(name=name, x=x, y=y, z=z) for name, (x, y, z) in atoms.items()],
        sequence=sequence,
        polymeric=True,
    )


def test_1s72_local_matches_online_fr3d_tsv_within_tolerance():
    online_counts = _fetch_online_1s72_counts()

    structure, messages = load_structure("1S72")
    if structure is None:
        pytest.skip("Could not load 1S72: %s" % "; ".join(messages))

    interaction_to_pairs, _, _, _ = annotate_nt_nt_in_structure(
        structure,
        categories={"basepair": []},
    )
    local_counts = _count_direct_basepair_families(interaction_to_pairs)

    compared = 0
    for family, online_count in sorted(online_counts.items()):
        local_count = local_counts[family]
        delta = local_count - online_count
        assert abs(delta) <= 5, "Family %s differs by %d (local=%d, online=%d)" % (
            family,
            delta,
            local_count,
            online_count,
        )
        compared += 1

    assert compared >= 10


def test_4tna_has_about_31_basepairs():
    structure, messages = load_structure("4TNA")
    if structure is None:
        pytest.skip("Could not load 4TNA: %s" % "; ".join(messages))

    interaction_to_pairs, _, _, _ = annotate_nt_nt_in_structure(
        structure,
        categories={"basepair": []},
    )

    unique_pairs = {
        tuple(sorted((unit_id_1, unit_id_2)))
        for interaction, pair_list in interaction_to_pairs.items()
        if _normalize_direct_family(interaction)
        for unit_id_1, unit_id_2, _crossing in pair_list
    }

    assert 26 <= len(unique_pairs) <= 36


def test_fr3d_symbolic_queries_smoke_do_not_crash():

    query_paths = [
        os.path.join("fr3d", "search", "queries", "symbolic_2_AU_cWW.json"),
        os.path.join("fr3d", "search", "queries", "symbolic_2_bSS_distance_5_6.json"),
        os.path.join("fr3d", "search", "queries", "symbolic_2_tSS_LR.json"),
    ]

    for query_path in query_paths:
        query = readQueryFromJSON(query_path)
        assert "errorStatus" not in query, "Could not load query file %s" % query_path

        with tempfile.TemporaryDirectory(prefix="fr3d_query_smoke_") as temp_dir:
            query["DATAPATHUNITS"] = os.path.join(temp_dir, "units")
            query["DATAPATHPAIRS"] = os.path.join(temp_dir, "pairs")
            query["OUTPUTPATH"] = os.path.join(temp_dir, "output")
            query["searchFiles"] = []

            # Prevent any representative-set or query-unit fetches in smoke mode.
            query.pop("repSetRelease", None)
            query.pop("repSetResolution", None)
            query.pop("repSetType", None)
            query.pop("unitID", None)
            query.pop("PDB_data_file", None)

            # Keep smoke tests fast and deterministic.
            query["printInputQuery"] = False
            query["printProcessedQuery"] = False
            query["printSearchFiles"] = 0
            query["printSearchingFile"] = False
            query["printTimer"] = False

            fr3d_runner.fr3d_search(query)


def test_hbond_table_contains_expected_a_u_cww_entries_from_csv():
    # Source: fr3d/classifiers/H_bonding_Atoms_from_Isostericity_Table.csv
    #   cis W/W,A,U,N6,H61,---,O4,-,
    #   cis W/W,A,U,-,N1,---,H3,N3,
    hbond = load_ideal_basepair_hydrogen_bonds()
    entries = {(a, b, c, order) for (a, b, c, order, *_rest) in hbond["A,U"]["cWW"]}

    assert ("N6", "H61", "O4", "12") in entries
    assert ("N3", "H3", "N1", "21") in entries


def test_hbond_table_propagates_rna_to_dna_pairs():
    hbond = load_ideal_basepair_hydrogen_bonds()
    assert "A,U" in hbond
    assert "DA,DT" in hbond
    assert "cWW" in hbond["A,U"]
    assert "cWW" in hbond["DA,DT"]
    assert len(hbond["A,U"]["cWW"]) == len(hbond["DA,DT"]["cWW"])


def test_check_hydrogen_bond_uses_heavy_atom_fallback_for_h2prime():
    nt1 = _build_component(
        "A",
        {
            "O2'": (0.0, 0.0, 0.0),
            "C2'": (1.0, 0.0, 0.0),
        },
    )
    nt2 = _build_component(
        "U",
        {
            "O2": (0.0, 0.0, 3.2),
        },
    )

    atom_set = ("O2'", "H2'", "O2", "12", 0.0, 3.6, 0.0, 180.0)
    result = check_hydrogen_bond(nt1, nt2, atom_set)

    assert result["bond_checked"] is True
    assert result["length"] != result["length"]  # NaN when H coordinates are unavailable
    assert result["donor_acceptor_distance"] == pytest.approx(3.2, abs=1e-6)
    assert result["heavy_donor_acceptor_angle"] == pytest.approx(90.0, abs=1e-6)


def test_check_hydrogen_bond_handles_missing_acceptor_atom_coordinates():
    nt1 = _build_component(
        "A",
        {
            "N6": (0.0, 0.0, 0.0),
            "H61": (0.0, 0.0, 1.0),
            "C6": (-1.0, 0.0, 0.0),
        },
    )
    nt2 = _build_component("U", {"N3": (0.0, 0.0, 3.0)})

    atom_set = ("N6", "H61", "O4", "12", 0.0, 3.6, 0.0, 180.0)
    result = check_hydrogen_bond(nt1, nt2, atom_set)

    assert result["bond_checked"] is False
    assert result["bond_made"] is False
    assert result["donor_acceptor_distance"] != result["donor_acceptor_distance"]  # NaN
    assert result["badness"] >= 100.0


def test_focus_basepair_cutoffs_limits_to_requested_family():
    focused = focus_basepair_cutoffs(nt_nt_cutoffs, ["cWW"])

    assert "A,U" in focused
    interactions = set(focused["A,U"][1].keys()) | set(focused["A,U"][-1].keys())

    assert "cWW" in interactions
    assert "tWW" not in interactions
    assert "tHS" not in interactions


def test_annotate_nt_nt_in_structure_handles_empty_structure():
    structure = Structure([], pdb="EMPTY")

    interaction_to_pairs, category_to_interactions, _timer_data, pair_to_data = annotate_nt_nt_in_structure(
        structure,
        categories={"basepair": []},
    )

    assert isinstance(interaction_to_pairs, dict)
    assert isinstance(category_to_interactions, dict)
    assert isinstance(pair_to_data, dict)
    assert len(interaction_to_pairs) == 0
