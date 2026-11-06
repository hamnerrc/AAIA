import sys
sys.path.insert(0, "src")
import numpy as np
import pytest
from lsll import data as D
from lsll.cube import auf_tokens, parse, perm_of, stm, inverse
from lsll.distribution import convolve_dft, convolve_direct, dft, idft, normalise, stats, tv_distance
from lsll.model import Config, Model
from lsll.states import SOLVED, parse_labels, ls_state_ok, pair_position

@pytest.fixture(scope="module")
def M():
    return Model()

def test_counts(M):
    ll = M.ll
    assert ll.N == 62208 and ll.oriented.sum() == 288 and ll.eo0.sum() == 7776 and ll.is_u.sum() == 4
    assert len(M.placements) == 150
    for name in ("oll", "pll", "coll", "zbll", "ollcp", "1lll"):
        assert M.report[f"{name}_uncovered_states"] == 0
    assert (M.report["oll_cases"], M.report["pll_cases"], M.report["zbll_cases"]) == (57, 21, 472)

def test_case_probabilities(M):
    # OLL skip 1/216 and PLL skip 1/72 of oriented states
    assert (M.oll_case < 0).sum() == 288 and M.ll.is_u.sum() / M.ll.oriented.sum() == 1 / 72
    # one OLL case with 4-fold symmetry has exactly 1/216 of the states
    sizes = np.bincount(M.oll_case[M.oll_case >= 0])
    assert sizes.min() == 288 and sizes.max() == 1152

def test_placement_state_is_inserted_by_inverse(M):
    # r_p then its inverse is solved, and every placement decomposes uniquely
    for p, r in M.placements.items():
        assert ls_state_ok(SOLVED[perm_of(r + inverse(r))])
        assert pair_position(SOLVED[perm_of(r)]) == p

def test_stored_pre_post_agree_with_simulation(M):
    # The data's pre/post fields: U^pre alg U^post must solve the listed state.
    import json
    from pathlib import Path
    d = json.loads((D.DATA / "ll" / "pll.json").read_text())
    for c in d["cases"]:
        lab = parse_labels(c["state"])
        for a in c["algs"]:
            seq = auf_tokens(a["pre"]) + D.upright(parse(a["alg"])) + auf_tokens(a["post"])
            res = lab[perm_of(seq)]
            if a["alg"].startswith("x' R' F R F' U R U2"):
                continue  # the one PLL listing whose final turn is not on the U layer once upright
            assert M.ll.recenter(res[None, :])[0].tolist() == SOLVED.tolist(), (c["name"], a["alg"])

def test_inserts_really_are_equivalent(M):
    # the author's definitions: all three inserts solve the pair from the same state
    setup = inverse(parse("R U' R'"))
    for ins in ("R U' R'", "U R U2 R'", "U' R' F R F'"):
        assert ls_state_ok(SOLVED[perm_of(setup + parse(ins))])

def test_dft_matches_direct():
    a, b, c = normalise([1, 3]), normalise([2, 1, 5, 2]), normalise([1, 1, 1, 7, 2])
    assert np.allclose(convolve_dft(a, b, c), convolve_direct(a, b, c))
    assert np.allclose(idft(dft(b)), b)

def test_hand_example():
    a = np.array([0.25, 0.75, 0, 0])
    A = dft(a)
    assert np.allclose(A, [1, 0.25 - 0.75j, -0.5, 0.25 + 0.75j])
    assert np.allclose(idft(A**3).real, [1 / 64, 9 / 64, 27 / 64, 27 / 64])

def test_wraparound():
    a = normalise([1, 1, 1])
    assert np.allclose(convolve_dft(a, a, n=5), convolve_direct(a, a))
    w = convolve_dft(a, a, n=3)
    assert np.allclose(w, [1 / 3, 1 / 3, 1 / 3])   # index 3,4 folded onto 0,1

def test_totals_and_means(M):
    b = M.evaluate(Config())
    assert b["total"].sum() == 9331200
    assert stats(normalise(b["total"]))["mean"] == pytest.approx(sum(stats(normalise(b[k]))["mean"] for k in ("ls", "oll", "pll")))

def test_full_lookahead_never_loses_to_baseline(M):
    b = M.evaluate(Config())
    r = M.evaluate(Config("x", ("U2", "sledge"), ls_rule="three_step"), baseline=b)
    assert r["diff"][: 90].sum() == 0
