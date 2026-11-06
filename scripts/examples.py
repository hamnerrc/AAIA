#!/usr/bin/env python3
"""Concrete solves the author can do by hand, one for each main finding.

For each example the scramble is the exact inverse of the baseline solution, applied to a
solved cube. Every line is replayed and asserted before it is written out.
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.argv = sys.argv[:1]
import numpy as np

from lsll import data as D
from lsll.cube import auf_tokens, inverse, join_savings, _frame_lists, perm_of, stm
from lsll.model import Config, Model
from lsll.states import LLSpace, SOLVED
from validate import Brute, block_recount, oriented, run, solved


def fmt(t):
    return " ".join(t) if t else "(nothing)"


def solve_from(B, M, lab, s1):
    a = run(lab, s1)
    s2, b = B.solve_oll(a)
    s3 = B.solve_pll(b)
    assert solved(run(b, s3))
    return s2, s3


def describe(M, B, lab, opts):
    rows = []
    for s in opts:
        s2, s3 = solve_from(B, M, lab, s)
        rows.append({"last_pair": fmt(s), "oll": fmt(s2), "pll": fmt(s3), "stm_pair": stm(s), "stm_oll": stm(s2), "stm_pll": stm(s3), "two_step_score": stm(s) + stm(s2), "total": stm(s) + stm(s2) + stm(s3)})
    return rows


def main():
    M = Model()
    B = Brute(M)
    ll = M.ll
    ex = {}
    rng = np.random.default_rng(7)
    pos = [p for p in M.pos_list if M.f2l_by_pos[p] is not None and M.ls_variants(M.f2l_cases[M.f2l_by_pos[p]]["body"])[1]]
    want = {"sledge_skips_oll", "u2_shorter_total", "two_step_loses", }
    tries = 0
    while want and tries < 200000:
        tries += 1
        p = pos[rng.integers(len(pos))]
        l = int(rng.integers(ll.N))
        lab = ll.labels[l][perm_of(M.placements[p])]
        opts = B.ls_candidates(lab, M, ("U2", "sledge", "std"))
        if len(opts) < 3:
            continue
        ll_flags = [oriented(run(lab, s)) for s in opts]
        rows = None
        if "sledge_skips_oll" in want and not ll_flags[0] and ll_flags[2]:
            rows = describe(M, B, lab, opts)
            s2, s3 = solve_from(B, M, lab, opts[0])
            if rows[2]["stm_pair"] - rows[0]["stm_pair"] <= 2 and rows[2]["total"] < rows[0]["total"]:
                ex["sledge_skips_oll"] = {"note": "The author's own example: a longer insert that skips OLL.", "scramble": fmt(inverse(opts[0] + s2 + s3)), "options": rows}
                want.discard("sledge_skips_oll")
        if "u2_shorter_total" in want and tries % 3 == 0:
            rows = rows or describe(M, B, lab, opts)
            if rows[1]["two_step_score"] < rows[0]["two_step_score"] - 2 and rows[1]["total"] < rows[0]["total"] and rows[1]["stm_pair"] > rows[0]["stm_pair"]:
                s2, s3 = solve_from(B, M, lab, opts[0])
                ex["u2_shorter_total"] = {"note": "A U2 insert that costs an extra move and gains more at OLL.", "scramble": fmt(inverse(opts[0] + s2 + s3)), "options": rows}
                want.discard("u2_shorter_total")
        if "two_step_loses" in want and tries % 3 == 0:
            rows = rows or describe(M, B, lab, opts)
            scores = [r["two_step_score"] for r in rows]
            pick = int(np.argmin(scores))
            if pick != 0 and rows[pick]["total"] >= rows[0]["total"] + 4:
                s2, s3 = solve_from(B, M, lab, opts[0])
                ex["two_step_loses"] = {"note": "The two-step rule picks the option best through OLL, and the PLL it leaves is much longer.", "scramble": fmt(inverse(opts[0] + s2 + s3)), "two_step_picks": rows[pick]["last_pair"], "options": rows}
                want.discard("two_step_loses")

    # a large cancellation: the OLL/PLL pair with the greatest saving
    best = None
    for oc in M.oll_cases:
        o = D.baseline_of(oc.algs).tokens
        for pc in M.pll_cases:
            q = D.baseline_of(pc.algs).tokens
            v = join_savings(_frame_lists(o)[1], _frame_lists(q)[0])
            if best is None or v > best[0]:
                best = (v, oc.name, o, pc.name, q)
    v, on, o, pn, q = best
    ex["largest_oll_pll_cancellation"] = {"oll_case": on, "oll_alg": fmt(o), "pll_case": pn, "pll_alg": fmt(q), "moves_cancelled": int(v), "stm_written": stm(o) + stm(q), "stm_after_merging": block_recount([o, q]),
                                          "probability_of_this_pair_at_this_angle": "1 in " + str(int(round(ll.N / ((M.oll_case == [c.name for c in M.oll_cases].index(on)).sum()) * 288 / ((M.pll_case == [c.name for c in M.pll_cases].index(pn)).sum()) * 16)))}

    (ROOT / "results" / "examples.json").write_text(json.dumps(ex, indent=1))
    print(json.dumps(ex, indent=1))
    print("remaining", want, "tries", tries)


if __name__ == "__main__":
    main()
