#!/usr/bin/env python3
"""Checks the model's assumptions against the reco.nz reconstructions.

Every solve is replayed. The last unsolved slot is identified from stickers, then the
state is viewed from the y rotation that makes that slot FR. Solves whose last slot was
already FR are reported separately, as the project is defined on FR.
"""
import json, re, sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
from scipy import stats as st

from lsll import data as D
from lsll.cube import parse, perm_of, stm
from lsll.model import Config, Model
from lsll.states import EDGE_F, CORNER_F, LLSpace, SOLVED, pair_position

SLOT = {"FR": (4, 8), "FL": (5, 9), "BL": (6, 10), "BR": (7, 11)}
Y_TO_FR = {"FR": (), "FL": ("y'",), "BL": ("y2",), "BR": ("y",)}


def slot_solved(lab, slot):
    c, e = SLOT[slot]
    return all(lab[i] == SOLVED[i] for i in CORNER_F[c]) and all(lab[i] == SOLVED[i] for i in EDGE_F[e])


def parse_file(path):
    solves = []
    cur = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("# reco.nz"):
            m = re.search(r"solve/(\d+) \| ([^|]+)\|.*\[method: (\w+)\]", line)
            cur = {"id": int(m.group(1)), "solver": m.group(2).strip(), "method": m.group(3), "lines": []}
            solves.append(cur)
        elif line.startswith("#") or not line.strip():
            continue
        elif cur is not None:
            cur["lines"].append(line)
    return solves


def main():
    M = Model()
    ll = M.ll
    solves = parse_file(ROOT / "reco_lsll.txt")
    recs = []
    errors = Counter()
    for s in solves:
        try:
            scr = parse(s["lines"][0])
            parts = [l.split("//") for l in s["lines"][1:]]
            insp = next((p[0] for p in parts if "inspection" in p[1]), "")
            path = next(p[0] for p in parts if "3-pair" in p[1])
            ls = parts[-1][0]
            lab = SOLVED[perm_of(scr + parse(insp) + parse(path))]
            lab = LLSpace.recenter(lab[None, :])[0]
            open_slots = [q for q in SLOT if not slot_solved(lab, q)]
            if len(open_slots) != 1:
                errors["not exactly one open slot"] += 1
                continue
            slot = open_slots[0]
            y = Y_TO_FR[slot]
            view = LLSpace.recenter(lab[perm_of(tuple(y))][None, :])[0]
            p = pair_position(view)
            after = LLSpace.recenter(lab[perm_of(parse(ls))][None, :])[0]
            idx = ll.index_of(after[None, :])[0]
            if idx < 0:
                errors["segment did not finish F2L"] += 1
                continue
            recs.append({"id": s["id"], "method": s["method"], "slot": slot, "pos": p, "case": M.f2l_by_pos[p], "ls_stm": stm(parse(ls)), "ll_idx": int(idx)})
        except Exception as e:
            errors[type(e).__name__] += 1
    out = {"solves_in_file": len(solves), "usable": len(recs), "errors": dict(errors)}
    out["last_slot_counts"] = dict(Counter(r["slot"] for r in recs))
    out["method_counts"] = dict(Counter(r["method"] for r in recs))

    cf = [r for r in recs if r["method"] == "OTHER"]
    fr = [r for r in cf if r["slot"] == "FR"]

    def uniform_test(rows):
        """F2L case frequencies against the model's uniform 150 placements."""
        counts = np.zeros(42)
        for r in rows:
            counts[41 if r["case"] is None else r["case"]] += 1
        exp = np.array([sum(M.f2l_by_pos[p] == i for p in M.pos_list) for i in range(41)] + [1]) / 150 * len(rows)
        groups = {"pair solved or in slot": [], "other": []}
        chi = st.chisquare(counts, exp)
        names = [c["name"] for c in M.f2l_cases] + ["Pair solved"]
        order = np.argsort(-(counts - exp))
        return {"n": len(rows), "chi2": float(chi.statistic), "p": float(chi.pvalue), "pair_already_solved": int(counts[41]),
                "pair_already_solved_expected": float(exp[41]),
                "most_over": [[names[i], int(counts[i]), round(float(exp[i]), 1)] for i in order[:5]],
                "most_under": [[names[i], int(counts[i]), round(float(exp[i]), 1)] for i in order[-5:]]}

    out["placement_uniformity_all_cfop"] = uniform_test(cf)
    out["placement_uniformity_FR_only"] = uniform_test(fr)

    def edge_test(rows):
        n = len(rows)
        eo = sum(bool(ll.eo0[r["ll_idx"]]) for r in rows)
        co = sum(bool(ll.co0[r["ll_idx"]]) for r in rows)
        oll = sum(bool(ll.oriented[r["ll_idx"]]) for r in rows)
        return {"n": n, "edges_oriented": eo, "edges_oriented_share": eo / n, "model_share": 1 / 8, "binomial_p": float(st.binomtest(eo, n, 1 / 8).pvalue),
                "corners_oriented_share": co / n, "corners_model_share": 1 / 27, "corners_binomial_p": float(st.binomtest(co, n, 1 / 27).pvalue), "oll_skips": oll, "oll_skips_expected": n / 216}

    out["after_ls_cfop_all"] = edge_test(cf)
    out["after_ls_cfop_FR_only"] = edge_test(fr)
    out["after_ls_method_zbls"] = edge_test([r for r in recs if r["method"] == "ZBLS"])

    def lsmove(rows):
        v = np.array([r["ls_stm"] for r in rows])
        return {"n": len(v), "mean": float(v.mean()), "sd": float(v.std(ddof=1)), "se": float(v.std(ddof=1) / len(v) ** 0.5), "pmf": (np.bincount(v, minlength=20) / len(v)).tolist()}

    base = M.evaluate(Config())
    from lsll.distribution import normalise, stats
    model_ls = normalise(base["ls"])
    out["last_pair_stm"] = {"model_mean_including_auf": stats(model_ls)["mean"], "model_sd": stats(model_ls)["sd"], "model_pmf": model_ls[:20].tolist(), "real_cfop_all": lsmove(cf), "real_cfop_FR_only": lsmove(fr)}

    # ------------------------------------------------- last-layer rows from pro_ll.json
    rows = json.loads((D.DATA / "ll" / "pro_ll.json").read_text())["rows"]
    res = {}
    for name, cases, case_of in (("oll", M.oll_cases, M.oll_case), ("pll", M.pll_cases, M.pll_case)):
        sel = [r for r in rows if r["set"] == name]
        moves = np.array([stm(parse(r["body"])) for r in sel])
        # model weights over cases, and body lengths
        base_alg = {D.baseline_of(c.algs).tokens for c in cases}
        same = sum(parse(r["body"]) in base_alg for r in sel)
        w = np.array([(case_of == i).sum() for i in range(len(cases))], float)
        w /= w.sum()
        mv = np.array([D.baseline_of(c.algs).moves for c in cases])
        # case frequencies, through the keys stored with each row
        counts = np.zeros(len(cases))
        unmatched = 0
        lookup = {}
        for l in np.nonzero(case_of >= 0)[0]:
            if name == "oll":
                k = "".join(map(str, ll.co[l])) + "|" + "".join(map(str, ll.eo[l]))
            else:
                k = "".join(map(str, ll.cp[l])) + "|" + "".join(map(str, ll.ep[l]))
            lookup.setdefault(k, case_of[l])
        for r in sel:
            c = lookup.get(r["key"])
            if c is None:
                unmatched += 1
            else:
                counts[c] += 1
        chi = st.chisquare(counts, w * counts.sum())
        d = counts / counts.sum() - w
        order = np.argsort(-d)
        res[name] = {"rows": len(sel), "real_mean_stm": float(moves.mean()), "real_se": float(moves.std(ddof=1) / len(moves) ** 0.5),
                     "model_mean_stm_over_non_skip": float(w @ mv), "baseline_algs_at_real_case_frequencies": float((counts / counts.sum()) @ mv),
                     "share_exactly_baseline_alg": same / len(sel), "unmatched_keys": unmatched, "chi2_case_uniform": float(chi.statistic), "p": float(chi.pvalue),
                     "most_over": [[cases[i].name, round(float(d[i]), 4)] for i in order[:5]], "most_under": [[cases[i].name, round(float(d[i]), 4)] for i in order[-5:]]}
    out["last_layer"] = res
    all_sets = Counter(r["set"] for r in rows)
    out["last_layer_set_counts"] = dict(all_sets)
    # edge-oriented share among first LL steps of solves, by what solvers did
    first = {}
    for r in rows:
        first.setdefault(r["id"], set()).add(r["set"])
    out["solves_by_ll_method"] = {"+".join(sorted(k)): v for k, v in Counter(tuple(sorted(v)) for v in first.values()).most_common()}

    # OLL edges pattern: share of OLL rows with all edges oriented
    ol = [r for r in rows if r["set"] == "oll"]
    eo = sum(r["key"].split("|")[1] == "0000" for r in ol)
    out["oll_rows_with_oriented_edges"] = {"n": len(ol), "share": eo / len(ol), "model_share_given_non_skip": 26 / 215}

    (ROOT / "results" / "real_solves.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:6000])


if __name__ == "__main__":
    main()
