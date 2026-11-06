#!/usr/bin/env python3
"""Runs the whole study: results/results.json, results/tables/*.csv. Figures: make_figures.py."""
import csv, json, sys, time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np

from lsll import data as D
from lsll.cube import auf_tokens, join_savings, stm, _frame_lists
from lsll.distribution import convolve_dft, convolve_direct, dft, exact_mean, idft, normalise, stats, tv_distance, trim
from lsll.model import MAXT, Config, Model

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
TAB = OUT / "tables"
TH = (20, 25, 30, 35, 40)
INS = ("U2", "sledge")
ALL_STATES = 9331200


def write_csv(name, rows):
    with open(TAB / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def trimmed(h):
    return trim(np.asarray(h, dtype=float) / np.asarray(h).sum())


def main():
    t0 = time.time()
    M = Model()
    ll = M.ll
    R = {"state_space": {"ll_states": ll.N, "placements": 150, "total_states": ALL_STATES}}

    # ---------------------------------------------------- data verification
    R["data_verification"] = {
        "oll": M.verify_set(M.oll_cases, M.oll_case, ll.oriented),
        "pll": M.verify_set(M.pll_cases, M.pll_case, ll.is_u),
        "coll": M.verify_set(M.coll_cases, M.coll_case, ll.corners_u & ll.eo0),
        "zbll": M.verify_set(M.zbll_cases, M.zbll_case, ll.is_u),
    }
    print("verified", R["data_verification"], flush=True)

    # ---------------------------------------------------- baseline and its steps
    base = M.evaluate(Config())
    lv = M.ll_level(Config())
    seg_key = {v: k for k, v in M.segs.ids.items()}
    oseg, onext = lv["options"][0]
    pseg = lv["pll_seg"]

    def seg_arrays(seg):
        k = np.array([seg_key[int(s)][0] if s >= 0 else 0 for s in seg])
        mv = np.array([stm(seg_key[int(s)][1]) if s >= 0 else 0 for s in seg])
        post = np.array([seg_key[int(s)][2] if s >= 0 else 0 for s in seg])
        return k, mv, post

    ok, omv, _ = seg_arrays(oseg)
    pk, pmv, ppost = seg_arrays(pseg)
    opts = M.ls_options(Config())
    ls_k = np.array([opts[p][0][1] for p in M.pos_list])
    ls_mv = np.array([stm(opts[p][0][2]) for p in M.pos_list])

    orient = np.where(ll.oriented)[0]
    steps = {
        "f2l_alg": np.bincount(ls_mv, minlength=13) / 150,
        "auf_before_pair": np.bincount((ls_k != 0).astype(int), minlength=2) / 150,
        "oll_alg": np.bincount(omv, minlength=20) / ll.N,
        "auf_before_oll": np.bincount((ok != 0).astype(int), minlength=2) / ll.N,
        "pll_alg": np.bincount(pmv[orient], minlength=25) / len(orient),
        "auf_before_pll": np.bincount((pk[orient] != 0).astype(int), minlength=2) / len(orient),
        "auf_after_pll": np.bincount((ppost[orient] != 0).astype(int), minlength=2) / len(orient),
    }
    exact_frac = {
        "auf_before_pair_P0": str(Fraction(int((ls_k == 0).sum()), 150)),
        "auf_before_oll_P0": str(Fraction(int((ok == 0).sum()), ll.N)),
        "auf_before_pll_P0": str(Fraction(int((pk[orient] == 0).sum()), len(orient))),
        "auf_after_pll_P0": str(Fraction(int((ppost[orient] == 0).sum()), len(orient))),
    }
    R["steps"] = {k: {"pmf": trim(v).tolist(), **{m: stats(v)[m] for m in ("mean", "sd")}} for k, v in steps.items()}
    R["auf_exact"] = exact_frac
    # joint AUF after/before PLL
    both = np.zeros((2, 2), np.int64)
    np.add.at(both, ((pk[orient] != 0).astype(int), (ppost[orient] != 0).astype(int)), 1)
    R["pll_auf_joint_counts_pre_post"] = both.tolist()

    # per-case tables
    f2l_rows = []
    for i, c in enumerate(M.f2l_cases):
        ps = [p for p in M.pos_list if M.f2l_by_pos[p] == i]
        ks = [opts[p][0][1] for p in ps]
        vs, tail = M.ls_variants(c["body"])
        f2l_rows.append({"case": c["name"], "group": c["group"], "votes": c["votes"], "as_listed": c["text"], "body_used": " ".join(c["body"]), "stm": c["moves"],
                         "placements": len(ps), "probability": f"{len(ps)}/150", "P_no_auf": f"{sum(k == 0 for k in ks)}/{len(ps)}", "insert_tail": tail or "", "has_U2_and_sledge_variants": bool(tail)})
    write_csv("f2l_baseline.csv", f2l_rows)
    R["inserts"] = {"f2l_cases_with_variants": sum(r["has_U2_and_sledge_variants"] for r in f2l_rows),
                    "placements_with_variants": sum(r["placements"] for r in f2l_rows if r["has_U2_and_sledge_variants"]),
                    "variants": {r["case"]: {k: " ".join(v) for k, v in M.ls_variants(next(c["body"] for c in M.f2l_cases if c["name"] == r["case"]))[0].items()} for r in f2l_rows if r["has_U2_and_sledge_variants"]}}

    def ll_rows(cases, case_of, kk, pp=None):
        rows = []
        for i, c in enumerate(cases):
            m = np.where(case_of == i)[0]
            b = D.baseline_of(c.algs)
            row = {"case": c.name, "group": c.group, "votes": b.votes, "alg": b.text, "stm": b.moves, "score": b.score,
                   "states": len(m), "P_no_pre_auf": str(Fraction(int((kk[m] == 0).sum()), len(m)))}
            if pp is not None:
                row["P_no_post_auf"] = str(Fraction(int((pp[m] == 0).sum()), len(m)))
            alts = D.alternates_of(c.algs, b)
            row["qualifying_alternates"] = len(alts)
            row["mean_alt_stm_minus_base"] = round(float(np.mean([a.moves - b.moves for a in alts])), 2) if alts else ""
            rows.append(row)
        return rows

    orow = ll_rows(M.oll_cases, M.oll_case, ok)
    prow = ll_rows(M.pll_cases, M.pll_case, pk, ppost)
    write_csv("oll_baseline.csv", orow)
    write_csv("pll_baseline.csv", prow)

    # alternates: is faster also shorter?
    alt = {}
    for name, cases in (("oll", M.oll_cases), ("pll", M.pll_cases)):
        d = []
        for c in cases:
            b = D.baseline_of(c.algs)
            d += [a.moves - b.moves for a in D.alternates_of(c.algs, b)]
        d = np.array(d)
        alt[name] = {"cases_with_alternates": sum(bool(D.alternates_of(c.algs, D.baseline_of(c.algs))) for c in cases), "alternates": len(d),
                     "mean_stm_minus_baseline": float(d.mean()), "share_longer": float((d > 0).mean()), "share_same": float((d == 0).mean()), "share_shorter": float((d < 0).mean())}
    R["alternates_faster_vs_shorter"] = alt

    # ---------------------------------------------------- independence: exact integer tests
    c1, c2, _ = M.ll_cost(Config(), None)
    # OLL stage vs PLL stage on the 62,208 LL states
    joint = np.zeros((c1.max() + 1, c2.max() + 1), np.int64)
    np.add.at(joint, (c1, c2), 1)
    prod = np.outer(joint.sum(1), joint.sum(0))
    oll_pll_exact = bool((joint * ll.N == prod).all())
    # LS vs LL: every placement must see the same LL distribution
    tables = M.ls_tables(Config())
    ref = np.bincount(c1 + c2, minlength=60)
    ls_exact = all((np.bincount((c1 + c2)[o[0][2]], minlength=60) == ref).all() for o in tables.values())
    # inside the OLL block: is the AUF independent of the alg length?
    j2 = np.zeros((2, omv.max() + 1), np.int64)
    np.add.at(j2, ((ok != 0).astype(int), omv), 1)
    auf_alg_oll = bool((j2 * ll.N == np.outer(j2.sum(1), j2.sum(0))).all())
    j3 = np.zeros((2, pmv.max() + 1), np.int64)
    np.add.at(j3, ((pk[orient] != 0).astype(int), pmv[orient]), 1)
    auf_alg_pll = bool((j3 * len(orient) == np.outer(j3.sum(1), j3.sum(0))).all())
    jj = np.zeros((2, pmv.max() + 1), np.int64)
    np.add.at(jj, ((ppost[orient] != 0).astype(int), pmv[orient]), 1)
    post_alg = bool((jj * len(orient) == np.outer(jj.sum(1), jj.sum(0))).all())
    R["independence_exact_tests"] = {
        "LS_block_vs_LL_block": ls_exact, "OLL_block_vs_PLL_block": oll_pll_exact,
        "AUF_before_OLL_vs_OLL_length": auf_alg_oll, "AUF_before_PLL_vs_PLL_length": auf_alg_pll, "AUF_after_PLL_vs_PLL_length": post_alg,
    }

    # ---------------------------------------------------- three ways to compute the baseline
    enumerated = trimmed(base["total"])
    blocks = [trimmed(base["ls"]), trimmed(base["oll"]), trimmed(base["pll"])]
    by_direct = convolve_direct(*blocks)
    by_dft = convolve_dft(*blocks)
    flat_auf = np.array([0.25, 0.75])
    bodies = [trim(steps["f2l_alg"]), trim(steps["oll_alg"]), trim(steps["pll_alg"])]
    flat = convolve_dft(flat_auf, bodies[0], flat_auf, bodies[1], flat_auf, bodies[2], flat_auf)
    exact_auf_flat = convolve_dft(trim(steps["auf_before_pair"]), bodies[0], trim(steps["auf_before_oll"]), bodies[1], trim(steps["auf_before_pll"]), bodies[2], trim(steps["auf_after_pll"]))
    R["baseline"] = {
        "pmf": enumerated.tolist(),
        "counts": base["total"][: len(enumerated)].tolist(),
        "exact_mean_fraction": str(exact_mean(base["total"])),
        "stats": stats(enumerated, TH),
        "stage_means": {k: stats(trimmed(base[k]))["mean"] for k in ("ls", "oll", "pll", "ll", "lsoll")},
        "convolution_vs_enumeration": {
            "tv_block_direct": tv_distance(by_direct, enumerated), "tv_block_dft": tv_distance(by_dft, enumerated),
            "tv_direct_dft": tv_distance(by_direct, by_dft),
            "tv_flat_auf_model": tv_distance(flat, enumerated), "flat_auf_stats": stats(flat, TH),
            "tv_exact_auf_but_steps_separate": tv_distance(exact_auf_flat, enumerated),
            "dft_padding_N": len(by_dft), "flat_pmf": flat.tolist(),
        },
    }

    R["block_pmfs"] = {"ls": blocks[0].tolist(), "oll": blocks[1].tolist(), "pll": blocks[2].tolist()}
    n_checked = sum(v["algs_checked"] for v in R["data_verification"].values())
    n_failed = sum(v["failed"] for v in R["data_verification"].values())
    R["data_verification"]["summary"] = {"algs_checked": n_checked, "failed": n_failed, "share_failed": n_failed / n_checked}

    # the alternates, listed case by case
    alt_rows = []
    for setname, cases in (("OLL", M.oll_cases), ("PLL", M.pll_cases)):
        for c in cases:
            b = D.baseline_of(c.algs)
            for a in sorted(D.alternates_of(c.algs, b), key=lambda q: q.score):
                alt_rows.append({"set": setname, "case": c.name, "baseline": b.text, "baseline_stm": b.moves, "baseline_score": b.score, "alternate": a.text, "alternate_stm": a.moves, "alternate_score": a.score, "stm_minus_baseline": a.moves - b.moves, "source": a.src})
    write_csv("alternates.csv", alt_rows)

    # DFT worked pieces
    a = np.array([0.25, 0.75, 0, 0])
    A = dft(a)
    R["dft_hand_example"] = {
        "input": a.tolist(), "dft": [[float(z.real), float(z.imag)] for z in A],
        "dft_cubed": [[float(z.real), float(z.imag)] for z in A**3], "inverse": idft(A**3).real.tolist(),
        "binomial_3_0.75": [0.015625, 0.140625, 0.421875, 0.421875],
    }
    need = len(by_dft)
    wraps = {}
    for n in (20, 32, 40, need):
        w = convolve_dft(*blocks, n=n)
        wraps[n] = {"pmf": w.tolist(), "tv_against_true": tv_distance(w, enumerated), "mean_as_read": float(w @ np.arange(n)), "max_prob_index": int(np.argmax(w))}
    R["wraparound"] = {"true_support_max": len(enumerated) - 1, "required_N": need, "tests": wraps}

    # ---------------------------------------------------- techniques
    tech = [
        ("Baseline", Config(), 0),
        ("U2 insert, always", Config("x", ("U2",), ls_rule="always_U2"), 0),
        ("Sledge insert, always", Config("x", ("sledge",), ls_rule="always_sledge"), 0),
        ("U2 insert, two-step", Config("x", ("U2",)), 0),
        ("Sledge insert, two-step", Config("x", ("sledge",)), 0),
        ("U2 + sledge, two-step", Config("x", INS), 0),
        ("U2 + sledge, full look-ahead", Config("x", INS, ls_rule="three_step"), 0),
        ("Cancellations", Config("x", cancel_auf=True, cancel_alg=True), 0),
        ("Angles OLL, speed-gated", Config("x", oll_angles="speed"), None),
        ("Angles PLL, speed-gated", Config("x", pll_angles="speed"), None),
        ("Angles both, speed-gated", Config("x", oll_angles="speed", pll_angles="speed"), None),
        ("Angles both, speed-gated, no awkward", Config("x", oll_angles="speed", pll_angles="speed", exclude_awk=True), None),
        ("Angles both, fewest moves", Config("x", oll_angles="min", pll_angles="min"), None),
        ("Angles OLL, two-step", Config("x", oll_angles="two_step"), None),
        ("COLL, always when available", Config("x", coll="always"), 40),
        ("COLL, two-step", Config("x", coll="two_step"), 40),
        ("OLLCP, two-step", Config("x", coll="ollcp_two_step"), 329),
        ("ZBLL only (no edge control)", Config("x", ll_mode="zb"), 472),
        ("ZBLS + ZBLL", Config("x", ll_mode="zb", ls_mode="zbls"), 302 + 472 - 41),
        ("ZBLS (shortest) + ZBLL", Config("x", ll_mode="zb", ls_mode="zbls", zb_pick="moves"), 302 + 472 - 41),
        ("ZBLS + ZBLL + cancellations", Config("x", ll_mode="zb", ls_mode="zbls", cancel_auf=True, cancel_alg=True), 302 + 472 - 41),
        ("OLS + PLL", Config("x", ls_mode="ols"), 8054 - 41),
        ("1LLL", Config("x", ll_mode="1lll"), 3915 - 57 - 21),
        ("All CFOP techniques, speed-gated", Config("x", INS, cancel_auf=True, cancel_alg=True, oll_angles="speed", pll_angles="speed", coll="two_step"), None),
        ("All CFOP techniques, fewest moves", Config("x", INS, cancel_auf=True, cancel_alg=True, oll_angles="min", pll_angles="min", coll="two_step"), None),
    ]

    def used_alternates(cfg):
        extra = 0
        if cfg.oll_angles != "off":
            rule = "speed" if cfg.oll_angles == "two_step" else cfg.oll_angles
            _, _, a_ = M.ll_table(M.oll_cases, M.oll_case, ll.oriented, False, rule, cfg.exclude_awk)
            for i, c in enumerate(M.oll_cases):
                u = np.unique(a_[M.oll_case == i])
                extra += sum(c.algs[j] != D.baseline_of(c.algs) for j in u)
        if cfg.pll_angles != "off":
            _, _, a_ = M.ll_table(M.pll_cases, M.pll_case, ll.is_u, True, cfg.pll_angles, cfg.exclude_awk)
            for i, c in enumerate(M.pll_cases):
                u = np.unique(a_[M.pll_case == i])
                extra += sum(c.algs[j] != D.baseline_of(c.algs) for j in u)
        return extra

    x = np.arange(-MAXT, MAXT + 1)
    results, rows = {}, []
    for name, cfg, added in tech:
        r = M.evaluate(cfg, baseline=base)
        d = r["diff"]
        p = trimmed(r["total"])
        s = stats(p, TH)
        blocks_t = [trimmed(r["ls"]), trimmed(r["oll"]), trimmed(r["pll"])]
        conv = convolve_dft(*blocks_t)
        cs = stats(conv, TH)
        if added is None:
            added = used_alternates(cfg) + (40 if cfg.coll != "off" else 0)
        total_algs = 119 + added
        saved = float(d @ x / d.sum())
        results[name] = {
            "pmf": p.tolist(), "exact_mean": str(exact_mean(r["total"])), **s,
            "saved": saved, "P_fewer": float(d[MAXT + 1:].sum() / d.sum()), "P_same": float(d[MAXT] / d.sum()), "P_more": float(d[:MAXT].sum() / d.sum()),
            "saved_pmf": {int(k): float(v / d.sum()) for k, v in zip(x, d) if v},
            "usage": {k: v / ALL_STATES for k, v in r["usage"].items()}, "ll_use": (r["ll_use"] / ALL_STATES).tolist(),
            "stage_means": {k: stats(trimmed(r[k]))["mean"] for k in ("ls", "oll", "pll", "ll", "lsoll")},
            "stage_pmf": {k: trimmed(r[k]).tolist() for k in ("ls", "lsoll", "ll")},
            "conv_of_marginals": {"tv": tv_distance(conv, p), "sd": cs["sd"], "P<=25": cs["P(T<=25)"], "P<=30": cs["P(T<=30)"], "mean": cs["mean"], "pmf": conv.tolist()},
            "algs_added": added, "algs_total": total_algs, "saved_per_100_added": (saved / added * 100) if added else None,
        }
        rows.append({"technique": name, "mean": round(s["mean"], 4), "sd": round(s["sd"], 3), "moves_saved": round(saved, 4), "P_fewer": round(results[name]["P_fewer"], 4), "P_same": round(results[name]["P_same"], 4), "P_more": round(results[name]["P_more"], 4),
                     "P(T<=25)": round(s["P(T<=25)"], 4), "P(T<=30)": round(s["P(T<=30)"], 4), "LS": round(results[name]["stage_means"]["ls"], 3), "LS+OLL": round(results[name]["stage_means"]["lsoll"], 3), "LL": round(results[name]["stage_means"]["ll"], 3),
                     "algs_added": added, "saved_per_100_algs": round(results[name]["saved_per_100_added"], 3) if added else "", "TV_convolution_error": f"{results[name]['conv_of_marginals']['tv']:.2e}"})
        print(name, rows[-1]["mean"], rows[-1]["moves_saved"], f"{time.time()-t0:.0f}s", flush=True)
    R["techniques"] = results
    cfop = ["U2 + sledge, two-step", "U2 + sledge, full look-ahead", "Cancellations", "Angles both, fewest moves", "COLL, two-step", "OLLCP, two-step", "All CFOP techniques, fewest moves"]
    worst = max(cfop, key=lambda k: results[k]["conv_of_marginals"]["tv"])
    R["independence_failure_example"] = {"technique": worst, "tv": results[worst]["conv_of_marginals"]["tv"], "conv_pmf": results[worst]["conv_of_marginals"]["pmf"]}
    write_csv("techniques.csv", rows)

    # ---------------------------------------------------- greedy gap
    two = M.evaluate(Config("x", INS))
    three = M.evaluate(Config("x", INS, ls_rule="three_step"), baseline=two)
    d = three["diff"]
    R["greedy_gap"] = {"mean_two_minus_full": float(d @ x / d.sum()), "P_two_step_suboptimal": float(d[MAXT + 1:].sum() / d.sum()), "P_two_step_more_than_baseline": results["U2 + sledge, two-step"]["P_more"], "pmf": {int(k): float(v / d.sum()) for k, v in zip(x, d) if v}}

    # ---------------------------------------------------- mixture demonstration (COLL always)
    r = M.evaluate(Config("x", coll="always"))
    c1m, c2m, used = M.ll_cost(Config("x", coll="always"), None)
    mix = np.zeros(80)
    parts = {}
    for label, mask in (("edges not oriented: OLL then PLL", ~ll.eo0), ("edges oriented, corners not: COLL then EPLL", ll.eo0 & ~ll.co0), ("oriented: PLL only", ll.oriented)):
        a_ = normalise(np.bincount(c1m[mask]))
        b_ = normalise(np.bincount(c2m[mask]))
        w = mask.sum() / ll.N
        q = convolve_dft(a_, b_)
        mix[: len(q)] += w * q
        j = np.zeros((c1m[mask].max() + 1, c2m[mask].max() + 1), np.int64)
        np.add.at(j, (c1m[mask], c2m[mask]), 1)
        parts[label] = {"weight": float(w), "step_independent_within_class": bool((j * mask.sum() == np.outer(j.sum(1), j.sum(0))).all()), "mean": float(np.bincount(c1m[mask] + c2m[mask]) @ np.arange(c1m[mask].max() + c2m[mask].max() + 1) / mask.sum())}
    mix = convolve_dft(trimmed(base["ls"]), trim(mix))
    R["mixture_demo"] = {"classes": parts, "tv_mix_vs_enumeration": tv_distance(mix, trimmed(r["total"])), "tv_plain_conv_vs_enumeration": results["COLL, always when available"]["conv_of_marginals"]["tv"]}

    # ---------------------------------------------------- cancellations: bilinear expected savings
    canc_cfg = Config("x", cancel_auf=True, cancel_alg=True)
    ls_w = {}
    for p, opt in tables.items():
        s = opt[0][1]
        ls_w[s] = ls_w.get(s, 0) + 1 / 150
    o_w = {int(s): int((oseg == s).sum()) / ll.N for s in np.unique(oseg[oseg >= 0])}
    p_all = np.where(ll.oriented)[0]
    p_w = {int(s): int((pseg[p_all] == s).sum()) / len(p_all) for s in np.unique(pseg[p_all])}
    skip = int(ll.oriented.sum()) / ll.N

    def J(a, b):
        return join_savings(_frame_lists(M.segs.tokens[a])[1], _frame_lists(M.segs.tokens[b])[0])

    pairs = {}
    e1 = 0.0
    for a, wa in ls_w.items():
        if a is None:
            continue
        for b, wb in o_w.items():
            v = J(a, b)
            e1 += wa * wb * v
            if v:
                key = (_frame_lists(M.segs.tokens[a])[1][-1], _frame_lists(M.segs.tokens[b])[0][0])
                pairs[key] = pairs.get(key, 0) + wa * wb * v
    e2 = 0.0
    pairs2 = {}
    for a, wa in o_w.items():
        for b, wb in p_w.items():
            v = J(a, b)
            e2 += wa * wb * v
            if v:
                key = (_frame_lists(M.segs.tokens[a])[1][-1], _frame_lists(M.segs.tokens[b])[0][0] if M.segs.tokens[b] else "")
                pairs2[key] = pairs2.get(key, 0) + wa * wb * v
    # LS straight into PLL when OLL skips
    e3 = sum(wa * skip * wb * J(a, b) for a, wa in ls_w.items() if a is not None for b, wb in p_w.items())
    enumerated_saving = results["Cancellations"]["saved"]
    R["cancellations"] = {
        "expected_LS_OLL": e1, "expected_OLL_PLL_given_oll": e2, "expected_LS_PLL_on_oll_skip": e3,
        "bilinear_total": e1 * (1) + e2 + e3, "enumerated_total": enumerated_saving,
        "top_LS_OLL": sorted(([f"{k[0]} + {k[1]}", v] for k, v in pairs.items()), key=lambda t: -t[1])[:10],
        "top_OLL_PLL": sorted(([f"{k[0]} + {k[1]}", v] for k, v in pairs2.items()), key=lambda t: -t[1])[:10],
    }
    canc = M.evaluate(canc_cfg, baseline=base)
    saved_hist = canc["diff"]
    R["cancellations"]["saved_pmf"] = {int(k): float(v / saved_hist.sum()) for k, v in zip(x, saved_hist) if v}

    # ---------------------------------------------------- factorial interactions
    fact = []
    means = {}
    for rule in ("speed", "min"):
        for I in (0, 1):
            for C in (0, 1):
                for A in (0, 1):
                    for L in (0, 1):
                        if rule == "min" and not A:
                            continue
                        key = (rule if A else "off", I, C, L)
                        if key in means:
                            continue
                        cfg = Config("x", INS if I else (), cancel_auf=bool(C), cancel_alg=bool(C), oll_angles=rule if A else "off", pll_angles=rule if A else "off", coll="two_step" if L else "off")
                        rr = M.evaluate(cfg)
                        means[key] = stats(trimmed(rr["total"]))["mean"]
    b0 = means[("off", 0, 0, 0)]
    for key, m in means.items():
        ang, I, C, L = key
        add = 0.0
        if ang != "off":
            add += b0 - means[(ang, 0, 0, 0)]
        add += I * (b0 - means[("off", 1, 0, 0)]) + C * (b0 - means[("off", 0, 1, 0)]) + L * (b0 - means[("off", 0, 0, 1)])
        fact.append({"angles": ang, "inserts_two_step": I, "cancellations": C, "coll_two_step": L, "mean": round(m, 4), "saved": round(b0 - m, 4), "sum_of_individual_savings": round(add, 4), "interaction": round((b0 - m) - add, 4)})
    write_csv("factorial.csv", fact)
    R["factorial"] = fact

    (OUT / "results.json").write_text(json.dumps(R, indent=1))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
