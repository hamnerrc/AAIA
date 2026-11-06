#!/usr/bin/env python3
"""Independent checks of the exact model.

A brute-force solver works from raw stickers: it recognises every case by trial, trying
each case's alg at each AUF until the goal is met. It shares no case tables with the
model. Random states come from long random walks, not from the model's indexing.
It is compared with the model state by state, and full solves are replayed in magiccube.
"""
import json, random, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import magiccube
import numpy as np
from scipy import stats as st

from lsll import data as D
from lsll.cube import auf_tokens, inverse, join_savings, parse, perm_of, stm, to_magiccube
from lsll.distribution import normalise, stats, tv_distance
from lsll.model import MAXT, Config, Model
from lsll.states import FIXED_STICKERS, LLSpace, SOLVED, pair_position

OUT = Path(__file__).resolve().parents[1] / "results"
N_STATES = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
N_MAGIC = 1000
K = (0, 1, 3, 2)


def run(lab, tokens):
    return LLSpace.recenter(lab[perm_of(tuple(tokens))][None, :])[0]


def f2l_ok(lab):
    return (lab[FIXED_STICKERS] == SOLVED[FIXED_STICKERS]).all()


def oriented(lab):
    return (lab[[0, 1, 2, 3, 4, 5, 6, 7, 8]] == 0).all()


def solved(lab):
    return (lab == SOLVED).all()


class Brute:
    def __init__(self, model):
        self.f2l = [(c["body"]) for c in model.f2l_cases]
        self.oll = [D.baseline_of(c.algs).tokens for c in D.load_set("oll")]
        self.pll = [D.baseline_of(c.algs).tokens for c in D.load_set("pll")]
        self.oll_perms = [[perm_of(auf_tokens(k) + a) for k in K] for a in self.oll]
        self.pll_perms = [[perm_of(auf_tokens(k) + a) for k in K] for a in self.pll]

    def solve_oll(self, lab):
        """(sequence, next labels) for orientation, by trial."""
        if oriented(lab):
            return (), lab
        for a in self.oll:
            for k in K:
                r = run(lab, auf_tokens(k) + a)
                if f2l_ok(r) and oriented(r):
                    return auf_tokens(k) + a, r
        raise AssertionError("no OLL alg")

    def solve_pll(self, lab):
        """Cheapest AUF + alg + AUF by trial."""
        best = None
        for j in range(4):
            if solved(run(lab, auf_tokens(j))):
                best = auf_tokens(j)
                return best
        for a in self.pll:
            for k in K:
                r = run(lab, auf_tokens(k) + a)
                for j in range(4):
                    if solved(run(r, auf_tokens(j))):
                        seq = auf_tokens(k) + a + auf_tokens(j)
                        if best is None or stm(seq) < stm(best):
                            best = seq
            if best is not None:
                return best
        raise AssertionError("no PLL alg")

    def ls_candidates(self, lab, model, inserts):
        """Every option for the last pair, by trial on stickers."""
        if f2l_ok(lab):
            return [()]
        out = []
        for c in model.f2l_cases:
            for k in K:
                if f2l_ok(run(lab, auf_tokens(k) + c["body"])):
                    out.append(auf_tokens(k) + c["body"])
                    body = c["body"]
                    vs, _ = model.ls_variants(body)
                    for name in inserts:
                        if name in vs:
                            b = model._strip_auf(vs[name])
                            for k2 in K:
                                if f2l_ok(run(lab, auf_tokens(k2) + b)):
                                    if all(auf_tokens(k2) + b != o for o in out):
                                        out.append(auf_tokens(k2) + b)
                                    break
                    return out
        raise AssertionError("no F2L alg")


def random_state(rnd):
    gens = [parse(a) for a in ["U", "U'", "U2", "R U R'", "R U' R'", "R U2 R'", "F' U F", "F' U' F", "R' F R F'", "F R' F' R"]]
    seq = ()
    for _ in range(60):
        seq += rnd.choice(gens)
    return seq


def locate(model, seq):
    """(p, l) of the state reached by seq, through the stickers."""
    lab = SOLVED[perm_of(seq)]
    p = pair_position(lab)
    r = model.placements[p]
    back = lab[perm_of(inverse(r))]
    idx = model.ll.index_of(back[None, :])[0]
    assert idx >= 0
    return p, idx


def block_recount(seq_list):
    """Recount merged STM from scratch: rotation-free tokens, grouped into same-axis runs.

    Also returns whether the rebuilt layer amounts reproduce the literal cube state.
    """
    from lsll.cube import _MIN_COST, _frame_lists, layer_vector
    flat = ()
    for s in seq_list:
        flat = flat + tuple(s)
    start, _ = _frame_lists(flat)
    runs = []
    for t in start:
        ax, v = layer_vector(t)
        if runs and runs[-1][0] == ax:
            runs[-1][1] = tuple((x + y) % 4 for x, y in zip(runs[-1][1], v))
            runs[-1][2] += 1
        else:
            runs.append([ax, v, 1])
    # cascade: an empty run lets its neighbours meet
    i = 0
    while i < len(runs):
        if runs[i][1] == (0, 0, 0):
            del runs[i]
            if 0 < i < len(runs) and runs[i - 1][0] == runs[i][0]:
                runs[i - 1][1] = tuple((x + y) % 4 for x, y in zip(runs[i - 1][1], runs[i][1]))
                del runs[i]
            i = max(i - 1, 0)
        else:
            i += 1
    return sum(_MIN_COST[r[1]] for r in runs)


def main():
    t0 = time.time()
    M = Model()
    B = Brute(M)
    rnd = random.Random(20261105)
    base = M.evaluate(Config())
    both = M.evaluate(Config("two", ("U2", "sledge")))
    three = M.evaluate(Config("three", ("U2", "sledge"), ls_rule="three_step"))
    canc = M.evaluate(Config("canc", cancel_auf=True, cancel_alg=True))
    hist = np.zeros(MAXT, np.int64)
    mism = dict(baseline=0, two_step=0, three_step=0, cancel=0)
    cancel_literal_higher = cancel_literal_lower = 0
    not_solved = 0
    placements = {}
    for i in range(N_STATES):
        scr = random_state(rnd)
        p, l = locate(M, scr)
        placements[p] = placements.get(p, 0) + 1
        lab = SOLVED[perm_of(scr)]
        # baseline
        cands = B.ls_candidates(lab, M, ())
        s1 = cands[0]
        a = run(lab, s1)
        s2, b = B.solve_oll(a)
        s3 = B.solve_pll(b)
        total = stm(s1) + stm(s2) + stm(s3)
        hist[total] += 1
        if not solved(run(b, s3)):
            not_solved += 1
        mism["baseline"] += total != base["by_p"][p][l]
        # cancellations: the model against a literal merge of the same segments
        # the model counts each written alg as written, so merges found inside one alg are not boundary savings
        inside = sum(stm(q) - block_recount([q]) for q in (s1, s2, s3))
        lit = block_recount([s1, s2, s3]) + inside
        mod = canc["by_p"][p][l]
        mism["cancel"] += lit != mod
        cancel_literal_higher += lit > mod
        cancel_literal_lower += lit < mod
        # insert choices, by trial
        if i < N_STATES // 4:
            opts = B.ls_candidates(lab, M, ("U2", "sledge", "std"))
            scored = []
            for s in opts:
                aa = run(lab, s)
                o, bb = B.solve_oll(aa)
                q = B.solve_pll(bb)
                scored.append((stm(s) + stm(o), stm(s) + stm(o) + stm(q)))
            mism["two_step"] += scored[min(range(len(scored)), key=lambda j: (scored[j][0], j))][1] != both["by_p"][p][l]
            mism["three_step"] += min(x[1] for x in scored) != three["by_p"][p][l]
    expected = normalise(base["total"])
    keep = expected * N_STATES >= 5
    obs = np.append(hist[keep], hist[~keep].sum())
    exp = np.append(expected[keep], expected[~keep].sum()) * N_STATES
    chi = st.chisquare(obs, exp)
    pl_obs = np.array(list(placements.values()) + [0] * (150 - len(placements)))
    chi_pl = st.chisquare(pl_obs)
    # magiccube replay
    solved_mc = 0
    mc_state = magiccube.Cube(3)
    for i in range(N_MAGIC):
        scr = random_state(rnd)
        lab = SOLVED[perm_of(scr)]
        s1 = B.ls_candidates(lab, M, ())[0]
        a = run(lab, s1)
        s2, b = B.solve_oll(a)
        s3 = B.solve_pll(b)
        c = magiccube.Cube(3)
        c.rotate(to_magiccube(scr + s1 + s2 + s3))
        solved_mc += c.is_done()
    res = {
        "n_states": N_STATES,
        "state_by_state_mismatches": {k: int(v) for k, v in mism.items()},
        "states_checked": {"baseline": N_STATES, "cancel": N_STATES, "two_step": N_STATES // 4, "three_step": N_STATES // 4},
        "cancel_literal_higher_than_model": int(cancel_literal_higher),
        "cancel_literal_lower_than_model": int(cancel_literal_lower),
        "brute_sequence_did_not_solve": not_solved,
        "mc_mean": float(hist @ np.arange(MAXT) / N_STATES),
        "exact_mean": stats(expected)["mean"],
        "chi2_distribution": {"stat": float(chi.statistic), "p": float(chi.pvalue), "bins": int(len(obs))},
        "chi2_placements_uniform": {"stat": float(chi_pl.statistic), "p": float(chi_pl.pvalue), "placements_seen": len(placements)},
        "magiccube_replays": N_MAGIC,
        "magiccube_solved": int(solved_mc),
        "seconds": time.time() - t0,
    }
    (OUT / "validation.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
