#!/usr/bin/env python3
"""Size of the one thing the pairwise cancellation model misses.

The model merges each pair of neighbouring segments. If a middle alg cancels away
completely, its two neighbours meet and can merge further. Here every (last pair, OLL,
PLL) segment triple where that could happen is recounted as one sequence.
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np

from lsll.cube import _MIN_COST, _frame_lists, join_savings, layer_vector, stm
from lsll.model import Config, Model


def blocks(tokens):
    out = []
    for t in tokens:
        ax, v = layer_vector(t)
        if out and out[-1][0] == ax:
            out[-1][1] = tuple((x + y) % 4 for x, y in zip(out[-1][1], v))
            out[-1][2] += 1
        else:
            out.append([ax, v, 1, False])
    return out


def absorb(stack, new):
    """Append blocks to the stack; only the first may merge, and an empty result cascades."""
    new = [list(b) for b in new]
    while stack and new and stack[-1][0] == new[0][0]:
        a, b = stack.pop(), new.pop(0)
        v = tuple((x + y) % 4 for x, y in zip(a[1], b[1]))
        if v != (0, 0, 0):
            stack.append([a[0], v, a[2] + b[2], True])
            break
    stack.extend(new)


def cost(stack):
    return sum(min(_MIN_COST[b[1]], b[2]) if b[3] else b[2] for b in stack)


def split_frames(a, b, c):
    """Rotation-free spelling of a + b + c in one frame, returned as three parts."""
    whole = _frame_lists(a + b + c)[0]
    na, nb = stm(a), stm(b)
    return whole[:na], whole[na:na + nb], whole[na + nb:]


def main():
    M = Model()
    ll = M.ll
    cfg = Config()
    L = M.ll_level(cfg)
    tables = M.ls_tables(cfg)
    oseg, _ = L["options"][0]
    pseg = L["pll_seg"]
    ori = np.where(ll.oriented)[0]
    w_ls = {}
    for p, opt in tables.items():
        s = opt[0][1]
        if s is not None:
            w_ls[s] = w_ls.get(s, 0) + 1 / 150
    w_o = {int(s): float((oseg == s).sum()) / ll.N for s in np.unique(oseg[oseg >= 0])}
    w_p = {int(s): float((pseg[ori] == s).sum()) / len(ori) for s in np.unique(pseg[ori])}
    tok = M.segs.tokens
    fr = {s: _frame_lists(tok[s]) for s in list(w_ls) + list(w_o) + list(w_p)}
    J = lambda a, b: join_savings(fr[a][1], fr[b][0])
    checked = differ = 0
    extra = 0.0
    prob = 0.0
    worst = None
    for o, wo in w_o.items():
        n = stm(tok[o])
        for q, wq in w_p.items():
            jq = J(o, q)
            for a, wa in w_ls.items():
                ja = J(a, o)
                if ja + jq < n - 1:  # the alg cannot be used up from both ends
                    continue
                checked += 1
                pa, pb, pc = split_frames(tok[a], tok[o], tok[q])
                st = blocks(pa)
                absorb(st, blocks(pb))
                absorb(st, blocks(pc))
                pair = stm(tok[a]) + n + stm(tok[q]) - ja - jq
                true = cost(st)
                if true != pair:
                    differ += 1
                    w = wa * wo * wq
                    extra += w * (pair - true)
                    prob += w
                    if worst is None or pair - true > worst[0]:
                        worst = (pair - true, " ".join(tok[a]), " ".join(tok[o]), " ".join(tok[q]), pair, true)
    pairwise = json.loads((ROOT / "results" / "results.json").read_text())["techniques"]["Cancellations"]["saved"]
    out = {"triples_recounted": checked, "triples_differing": differ, "probability_of_such_a_state": prob,
           "mean_extra_moves_saved": extra, "pairwise_mean_saving": pairwise, "corrected_mean_saving": pairwise + extra, "largest_case": worst}
    (ROOT / "results" / "triple_cancellation.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
