"""Exact movecount model of FR-last-slot + last-layer.

Every technique is evaluated on the same 9,331,200 states (placement p, LL state l),
so totals are exact integer counts, not samples, and comparisons are paired.

A solve is a chain of segments:  [AUF LS-alg] [AUF OLL] [AUF PLL AUF]
A segment's cost is its STM length. With cancellations on, an AUF merges into its alg
and the end of each segment merges with the start of the next one.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import data as D
from .cube import ROTATIONS, _frame_lists, auf_tokens, join_savings, merge_tokens, parse, perm_of, stm
from .states import (
    LLSpace, SOLVED, build_placements, ls_corner_pattern, ls_edge_pattern, ls_state_ok, pair_position, parse_labels,
)

K_ORDER = (0, 1, 3, 2)  # prefer no AUF; U2 is last among equal-cost AUFs
MAXT = 90


@dataclass(frozen=True)
class Config:
    name: str = "baseline"
    inserts: tuple = ()              # subset of ("U2", "sledge")
    ls_rule: str = "two_step"        # two_step | three_step | always_U2 | always_sledge
    cancel_auf: bool = False         # AUF merges into the alg next to it
    cancel_alg: bool = False         # consecutive algs merge at their boundary
    oll_angles: str = "off"          # off | speed | min
    pll_angles: str = "off"          # off | speed | min
    exclude_awk: bool = False
    coll: str = "off"                # off | two_step | always | ollcp_two_step | ollcp_always
    ll_mode: str = "cfop"            # cfop | zb (ZBLL on oriented edges) | 1lll
    ls_mode: str = "f2l"             # f2l | zbls | ols
    zb_pick: str = "standard"        # standard (votes, then speed) | moves (shortest listed)


class Segs:
    """Registry of distinct AUF + alg + AUF segments."""

    def __init__(self):
        self.ids = {}
        self.tokens, self.raw, self.internal = [], [], []

    def add(self, k, body, post=0):
        key = (k, body, post)
        i = self.ids.get(key)
        if i is None:
            i = len(self.tokens)
            self.ids[key] = i
            pre, post_t = auf_tokens(k), auf_tokens(post)
            self.tokens.append(pre + body + post_t)
            self.raw.append(len(pre) + stm(body) + len(post_t))
            saved = join_savings(pre, body) if pre else 0
            if post_t:
                saved += join_savings(pre + body, post_t)
            self.internal.append(saved)
        return i


class Model:
    def __init__(self):
        self.ll = LLSpace()
        ll = self.ll
        self.placements = build_placements(ll)
        self.pos_list = list(self.placements.keys())
        self.segs = Segs()
        self.report = {}
        self.ru = [ll.apply(auf_tokens(a)) for a in range(4)]
        self.lu = [ll.left_u(b) for b in range(4)]
        self._load_sets()
        self._build_f2l()
        self._end_sig, self._start_sig = [], []
        self._end_ids, self._start_ids = {}, {}
        self._join_memo = {}
        self._ll_cache = {}
        self._ls_cache = {}

    # ------------------------------------------------ cases on the LL space
    def _assign(self, cases, key, accept, label):
        """Map each accepted LL state to its case by orbit under U^b . l . U^a."""
        ll = self.ll
        case_of = np.full(ll.N, -1, dtype=np.int32)
        rows = np.array([ll.index_of(parse_labels(c.state)[None, :])[0] for c in cases])
        if (rows < 0).any():
            raise ValueError(f"{label}: case state outside the LL space")
        orbits = 0
        for i, r in enumerate(rows):
            orbit = {key[self.ru[a][self.lu[b][r]]] for a in range(4) for b in range(4)}
            mask = np.isin(key, list(orbit)) & accept
            case_of[mask] = i
        missing = int((accept & (case_of < 0)).sum())
        self.report[f"{label}_cases"] = len(cases)
        self.report[f"{label}_uncovered_states"] = missing
        return case_of

    def _load_sets(self):
        ll = self.ll
        full_key = ll._code
        eo_key = (ll.co @ 3 ** np.arange(4)) * 16 + ll.eo @ 2 ** np.arange(4)
        self.oll_cases = D.load_set("oll")
        self.oll_case = self._assign(self.oll_cases, eo_key, ~ll.oriented, "oll")
        self.pll_cases = D.load_set("pll")
        self.pll_case = self._assign(self.pll_cases, full_key, ll.oriented & ~ll.is_u, "pll")
        cp_key = (ll.cp @ 4 ** np.arange(4)) * 216 + eo_key
        coll_all = D.load_set("coll")
        # true COLL: the 329-case file restricted to oriented edges (all corner twists)
        self.coll_cases = [c for c in coll_all if ll.eo0[ll.index_of(parse_labels(c.state)[None, :])[0]]]
        self.coll_case = self._assign(self.coll_cases, cp_key, ll.eo0 & ~ll.co0, "coll")
        self.ollcp_cases = coll_all
        self.ollcp_case = self._assign(coll_all, cp_key, ~(ll.eo0 & ll.co0), "ollcp")
        self.zbll_cases = D.load_set("zbll")
        self.zbll_case = self._assign(self.zbll_cases, full_key, ll.eo0 & ~ll.co0, "zbll")
        self.lll_cases = D.load_set("1lll")
        self.lll_case = self._assign(self.lll_cases, full_key, ~ll.is_u, "1lll")

    # ------------------------------------------------ one alg on a set of states
    def _reach(self, rows, body, goal):
        """For alg `body`, per state: ok[k], after[k] for each pre-AUF k."""
        out_ok, out_after = [], []
        for k in range(4):
            res = self._apply_rows(rows, auf_tokens(k) + body)
            if res is None:
                out_ok.append(np.zeros(len(rows), bool))
                out_after.append(np.zeros(len(rows), np.int32))
            else:
                out_ok.append(goal[res])
                out_after.append(res)
        return np.array(out_ok), np.array(out_after)

    def _apply_rows(self, rows, tokens):
        ll = self.ll
        lab = ll.recenter(ll.labels[rows][:, perm_of(tokens)])
        idx = ll.index_of(lab)
        return None if (idx < 0).any() else idx.astype(np.int32)

    def _best_auf(self, ok, after, finishing):
        """Cheapest valid pre-AUF per state; with a post-AUF cost when `finishing`."""
        ll = self.ll
        n = ok.shape[1]
        best_c = np.full(n, 99)
        best_k = np.full(n, -1, dtype=np.int8)
        for k in K_ORDER:
            c = np.full(n, 99)
            m = ok[k]
            c[m] = (k != 0)
            if finishing:
                c[m] += ll.post_of[after[k][m]] != 0
            better = c < best_c
            best_c[better], best_k[better] = c[better], k
        return best_k, best_c

    def verify_set(self, cases, case_of, goal):
        """Simulate every listed alg on every state of its case."""
        checked = failed = 0
        bad = []
        for i, c in enumerate(cases):
            rows = np.nonzero(case_of == i)[0]
            for a in c.algs:
                ok, _ = self._reach(rows, a.tokens, goal)
                checked += 1
                if not ok.any(axis=0).all():
                    failed += 1
                    bad.append((c.name, a.text))
        return {"algs_checked": checked, "failed": failed, "examples": bad[:5]}

    def ll_table(self, cases, case_of, goal, finishing, angles="off", exclude_awk=False):
        """Selected segment and next state for every state of a set, under an angle rule."""
        ll = self.ll
        seg = np.full(ll.N, -1, dtype=np.int32)
        nxt = np.full(ll.N, -1, dtype=np.int32)
        alg_of = np.full(ll.N, -1, dtype=np.int32)
        for i, c in enumerate(cases):
            rows = np.nonzero(case_of == i)[0]
            if not len(rows):
                continue
            base = D.baseline_of(c.algs)
            okb, afterb = self._reach(rows, base.tokens, goal)
            k_b, cost_b = self._best_auf(okb, afterb, finishing)
            assert (k_b >= 0).all(), (c.name, base.text)
            choice_alg = np.full(len(rows), -1, dtype=np.int32)
            choice_k = k_b.copy()
            choice_after = afterb[k_b, np.arange(len(rows))]
            base_idx = c.algs.index(base)
            choice_alg[:] = base_idx
            if angles != "off":
                alts = D.alternates_of(c.algs, base)
                if exclude_awk:
                    alts = [a for a in alts if not a.awk]
                alts = sorted(alts, key=lambda a: (a.score, a.moves))
                # cost of the current choice, with the post-AUF when this step finishes the cube
                cur = cost_b + base.moves
                for a in alts:
                    ok, after = self._reach(rows, a.tokens, goal)
                    usable = ok[0]  # "gotten from that angle": no AUF
                    post = ll.post_of[after[0]] != 0 if finishing else np.zeros(len(rows), bool)
                    this = a.moves + post
                    if angles == "speed":
                        # at the baseline's own angle keep the baseline; elsewhere the
                        # fastest qualifying alg from this angle replaces AUF + baseline
                        take = usable & (k_b != 0) & (choice_alg == base_idx)
                    else:  # "min": the move-minimising choice among qualifying algs
                        take = usable & (this < cur)
                    choice_alg[take] = c.algs.index(a)
                    choice_k[take] = 0
                    choice_after[take] = after[0][take]
                    cur = np.where(take, this, cur)
            nxt[rows] = choice_after
            alg_of[rows] = choice_alg
            for j in np.unique(choice_alg):
                body = c.algs[j].tokens
                for k in np.unique(choice_k[choice_alg == j]):
                    m = (choice_alg == j) & (choice_k == k)
                    if finishing:
                        posts = ll.post_of[choice_after[m]]
                        for pj in np.unique(posts):
                            mm = np.nonzero(m)[0][posts == pj]
                            seg[rows[mm]] = self.segs.add(int(k), body, int(pj))
                    else:
                        seg[rows[m]] = self.segs.add(int(k), body)
        return seg, nxt, alg_of

    # ------------------------------------------------ F2L
    def _pair_ok(self, tokens):
        r = LLSpace.recenter(SOLVED[perm_of(tokens)][None, :])[0]
        return ls_state_ok(r)

    def _strip_auf(self, tokens):
        """Drop U-layer turns that are really a setup AUF, also behind a y rotation."""
        toks = list(tokens)
        i = 0
        while i < len(toks):
            if toks[i][0] == "U":
                del toks[i]
            elif toks[i][0] == "y":
                i += 1
            else:
                break
        return tuple(toks)

    def _build_f2l(self):
        f2l = D.load_scdb_f2l()
        self.f2l_cases = []
        orbit_of = {}
        rejected = 0
        total = 0
        for idx, c in enumerate(f2l):
            setup = parse(c["setup"])
            pos0 = pair_position(SOLVED[perm_of(setup)])
            valid = []
            for a in c["algs"]:
                total += 1
                toks = parse(a["alg"])
                if any(self._pair_ok(setup + auf_tokens(k) + toks) for k in range(4)):
                    body = self._strip_auf(toks)
                    if any(self._pair_ok(setup + auf_tokens(k) + body) for k in range(4)):
                        valid.append((a["votes"], a["alg"], body))
                        continue
                rejected += 1
            top = max(valid, key=lambda v: v[0])
            self.f2l_cases.append({"name": c["name"], "group": c.get("subgroup", ""), "text": top[1], "votes": top[0], "body": top[2], "moves": stm(top[2]), "pos": pos0, "n_valid": len(valid), "valid": valid})
            for a in range(4):
                orbit_of[pair_position(SOLVED[perm_of(setup + auf_tokens(a))])] = idx
        self.report["f2l_algs_listed"] = total
        self.report["f2l_algs_rejected_other_slot"] = rejected
        self.orbit_of = orbit_of
        solved = pair_position(SOLVED)
        self.f2l_by_pos = {}
        for p in self.pos_list:
            if p == solved:
                self.f2l_by_pos[p] = None
            else:
                self.f2l_by_pos[p] = orbit_of[p]
        self.report["f2l_cases"] = len(f2l)
        self.report["placements_covered"] = sum(v is not None for v in self.f2l_by_pos.values()) + 1

    def _auf_for(self, p, body):
        r = self.placements[p]
        for k in K_ORDER:
            if self._pair_ok(r + auf_tokens(k) + body):
                return k
        return None

    def ls_variants(self, body):
        """The insert family from the author's definitions, with tails swapped."""
        forms = {"std": parse("R U' R'"), "U2": parse("U R U2 R'"), "sledge": parse("U' R' F R F'")}
        out = {"base": body}
        for name, f in forms.items():
            if len(body) >= len(f) and body[-len(f):] == f:
                prefix = body[: -len(f)]
                for other, g in forms.items():
                    if other != name:
                        out[other] = merge_tokens(prefix, g)
                return out, name
        return out, None

    def ls_options(self, cfg):
        """Per placement: list of (label, k, body) options; first is the baseline alg."""
        key = cfg.inserts
        if key in self._ls_cache:
            return self._ls_cache[key]
        out = {}
        for p in self.pos_list:
            idx = self.f2l_by_pos[p]
            if idx is None:
                out[p] = [("skip", 0, ())]
                continue
            body = self.f2l_cases[idx]["body"]
            opts = [("base", self._auf_for(p, body), body)]
            vs, _ = self.ls_variants(body)
            allowed = []
            if "U2" in cfg.inserts:
                allowed += ["U2", "std"]
            if "sledge" in cfg.inserts:
                allowed += ["sledge", "std"]
            for name in allowed:
                if name in vs:
                    b = self._strip_auf(vs[name])
                    if all(b != o[2] for o in opts):
                        k = self._auf_for(p, b)
                        if k is not None:
                            opts.append((name, k, b))
            out[p] = opts
        self._ls_cache[key] = out
        return out

    # ------------------------------------------------ joins
    def _sync_sigs(self):
        """Boundary signatures (rotation-free ends and starts) of any new segments."""
        toks = self.segs.tokens
        for t in toks[len(self._end_sig):]:
            st, en = _frame_lists(t)
            self._end_sig.append(self._end_ids.setdefault(en[-16:], len(self._end_ids)))
            self._start_sig.append(self._start_ids.setdefault(st[:16], len(self._start_ids)))
        self._end_arr = np.array(self._end_sig + [-1], dtype=np.int64)
        self._start_arr = np.array(self._start_sig, dtype=np.int64)
        self._end_rev = {v: k for k, v in self._end_ids.items()}
        self._start_rev = {v: k for k, v in self._start_ids.items()}

    def joins(self, end_segs, start_segs):
        """Moves saved at each boundary; an end of len(tokens) means nothing came before."""
        end_segs = np.where(np.asarray(end_segs) >= len(self.segs.tokens), -1, end_segs)
        e = self._end_arr[end_segs]
        e = np.where(np.asarray(end_segs) < 0, -1, e)
        s = self._start_arr[start_segs]
        e, s = np.broadcast_arrays(e, s)
        pairs = (e + 1) * 10_000_000 + s
        uniq, inv = np.unique(pairs, return_inverse=True)
        vals = np.empty(len(uniq), dtype=np.int32)
        for i, q in enumerate(uniq):
            q = int(q)
            v = self._join_memo.get(q)
            if v is None:
                ei, si = q // 10_000_000 - 1, q % 10_000_000
                v = 0 if ei < 0 else join_savings(self._end_rev[ei], self._start_rev[si])
                self._join_memo[q] = v
            vals[i] = v
        return vals[inv].reshape(e.shape)

    def cost_of(self, cfg, seg):
        raw = np.asarray(self.segs.raw)
        if cfg.cancel_auf:
            raw = raw - np.asarray(self.segs.internal)
        return raw[seg]

    # ------------------------------------------------ LL level
    def ll_level(self, cfg):
        key = (cfg.oll_angles, cfg.pll_angles, cfg.exclude_awk, cfg.coll, cfg.ll_mode)
        if key in self._ll_cache:
            return self._ll_cache[key]
        ll = self.ll
        # PLL for oriented states; a pure AUF state is a skip
        pseg, pnext, palg = self.ll_table(self.pll_cases, self.pll_case, ll.is_u, True, cfg.pll_angles, cfg.exclude_awk)
        for u in ll.u_power:
            pseg[u] = self.segs.add(0, (), int(ll.post_of[u]))
            pnext[u] = u
        oseg, onext, oalg = self.ll_table(self.oll_cases, self.oll_case, ll.oriented, False, cfg.oll_angles, cfg.exclude_awk)
        options = [(oseg, onext)]
        if cfg.oll_angles == "two_step":
            oseg, onext, _ = self.ll_table(self.oll_cases, self.oll_case, ll.oriented, False)
            options = [(oseg, onext)] + self._alt_option_tables(cfg)
        forced = None
        if cfg.ll_mode == "zb":
            zseg, znext, _ = self.ll_table(self.zbll_cases, self.zbll_case, ll.is_u, True)
            forced = (zseg, znext)
        elif cfg.ll_mode == "1lll":
            zseg, znext, _ = self.ll_table(self.lll_cases, self.lll_case, ll.is_u, True)
            forced = (zseg, znext)
        elif cfg.coll != "off":
            cases, case_of = (self.coll_cases, self.coll_case) if not cfg.coll.startswith("ollcp") else (self.ollcp_cases, self.ollcp_case)
            goal = ll.corners_u & ll.eo0
            cseg, cnext, _ = self.ll_table(cases, case_of, goal, False)
            options.append((cseg, cnext))
        res = {"pll_seg": pseg, "options": options, "forced": forced}
        self._ll_cache[key] = res
        return res

    def _alt_option_tables(self, cfg):
        """The i-th qualifying OLL alternate of each case, usable with no AUF only."""
        ll = self.ll
        tables = []
        for i, c in enumerate(self.oll_cases):
            rows = np.nonzero(self.oll_case == i)[0]
            base = D.baseline_of(c.algs)
            alts = sorted(D.alternates_of(c.algs, base), key=lambda a: (a.score, a.moves))
            for j, a in enumerate(alts):
                while len(tables) <= j:
                    tables.append((np.full(ll.N, -1, np.int32), np.full(ll.N, -1, np.int32)))
                ok, after = self._reach(rows, a.tokens, ll.oriented)
                m = ok[0]
                tables[j][0][rows[m]] = self.segs.add(0, a.tokens)
                tables[j][1][rows[m]] = after[0][m]
        return tables

    def ll_cost(self, cfg, e_seg):
        """(oll-stage cost, pll-stage cost, used option) per LL state, after segment e_seg."""
        L = self.ll_level(cfg)
        ll = self.ll
        pseg = L["pll_seg"]
        n_none = len(self.segs.tokens)
        e = n_none if e_seg is None else e_seg

        def pll_after(end_seg, x):
            s = pseg[x]
            out = self.cost_of(cfg, s).astype(np.int32)
            if cfg.cancel_alg:
                out = out - self.joins(end_seg, s)
            return out

        if L["forced"] is not None:
            zseg, znext = L["forced"]
            mask = zseg >= 0
            first = np.zeros(ll.N, np.int32)
            zc = self.cost_of(cfg, zseg[mask]).astype(np.int32)
            if cfg.cancel_alg:
                zc = zc - self.joins(np.full(mask.sum(), e), zseg[mask])
            first[mask] = zc
            # no alg needed at the "orientation" stage for a 1-look case; a U-only state uses PLL
        results = []
        for oseg, onext in L["options"]:
            usable = oseg >= 0
            s = np.where(usable, oseg, 0)
            c1 = np.where(usable, self.cost_of(cfg, s), 0).astype(np.int32)
            ends = np.where(usable, s, e)
            if cfg.cancel_alg:
                c1 = c1 - np.where(usable, self.joins(np.full(ll.N, e), s), 0)
            nx = np.where(usable, onext, np.arange(ll.N))
            # skip of the orientation step is only valid when already oriented
            valid = usable | ll.oriented
            c2 = np.zeros(ll.N, np.int32)
            ok = ll.oriented[nx] & valid
            c2[ok] = pll_after(ends[ok], nx[ok])
            results.append((c1, c2, valid & ok))
        c1, c2, valid = results[0]
        used = np.zeros(ll.N, np.int8)
        assert valid.all()
        for i, (d1, d2, v) in enumerate(results[1:], start=1):
            has_alg = L["options"][i][0] >= 0
            if cfg.coll in ("always", "ollcp_always"):
                take = v & has_alg
            else:  # two-step: this alg plus the rest of the solve, strictly better only
                take = v & has_alg & (d1 + d2 < c1 + c2)
            c1, c2 = np.where(take, d1, c1), np.where(take, d2, c2)
            used[take] = 1
        if L["forced"] is not None:
            zseg, _ = L["forced"]
            mask = zseg >= 0
            c1 = np.where(mask, first, c1)
            c2 = np.where(mask, 0, c2)
            used[mask] = 2
        return c1, c2, used

    # ------------------------------------------------ whole solve
    def ls_tables(self, cfg):
        """Per placement: option arrays (label, seg per l or scalar, next per l)."""
        ll = self.ll
        if cfg.ls_mode in ("zbls", "ols"):
            return self._zb_tables(cfg)
        opts = self.ls_options(cfg)
        out = {}
        for p, lst in opts.items():
            r = self.placements[p]
            row = []
            for label, k, body in lst:
                if label == "skip":
                    row.append((label, None, np.arange(ll.N, dtype=np.int32)))
                    continue
                seg = self.segs.add(k, body)
                nxt = ll.apply(r + auf_tokens(k) + body)
                assert nxt is not None, (p, body)
                row.append((label, seg, nxt))
            out[p] = row
        return out

    def _zb_tables(self, cfg):
        ll = self.ll
        cases = D.load_set("zbls") if cfg.ls_mode == "zbls" else D.load_set("ols") + D.load_set("ols-more")
        with_twist = cfg.ls_mode == "ols"
        goal = ll.oriented if with_twist else ll.eo0

        def pat(lab):
            e = ls_edge_pattern(lab).astype(np.int64)
            return e * 243 + ls_corner_pattern(lab) if with_twist else e

        table = {}
        for ci, c in enumerate(cases):
            base = parse_labels(c.state)
            for a in range(4):
                lab = base[perm_of(auf_tokens(a))][None, :]
                table.setdefault((pair_position(lab[0]), int(pat(lab)[0])), ci)
        self.report[f"{cfg.ls_mode}_cases"] = len(cases)
        out = {}
        uncovered = 0
        rejected = 0
        for p in self.pos_list:
            r = self.placements[p]
            if self.f2l_by_pos[p] is None:
                out[p] = [("skip", None, np.arange(ll.N, dtype=np.int32))]
                continue
            lab = ll.labels[:, perm_of(r)]
            pats = pat(lab)
            seg = np.full(ll.N, -1, dtype=np.int32)
            nxt = np.full(ll.N, -1, dtype=np.int32)
            for v in np.unique(pats):
                rows = np.nonzero(pats == v)[0]
                ci = table.get((p, int(v)))
                if ci is None:
                    uncovered += len(rows)
                    continue
                order = sorted(cases[ci].algs, key=D.baseline_key if cfg.zb_pick == "standard" else (lambda a: (a.moves, a.score)))
                done = False
                for a in order:
                    body = self._strip_auf(a.tokens)
                    for k in K_ORDER:
                        res = self._apply_rows(rows, r + auf_tokens(k) + body)
                        if res is not None and goal[res].all():
                            seg[rows] = self.segs.add(k, body)
                            nxt[rows] = res
                            done = True
                            break
                    if done:
                        break
                    rejected += 1
                assert done, (c.name)
            out[p] = [("zb", seg, nxt)]
        self.report[f"{cfg.ls_mode}_uncovered_states"] = uncovered
        self.report[f"{cfg.ls_mode}_algs_rejected_while_selecting"] = rejected
        return out

    def evaluate(self, cfg: Config, baseline=None):
        ll = self.ll
        L = self.ll_level(cfg)
        tables = self.ls_tables(cfg)
        if cfg.cancel_alg:
            self._sync_sigs()
        H_total, H_ls, H_lsoll, H_ll, H_oll, H_pll = (np.zeros(MAXT, np.int64) for _ in range(6))
        H_diff = np.zeros(2 * MAXT + 1, np.int64)
        usage = {}
        ll_use = np.zeros(3, np.int64)
        by_p = {}
        cache = {}
        for p, opts in tables.items():
            per = []
            for label, seg, nxt in opts:
                if isinstance(seg, np.ndarray):
                    # state-dependent LS segment (ZBLS, OLS)
                    ls_c = self.cost_of(cfg, seg).astype(np.int32)
                    a1, a2, u = np.zeros(ll.N, np.int32), np.zeros(ll.N, np.int32), np.zeros(ll.N, np.int8)
                    for s in np.unique(seg):
                        m = seg == s
                        c1, c2, used = self._ll_cost_cached(cfg, int(s), cache)
                        a1[m], a2[m], u[m] = c1[nxt[m]], c2[nxt[m]], used[nxt[m]]
                else:
                    ls_c = np.full(ll.N, 0 if seg is None else int(self.cost_of(cfg, seg)), np.int32)
                    c1, c2, used = self._ll_cost_cached(cfg, seg, cache)
                    a1, a2, u = c1[nxt], c2[nxt], used[nxt]
                per.append((label, ls_c, a1, a2, u))
            if len(per) == 1:
                ls_c, a1, a2, u = per[0][1:]
                pick = np.zeros(ll.N, np.int8)
            else:
                rule = cfg.ls_rule
                if rule == "two_step":
                    score = np.array([q[1] + q[2] for q in per])
                elif rule == "three_step":
                    score = np.array([q[1] + q[2] + q[3] for q in per])
                elif rule.startswith("always_"):
                    want = rule.split("_", 1)[1]
                    idx = next((i for i, q in enumerate(per) if q[0] == want), 0)
                    score = np.array([np.full(ll.N, 0 if i == idx else 1) for i in range(len(per))])
                else:
                    raise ValueError(rule)
                pick = np.argmin(score, axis=0).astype(np.int8)  # first (baseline) wins ties
                rows = np.arange(ll.N)
                ls_c = np.array([q[1] for q in per])[pick, rows]
                a1 = np.array([q[2] for q in per])[pick, rows]
                a2 = np.array([q[3] for q in per])[pick, rows]
                u = np.array([q[4] for q in per])[pick, rows]
            for i, q in enumerate(per):
                usage[q[0]] = usage.get(q[0], 0) + int((pick == i).sum())
            ll_use += np.bincount(u, minlength=3)
            total = ls_c + a1 + a2
            # A step cancelled away entirely can save more than its own length; charge
            # the excess to the step before it, so no stage is negative.
            a1 = a1 + np.minimum(a2, 0)
            a2 = np.maximum(a2, 0)
            ls_c = ls_c + np.minimum(a1, 0)
            a1 = np.maximum(a1, 0)
            H_total += np.bincount(total, minlength=MAXT)[:MAXT]
            H_ls += np.bincount(ls_c, minlength=MAXT)[:MAXT]
            H_lsoll += np.bincount(ls_c + a1, minlength=MAXT)[:MAXT]
            H_ll += np.bincount(a1 + a2, minlength=MAXT)[:MAXT]
            H_oll += np.bincount(a1, minlength=MAXT)[:MAXT]
            H_pll += np.bincount(a2, minlength=MAXT)[:MAXT]
            by_p[p] = total.astype(np.int16)
            if baseline is not None:
                H_diff += np.bincount(baseline["by_p"][p] - total + MAXT, minlength=2 * MAXT + 1)
        assert H_total.sum() == 150 * ll.N
        return {"cfg": cfg, "total": H_total, "ls": H_ls, "lsoll": H_lsoll, "ll": H_ll, "oll": H_oll, "pll": H_pll, "diff": H_diff if baseline is not None else None,
                "usage": usage, "ll_use": ll_use, "by_p": by_p}

    def _ll_cost_cached(self, cfg, seg, cache):
        if not cfg.cancel_alg:
            seg = None
        if seg not in cache:
            cache[seg] = self.ll_cost(cfg, seg)
        return cache[seg]
