"""The exact state space: all last-layer states and all last-slot placements.

LL   : the 62,208 states with F2L solved (4!*3^3*4!*2^3/2, AUF offsets included).
LSLL : LL state l followed by a placement move r_p: state = l . r_p, 150 placements
       of the FR pair, so 150 * 62,208 = 9,331,200 states.
"""
from __future__ import annotations

from itertools import permutations

import numpy as np

from .cube import IDENTITY, PERM, auf_tokens, compose, parse, perm_of

CORNER_F = [[8, 9, 20], [6, 18, 38], [0, 36, 47], [2, 45, 11], [29, 26, 15], [27, 44, 24], [33, 53, 42], [35, 17, 51]]
EDGE_F = [[5, 10], [7, 19], [3, 37], [1, 46], [32, 16], [28, 25], [30, 43], [34, 52], [23, 12], [21, 41], [50, 39], [48, 14]]
CORNER_N = ["URF", "UFL", "ULB", "UBR", "DFR", "DLF", "DBL", "DRB"]
EDGE_N = ["UR", "UF", "UL", "UB", "DR", "DF", "DL", "DB", "FR", "FL", "BL", "BR"]
FACE = "URFDLB"
CENTERS = np.array([4, 13, 22, 31, 40, 49])
SOLVED = np.repeat(np.arange(6, dtype=np.uint8), 9)
LL_STICKERS = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 18, 19, 20, 36, 37, 38, 45, 46, 47])
FIXED_STICKERS = np.array([i for i in range(54) if i not in set(LL_STICKERS.tolist())])
_POW6 = 6 ** np.arange(len(LL_STICKERS), dtype=np.int64)


def parse_labels(state: str) -> np.ndarray:
    return np.array([FACE.index(c) for c in state], dtype=np.uint8)


def _parity(p):
    return sum(p[j] < p[i] for i in range(len(p)) for j in range(i + 1, len(p))) % 2


def _labels(cp, co, ep, eo):
    out = SOLVED.copy()
    for i in range(4):
        col = [FACE.index(c) for c in CORNER_N[cp[i]]]
        for k in range(3):
            out[CORNER_F[i][(co[i] + k) % 3]] = col[k]
        col = [FACE.index(c) for c in EDGE_N[ep[i]]]
        for k in range(2):
            out[EDGE_F[i][(eo[i] + k) % 2]] = col[k]
    return out


class LLSpace:
    """All 62,208 LL states with fast right-multiplication by any move sequence."""

    N = 62208

    def __init__(self):
        cps, cos, eps, eos, labels = [], [], [], [], []
        perms = list(permutations(range(4)))
        twists = []
        for t in range(27):
            a = [t % 3, (t // 3) % 3, (t // 9) % 3]
            twists.append(a + [(6 - sum(a)) % 3])
        flips = []
        for t in range(8):
            a = [t & 1, (t >> 1) & 1, (t >> 2) & 1]
            flips.append(a + [sum(a) % 2])
        for cp in perms:
            for co in twists:
                for ep in perms:
                    if _parity(cp) != _parity(ep):
                        continue
                    for eo in flips:
                        cps.append(cp); cos.append(co); eps.append(ep); eos.append(eo)
                        labels.append(_labels(cp, co, ep, eo))
        self.cp, self.co = np.array(cps, np.int8), np.array(cos, np.int8)
        self.ep, self.eo = np.array(eps, np.int8), np.array(eos, np.int8)
        self.labels = np.array(labels, dtype=np.uint8)
        assert len(self.labels) == self.N
        keys = self.labels[:, LL_STICKERS].astype(np.int64) @ _POW6
        self._order = np.argsort(keys)
        self._sorted = keys[self._order]
        assert len(np.unique(keys)) == self.N
        self._cache = {}

        # U powers and goal predicates
        self.solved = self.index_of(SOLVED[None, :])[0]
        self.u_power = np.array([self.index_of(SOLVED[perm_of(auf_tokens(k))][None, :])[0] for k in range(4)])
        self.oriented = (self.co == 0).all(1) & (self.eo == 0).all(1)
        self.eo0 = (self.eo == 0).all(1)
        self.co0 = (self.co == 0).all(1)
        self.is_u = np.isin(np.arange(self.N), self.u_power)
        # j with state == U^(-j) needed... `post_of[l]` = k with l . U^k solved, for U-layer-only states
        self.post_of = np.full(self.N, -1, dtype=np.int8)
        for k in range(4):
            inv = self.u_power[(4 - k) % 4]
            self.post_of[inv] = k
        # corners solved up to U (with whatever edges)
        rel = (self.cp - np.arange(4)) % 4
        self.corners_u = self.co0 & (rel == rel[:, :1]).all(1)
        self.ori_pattern = (self.co @ 3 ** np.arange(4)) * 16 + self.eo @ 2 ** np.arange(4)

    def index_of(self, labels: np.ndarray):
        """Index of each row in the LL space, or -1 when F2L is not intact."""
        keys = labels[:, LL_STICKERS].astype(np.int64) @ _POW6
        pos = np.searchsorted(self._sorted, keys)
        pos[pos >= self.N] = self.N - 1
        hit = self._sorted[pos] == keys
        hit &= (labels[:, FIXED_STICKERS] == SOLVED[FIXED_STICKERS]).all(1)
        out = np.where(hit, self._order[pos], -1)
        return out

    @staticmethod
    def recenter(labels: np.ndarray):
        """Relabel stickers by the current centers, removing any net rotation."""
        centers = labels[0, CENTERS]
        if (centers == np.arange(6)).all():
            return labels
        mp = np.zeros(6, dtype=np.uint8)
        mp[centers] = np.arange(6)
        return mp[labels]

    def apply(self, tokens: tuple[str, ...], rows=None):
        """Index map l -> l . tokens, or None unless every result is an LL state."""
        if rows is None:
            key = tokens
            if key in self._cache:
                return self._cache[key]
        labels = self.labels if rows is None else self.labels[rows]
        res = self.recenter(labels[:, perm_of(tokens)])
        idx = self.index_of(res)
        out = None if (idx < 0).any() else idx.astype(np.int32)
        if rows is None:
            self._cache[key] = out
        return out

    def left_u(self, b: int) -> np.ndarray:
        """Index map l -> U^b . l (a change of colour labels)."""
        u = self.labels[self.u_power[b % 4]]
        # applying the state l's permutation on top of U^b: compute by cubie composition
        sh = {0: [0, 1, 2, 3], 1: [3, 0, 1, 2], 2: [2, 3, 0, 1], 3: [1, 2, 3, 0]}[b % 4]
        # U^b as a cubie state, from its labels: position i holds piece ucp[i]
        ucp, uep = self.cp[self.u_power[b % 4]], self.ep[self.u_power[b % 4]]
        cp, ep = ucp[self.cp], uep[self.ep]
        lab = np.array(self.labels)
        out = np.empty(self.N, dtype=np.int32)
        code = self._cubie_key(cp, self.co, ep, self.eo)
        order = np.argsort(self._code)
        out = order[np.searchsorted(self._code[order], code)]
        return out

    @property
    def _code(self):
        if not hasattr(self, "_code_cache"):
            self._code_cache = self._cubie_key(self.cp, self.co, self.ep, self.eo)
        return self._code_cache

    @staticmethod
    def _cubie_key(cp, co, ep, eo):
        return (cp @ 4 ** np.arange(4)) * 4**4 * 3**4 * 2**4 + (ep @ 4 ** np.arange(4)) * 3**4 * 2**4 + (co @ 3 ** np.arange(4)) * 2**4 + eo @ 2 ** np.arange(4)


# ------------------------------------------------------------------ last slot

def pair_position(labels_row):
    """(corner position, twist, edge position, flip) of the FR pair, on a labelled state."""
    pos = None
    for i, f in enumerate(CORNER_F):
        s = labels_row[f]
        if set(s.tolist()) == {3, 1, 2}:  # D, R, F stickers
            pos = (i, int(np.argmax(s == 3)))
    assert pos is not None
    for i, f in enumerate(EDGE_F):
        s = labels_row[f]
        if set(s.tolist()) == {2, 1}:
            return (pos[0], pos[1], i, 0 if s[0] == 2 else 1)
    raise AssertionError


def build_placements(ll: LLSpace):
    """BFS over pair-moving triggers to a move sequence reaching each of the 150 placements."""
    gens = [parse(a) for a in ["U", "U'", "U2", "R U R'", "R U' R'", "R U2 R'", "F' U F", "F' U' F", "F' U2 F", "R' F R F'", "F R' F' R"]]
    start = pair_position(SOLVED)
    found = {start: ()}
    frontier = [()]
    while frontier and len(found) < 150:
        nxt = []
        for seq in frontier:
            for g in gens:
                s = seq + g
                pos = pair_position(SOLVED[perm_of(s)])
                if pos not in found:
                    found[pos] = s
                    nxt.append(s)
        frontier = nxt
    assert len(found) == 150, len(found)
    return found


def ls_state_ok(labels_row) -> bool:
    """Pair FR solved and the other three pairs and cross intact, after recentring."""
    return (labels_row[FIXED_STICKERS] == SOLVED[FIXED_STICKERS]).all()


EDGE_LS_POS = [0, 1, 2, 3, 8]


def ls_edge_pattern(labels: np.ndarray) -> np.ndarray:
    """5-bit orientation of the edges at UR, UF, UL, UB, FR positions on labelled rows."""
    out = np.zeros(len(labels), dtype=np.int64)
    for b, pos in enumerate(EDGE_LS_POS):
        a0, a1 = EDGE_F[pos]
        s0, s1 = labels[:, a0], labels[:, a1]
        # eo=0 if the sticker pair reads as a named edge in its name's order:
        # U-edges read (U, x); the FR edge reads (F, R).
        ok = (s0 == 0) | ((s0 == 2) & (s1 == 1))
        out |= (~ok).astype(np.int64) << b
    return out


CORNER_LS_POS = [0, 1, 2, 3, 4]


def ls_corner_pattern(labels: np.ndarray) -> np.ndarray:
    """Base-3 twist of the corners at the 4 U positions and DFR, on labelled rows."""
    out = np.zeros(len(labels), dtype=np.int64)
    for b, pos in enumerate(CORNER_LS_POS):
        f = CORNER_F[pos]
        st = labels[:, f]
        tw = np.argmax((st == 0) | (st == 3), axis=1)
        out += tw * 3**b
    return out
