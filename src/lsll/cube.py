"""54-facelet cube engine, the STM metric, and move-sequence simplification.

Facelets are numbered U 0-8, R 9-17, F 18-26, D 27-35, L 36-44, B 45-53.
A move is a permutation p with after[i] = before[p[i]].

STM (slice turn metric): every face, wide or slice turn counts as one move,
whatever its angle (so U2 = 1). Whole-cube rotations x, y, z count as zero.
"""
from __future__ import annotations

import re
from functools import lru_cache

import numpy as np

_BASE = {
    "U": [6, 3, 0, 7, 4, 1, 8, 5, 2, 45, 46, 47, 12, 13, 14, 15, 16, 17, 9, 10, 11, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 18, 19, 20, 39, 40, 41, 42, 43, 44, 36, 37, 38, 48, 49, 50, 51, 52, 53],
    "D": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 24, 25, 26, 18, 19, 20, 21, 22, 23, 42, 43, 44, 33, 30, 27, 34, 31, 28, 35, 32, 29, 36, 37, 38, 39, 40, 41, 51, 52, 53, 45, 46, 47, 48, 49, 50, 15, 16, 17],
    "R": [0, 1, 20, 3, 4, 23, 6, 7, 26, 15, 12, 9, 16, 13, 10, 17, 14, 11, 18, 19, 29, 21, 22, 32, 24, 25, 35, 27, 28, 51, 30, 31, 48, 33, 34, 45, 36, 37, 38, 39, 40, 41, 42, 43, 44, 8, 46, 47, 5, 49, 50, 2, 52, 53],
    "L": [53, 1, 2, 50, 4, 5, 47, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 19, 20, 3, 22, 23, 6, 25, 26, 18, 28, 29, 21, 31, 32, 24, 34, 35, 42, 39, 36, 43, 40, 37, 44, 41, 38, 45, 46, 33, 48, 49, 30, 51, 52, 27],
    "F": [0, 1, 2, 3, 4, 5, 44, 41, 38, 6, 10, 11, 7, 13, 14, 8, 16, 17, 24, 21, 18, 25, 22, 19, 26, 23, 20, 15, 12, 9, 30, 31, 32, 33, 34, 35, 36, 37, 27, 39, 40, 28, 42, 43, 29, 45, 46, 47, 48, 49, 50, 51, 52, 53],
    "B": [11, 14, 17, 3, 4, 5, 6, 7, 8, 9, 10, 35, 12, 13, 34, 15, 16, 33, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 36, 39, 42, 2, 37, 38, 1, 40, 41, 0, 43, 44, 51, 48, 45, 52, 49, 46, 53, 50, 47],
    "x": [18, 19, 20, 21, 22, 23, 24, 25, 26, 15, 12, 9, 16, 13, 10, 17, 14, 11, 27, 28, 29, 30, 31, 32, 33, 34, 35, 53, 52, 51, 50, 49, 48, 47, 46, 45, 38, 41, 44, 37, 40, 43, 36, 39, 42, 8, 7, 6, 5, 4, 3, 2, 1, 0],
    "y": [6, 3, 0, 7, 4, 1, 8, 5, 2, 45, 46, 47, 48, 49, 50, 51, 52, 53, 9, 10, 11, 12, 13, 14, 15, 16, 17, 29, 32, 35, 28, 31, 34, 27, 30, 33, 18, 19, 20, 21, 22, 23, 24, 25, 26, 36, 37, 38, 39, 40, 41, 42, 43, 44],
    "z": [42, 39, 36, 43, 40, 37, 44, 41, 38, 6, 3, 0, 7, 4, 1, 8, 5, 2, 24, 21, 18, 25, 22, 19, 26, 23, 20, 15, 12, 9, 16, 13, 10, 17, 14, 11, 33, 30, 27, 34, 31, 28, 35, 32, 29, 47, 50, 53, 46, 49, 52, 45, 48, 51],
}


def compose(first, second):
    """Permutation for: apply `first`, then `second`."""
    return np.asarray(first)[np.asarray(second)]


def _invert(p):
    q = np.empty(len(p), dtype=np.int16)
    q[np.asarray(p)] = np.arange(len(p), dtype=np.int16)
    return q


PERM: dict[str, np.ndarray] = {}
for _n, _p in _BASE.items():
    _a = np.array(_p, dtype=np.int16)
    PERM[_n], PERM[_n + "'"], PERM[_n + "2"] = _a, _invert(_a), compose(_a, _a)

_DERIVED = {
    "r": ["L", "x"], "l": ["R", "x'"], "u": ["D", "y"], "d": ["U", "y'"], "f": ["B", "z"], "b": ["F", "z'"],
    "M": ["R", "L'", "x'"], "E": ["U", "D'", "y'"], "S": ["F'", "B", "z"],
}
IDENTITY = np.arange(54, dtype=np.int16)
for _n, _parts in _DERIVED.items():
    _a = IDENTITY
    for _t in _parts:
        _a = compose(_a, PERM[_t])
    PERM[_n], PERM[_n + "'"], PERM[_n + "2"] = _a, _invert(_a), compose(_a, _a)

ROTATIONS = {"x", "y", "z"}
_TOKEN = re.compile(r"^[URFDLBurfdlbMESxyz]['2]?$")


def parse(alg: str) -> tuple[str, ...]:
    """Normalise algorithm text to tokens: brackets dropped, Rw->r, R2'->R2, R3->R'."""
    out = []
    for t in re.sub(r"[()\[\],]", " ", alg).split():
        m = re.match(r"^([URFDLB])w(.*)$", t)
        if m:
            t = m.group(1).lower() + m.group(2)
        m = re.match(r"^([A-Za-z])3('?)$", t)
        if m:
            t = m.group(1) + ("" if m.group(2) else "'")
        if re.match(r"^[A-Za-z]2'$", t):
            t = t[:2]
        if not _TOKEN.match(t):
            raise ValueError(f"unrecognised move token {t!r} in {alg!r}")
        out.append(t)
    return tuple(out)


def stm(tokens) -> int:
    """STM length: number of tokens that are not whole-cube rotations."""
    if isinstance(tokens, str):
        tokens = parse(tokens)
    return sum(t[0] not in ROTATIONS for t in tokens)


@lru_cache(maxsize=None)
def perm_of(tokens: tuple[str, ...]) -> np.ndarray:
    p = IDENTITY
    for t in tokens:
        p = compose(p, PERM[t])
    p = np.array(p, dtype=np.int16)
    p.setflags(write=False)
    return p


def inverse(tokens) -> tuple[str, ...]:
    return tuple(t if t.endswith("2") else (t[:-1] if t.endswith("'") else t + "'") for t in reversed(tokens))


AUF = ("", "U", "U2", "U'")


def auf_tokens(k: int) -> tuple[str, ...]:
    return (AUF[k % 4],) if k % 4 else ()


# ---------------------------------------------------------------- simplification
# Every turn acts on an axis, as amounts (mod 4) on its three layers.
#   x axis, direction of R: R=(1,0,0) on (R, M, L); M=(0,-1,0); L=(0,0,-1);
#   r=(1,1,0); l=(0,-1,-1); x=(1,1,1). Similarly y (U,E,D) and z (F,S,B).
_AXIS = {
    "R": ("x", (1, 0, 0)), "M": ("x", (0, 3, 0)), "L": ("x", (0, 0, 3)), "r": ("x", (1, 1, 0)), "l": ("x", (0, 3, 3)),
    "U": ("y", (1, 0, 0)), "E": ("y", (0, 3, 0)), "D": ("y", (0, 0, 3)), "u": ("y", (1, 1, 0)), "d": ("y", (0, 3, 3)),
    "F": ("z", (1, 0, 0)), "S": ("z", (0, 1, 0)), "B": ("z", (0, 0, 3)), "f": ("z", (1, 1, 0)), "b": ("z", (0, 3, 3)),
}
_POWER = {"": 1, "2": 2, "'": 3}


def layer_vector(token: str):
    axis, v = _AXIS[token[0]]
    k = _POWER[token[1:]]
    return axis, tuple((a * k) % 4 for a in v)


def _min_cost_table():
    """Fewest STM turns reaching each layer vector on one axis (BFS)."""
    gens = []
    for base in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)]:
        for k in (1, 2, 3):
            gens.append(tuple((a * k) % 4 for a in base))
    table = {(0, 0, 0): 0}
    frontier = [(0, 0, 0)]
    while frontier:
        nxt = []
        for v in frontier:
            for g in gens:
                w = tuple((a + b) % 4 for a, b in zip(v, g))
                if w not in table:
                    table[w] = table[v] + 1
                    nxt.append(w)
        frontier = nxt
    # An exact whole-cube rotation is free.
    for k in (1, 2, 3):
        table[(k, k, k)] = 0
    return table


_MIN_COST = _min_cost_table()


def _frame_lists(tokens):
    """Rewrite `tokens` without rotations, in the start frame and in the end frame."""
    if not any(t[0] in ROTATIONS for t in tokens):
        return tuple(tokens), tuple(tokens)
    rot = IDENTITY
    start = []
    moves = []
    for t in tokens:
        if t[0] in ROTATIONS:
            rot = compose(rot, PERM[t])
            continue
        start.append(_conjugate_to_start(rot, t))
        moves.append((t, rot))
    final = rot
    end = []
    for t, r in moves:
        # (rotation r) then t, followed later by the remaining rotation:
        # express t in the frame after all rotations.
        rel = compose(_invert(r), final)
        end.append(_conjugate_to_start(_invert(rel), t))
    return tuple(start), tuple(end)


_BY_PERM = {PERM[t].tobytes(): t for t in PERM if t[0] not in ROTATIONS}


def _conjugate_to_start(rot, token):
    """The token T with: `rot` then `token` == T then `rot`."""
    p = compose(compose(rot, PERM[token]), _invert(rot))
    return _BY_PERM[np.array(p, dtype=np.int16).tobytes()]


def _blocks(tokens):
    """Maximal runs on one axis: list of [axis, vector, token_count]."""
    out = []
    for t in tokens:
        axis, v = layer_vector(t)
        if out and out[-1][0] == axis:
            out[-1][1] = tuple((a + b) % 4 for a, b in zip(out[-1][1], v))
            out[-1][2] += 1
        else:
            out.append([axis, v, 1])
    return out


@lru_cache(maxsize=None)
def join_savings(first: tuple[str, ...], second: tuple[str, ...]) -> int:
    """Moves saved by performing `second` straight after `first`.

    Only the boundary may merge: blocks inside either sequence stay as written
    until the blocks in front of them have cancelled away completely, in which
    case the merge cascades (R U + U' R' -> nothing saves 4).
    """
    tail = [list(b) for b in _blocks(_frame_lists(first)[1])]
    head = [list(b) for b in _blocks(_frame_lists(second)[0])]
    saved = 0
    while tail and head and tail[-1][0] == head[0][0]:
        a, b = tail.pop(), head.pop(0)
        v = tuple((x + y) % 4 for x, y in zip(a[1], b[1]))
        cost = _MIN_COST[v]
        saved += a[2] + b[2] - min(cost, a[2] + b[2])
        if v != (0, 0, 0):
            break
    return saved


def merge_tokens(first, second) -> tuple[str, ...]:
    """A literal move sequence for `first` then `second` with the boundary merged.

    Used to build real variant algorithms (eg swapping an insert behind a prefix);
    rotations must not sit at the boundary for this.
    """
    a, b = list(first), list(second)
    while a and b and a[-1][0] not in ROTATIONS and b[0][0] not in ROTATIONS:
        ax, va = layer_vector(a[-1])
        bx, vb = layer_vector(b[0])
        if ax != bx or a[-1][0] != b[0][0]:
            break
        k = (_POWER[a[-1][1:]] + _POWER[b[0][1:]]) % 4
        f = a[-1][0]
        a.pop()
        b.pop(0)
        if k:
            a.append(f + {1: "", 2: "2", 3: "'"}[k])
            break
    return tuple(a) + tuple(b)


def to_magiccube(tokens) -> str:
    """Spell tokens in magiccube's notation (Rw for r, etc.)."""
    out = []
    for t in tokens:
        if t[0] in "urfdlb":
            out.append(t[0].upper() + "w" + t[1:])
        else:
            out.append(t)
    return " ".join(out)
