import sys, random
sys.path.insert(0, "src")
import numpy as np
import magiccube
from lsll.cube import *
from lsll import cube as C

ALL = [b + s for b in "URFDLBurfdlbMESxyz" for s in ("", "'", "2")]
FACE = "URFDLB"
SOLVED = np.repeat(np.arange(6), 9)

def lab(tokens):
    return "".join(FACE[i] for i in SOLVED[perm_of(tuple(tokens))])

def mc(tokens):
    c = magiccube.Cube(3)
    if tokens:
        c.rotate(to_magiccube(tokens))
    return c

def norm_mc(c):
    # map magiccube's colors to labels through its center stickers of a solved cube
    s = c.get_kociemba_facelet_colors()
    return s

def test_stm():
    assert stm("R U2 R' F R F' U2") == 7
    assert stm("x R' U R' D2 R U' R' D2 R2 x'") == 9
    assert stm("r U R' U' M") == 5
    assert stm("y2 Rw U") == 2
    assert stm("(R U R') [U' R']") == 5

def test_group_orders():
    for t, o in [("U", 4), ("M", 4), ("r", 4), ("x", 4), ("S", 4), ("E", 4), ("d", 4)]:
        p = IDENTITY
        for _ in range(o):
            p = compose(p, PERM[t])
        assert (p == IDENTITY).all(), t
    assert (perm_of(parse("R U R' U'") * 6) == IDENTITY).all()

def test_derived_definitions():
    # wide = face + slice, in standard notation
    assert (perm_of(("r",)) == perm_of(("R", "M'"))).all()
    assert (perm_of(("l",)) == perm_of(("L", "M"))).all()
    assert (perm_of(("u",)) == perm_of(("U", "E'"))).all()
    assert (perm_of(("d",)) == perm_of(("D", "E"))).all()
    assert (perm_of(("f",)) == perm_of(("F", "S"))).all()
    assert (perm_of(("b",)) == perm_of(("B", "S'"))).all()
    assert (perm_of(("x",)) == perm_of(("R", "M'", "L'"))).all()
    assert (perm_of(("y",)) == perm_of(("U", "E'", "D'"))).all()
    assert (perm_of(("z",)) == perm_of(("F", "S", "B'"))).all()

def test_layer_vectors_match_perms():
    # the layer-vector algebra used for cancellations must agree with the permutations
    layer = {"x": ("R", "M'", "L'"), "y": ("U", "E'", "D'"), "z": ("F", "S", "B'")}
    for t in ALL:
        if t[0] in "xyz":
            continue
        axis, v = layer_vector(t)
        p = IDENTITY
        for gen, k in zip(layer[axis], v):
            for _ in range(k):
                p = compose(p, PERM[gen])
        assert (p == PERM[t]).all(), t

def test_against_magiccube():
    """Independent oracle: every token and random sequences match magiccube."""
    rnd = random.Random(1)
    # learn the fixed facelet mapping from a solved cube: identical ordering, color names differ
    solved = mc([]).get_kociemba_facelet_colors()
    cmap = {}
    for i, ch in enumerate(solved):
        cmap[ch] = FACE[i // 9]
    for seq in [[t] for t in ALL] + [[rnd.choice(ALL) for _ in range(40)] for _ in range(300)]:
        got = "".join(cmap[ch] for ch in mc(seq).get_kociemba_facelet_colors())
        assert got == lab(seq), (seq, got, lab(seq))

def test_join_savings():
    p = parse
    assert join_savings(p("R U"), p("U R")) == 1          # U U -> U2
    assert join_savings(p("R U"), p("U' R")) == 3         # U U' -> nothing; R R -> R2
    assert join_savings(p("R U"), p("U' R'")) == 4        # all gone
    assert join_savings(p("R U2"), p("U R")) == 1
    assert join_savings(p("F R"), p("U R")) == 0
    assert join_savings(p("R'"), p("r U")) == 1           # R' r = M'
    assert join_savings(p("R y"), p("F'")) == 2           # y moves the F' onto the R
    assert join_savings(p("R y'"), p("F'")) == 0
    assert join_savings(p("U x"), p("U")) == 0

def test_join_is_physically_correct():
    # when savings are claimed via a rotation, the literal perms really do cancel
    a, b = parse("R y"), parse("F'")
    assert (perm_of(a + b) == perm_of(("y",))).all()
