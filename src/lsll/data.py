"""Loads the algorithm files, and applies the project's selection rules.

Baseline alg of a case: most votes; ties (including no votes at all) go to database
algs before generated ones, then the best (lowest) speed score.
Alternates of a case: every alg ranked at least as fast as the baseline, score <= base.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .cube import ROTATIONS, parse, perm_of, stm

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "cube-simulator-and-algs" / "data"


@dataclass(frozen=True)
class Alg:
    text: str
    tokens: tuple
    moves: int
    score: float
    votes: int
    src: str
    awk: bool


@dataclass
class Case:
    name: str
    group: str
    state: str
    algs: list


def upright(tokens):
    """Add a free rotation if the alg leaves the U centre off the top.

    A few listed algs finish with the cube on its side (a net x, z, or wide moves). The
    solver holds the cube upright again; rotations cost nothing, so the count is unchanged.
    """
    tokens = tuple(tokens)
    if perm_of(tokens)[4] == 4:
        return tokens
    for r in ("x", "x'", "x2", "z", "z'"):
        if perm_of(tokens + (r,))[4] == 4:
            return tokens + (r,)
    raise ValueError(tokens)


def _alg(a: dict) -> Alg:
    toks = upright(parse(a["alg"]))
    return Alg(a["alg"], toks, stm(toks), float(a.get("score", 1e9)), int(a.get("votes", 0) or 0), a.get("src", "db"), bool(a.get("awk")))


def load_set(name: str) -> list[Case]:
    d = json.loads((DATA / "ll" / f"{name}.json").read_text())
    return [Case(c["name"], c.get("group", ""), c["state"], [_alg(a) for a in c["algs"]]) for c in d["cases"]]


def baseline_key(a: Alg):
    return (-a.votes, a.src != "db", a.score)


def baseline_of(algs):
    return min(algs, key=baseline_key)


def alternates_of(algs, base: Alg):
    return [a for a in algs if a.score <= base.score and a.tokens != base.tokens]


def load_scdb_f2l():
    d = json.loads((DATA / "ll" / "scdb.json").read_text())
    return d["sets"]["f2l"]
