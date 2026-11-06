"""Distributions on movecounts: direct convolution, the DFT (hand-written), statistics."""
from __future__ import annotations

from fractions import Fraction

import numpy as np


def normalise(counts) -> np.ndarray:
    c = np.asarray(counts, dtype=float)
    return c / c.sum()


def trim(p, eps=0.0):
    p = np.asarray(p)
    nz = np.nonzero(np.abs(p) > eps)[0]
    return p[: nz[-1] + 1] if len(nz) else p[:1]


def convolve_direct(*dists) -> np.ndarray:
    out = np.array([1.0])
    for d in dists:
        res = np.zeros(len(out) + len(d) - 1)
        for i, a in enumerate(out):
            for j, b in enumerate(d):
                res[i + j] += a * b
        out = res
    return out


def dft(a, sign=-1) -> np.ndarray:
    """Evaluate the polynomial with coefficients `a` at the N-th roots of unity."""
    a = np.asarray(a, dtype=complex)
    n = len(a)
    w = np.exp(sign * 2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n)
    return w @ a


def idft(v) -> np.ndarray:
    """Inverse DFT: same sum with the conjugate root, then divide by N."""
    return dft(v, sign=+1) / len(v)


def convolve_dft(*dists, n=None) -> np.ndarray:
    """All distributions at once: pad to N, evaluate, multiply pointwise, invert."""
    need = sum(len(d) - 1 for d in dists) + 1
    n = need if n is None else n
    prod = np.ones(n, dtype=complex)
    for d in dists:
        pad = np.zeros(n)
        if len(d) > n:  # a too-short N: coefficients wrap around
            for i, x in enumerate(d):
                pad[i % n] += x
        else:
            pad[: len(d)] = d
        prod = prod * dft(pad)
    return idft(prod).real


def stats(p, thresholds=()) -> dict:
    p = np.asarray(p, dtype=float)
    x = np.arange(len(p))
    mean = float(p @ x)
    var = float(p @ (x - mean) ** 2)
    cdf = np.cumsum(p)
    out = {
        "mean": mean,
        "sd": var**0.5,
        "min": int(np.nonzero(p > 0)[0][0]),
        "max": int(np.nonzero(p > 0)[0][-1]),
        "median": int(np.searchsorted(cdf, 0.5)),
        "q10": int(np.searchsorted(cdf, 0.1)),
        "q90": int(np.searchsorted(cdf, 0.9)),
    }
    for m in thresholds:
        out[f"P(T<={m})"] = float(cdf[m]) if m < len(cdf) else 1.0
    return out


def tv_distance(p, q) -> float:
    n = max(len(p), len(q))
    a, b = np.zeros(n), np.zeros(n)
    a[: len(p)] = p
    b[: len(q)] = q
    return float(0.5 * np.abs(a - b).sum())


def exact_mean(counts) -> Fraction:
    total = sum(int(c) for c in counts)
    return Fraction(sum(i * int(c) for i, c in enumerate(counts)), total)
