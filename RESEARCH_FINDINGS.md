# Research Findings: CFOP Last Slot + Last Layer

Everything here comes from the code in this repository. Rerun it with the commands in section 12. Numbers are exact unless marked as a sample.

**Metric:** STM. Every face, wide or slice turn counts 1 at any angle (U2 = 1), and rotations count 0.
**Scope:** FR last slot, then last layer.
**Method:** full enumeration of all 9,331,200 states, so results have no sampling error. Techniques are compared on identical states.

---

## 1. Headline results

| Question | Answer |
|---|---|
| Baseline movecount | Mean **32.43** (exactly 350237/10800), SD 4.83, median 33, range 0 to 49 |
| Where the moves go | Last pair 7.81, OLL 10.31, PLL 14.31 (AUFs included) |
| No new algorithms, no decision | Cancellations save 0.52 moves. They happen automatically |
| No new algorithms, with a decision | U2 + sledge inserts save 0.49 by the two-step rule, and 1.20 with full look-ahead (a ceiling; see section 11) |
| Everything in CFOP combined | 1.47 to 2.06 moves saved |
| ZBLS + ZBLL | **7.22 moves saved**, for 733 more algorithms |
| Most surprising result | Using "faster" PLL algs from other angles **adds** 0.39 moves |

**Central finding for the paper:** convolution is exact for the baseline, but only when each AUF stays grouped with its own algorithm. It stops being exact as soon as a decision or a cancellation crosses a step boundary. The mean survives that failure, but the spread and tails do not.

---

## 2. Data decisions

| Topic | Decision | Reason |
|---|---|---|
| Baseline alg | Most votes in the JSON. Ties go to database algs before generated ones, then to the lowest speed score | Your Q1 and Q2 |
| F2L source | `scdb.json`, set `f2l` | The only F2L data with votes |
| FR filter | An alg counts only if simulation shows it solves FR and keeps the other pairs | 458 of the 622 listed F2L algs are for other slots. Exactly 4 per case survive, which leaves 164 |
| AUF built into F2L algs | **Stripped.** `U R U' R'` becomes the body `R U' R'` plus a separately counted AUF | Otherwise AUFs would be double-counted. This was the trap flagged in `CLAUDE.md` |
| Rotations | Only those written in an alg; none added as setup | Your Q15 |
| Net rotation | If an alg ends with the cube on its side, it is turned upright for 0 moves | Rotations are free, so the count does not change |
| COLL | The edges-oriented subset of `coll.json`: **40 cases** | The file's 329 cases are really OLLCP. Those 40 cover every edges-oriented state with unoriented corners. The earlier "42 cases" note does not match this data |
| ZBLL convention | 472 ZBLL + 21 PLL + 1 skip = 494 | State this in the paper, as `CLAUDE.md` asks |
| OLS | 432 + 7,622 = 8,054 cases | Reference point only |

**Data check:** 7,488 listed OLL, PLL, COLL and ZBLL algs were simulated on every state of their case.
- 12 (0.16%) failed: 1 PLL, 2 COLL, 9 ZBLL. They finish with a turn that is not on the U layer once the cube is upright.
- None is a baseline alg, and the model never selects them.
- All 741 OLL algs passed.

---

## 3. The baseline

### Step distributions (`results/tables/*_baseline.csv` has every case)

| Step | Mean | SD |
|---|---:|---:|
| F2L alg | 7.087 | 1.908 |
| OLL alg | 9.579 | 1.987 |
| PLL alg | 12.875 | 3.855 |

| Block (AUFs inside) | Mean |
|---|---:|
| Last pair | 7.807 |
| OLL | 10.310 |
| PLL | 14.313 |
| **Total** | **32.429** |

**Thresholds:** P(T ≤ 20) = 0.0152, P(T ≤ 25) = 0.0818, P(T ≤ 30) = 0.3270, P(T ≤ 35) = 0.7177, P(T ≤ 40) = 0.9730.

### The AUF distribution is not exactly 0.25 / 0.75

You gave P(0 moves) = 0.25 for every AUF. Enumeration gives:

| AUF | Exact P(0 moves) | Decimal | Why it differs from 1/4 |
|---|---:|---:|---|
| Before last pair | 7/25 | 0.2800 | The 6 placements with both pieces in the slot never need one |
| Before OLL | 29/108 | 0.2685 | OLL 20 never needs one. OLL 1, 21, 55, 56, 57 need none half the time. An OLL skip needs none |
| Before PLL | 11/36 | 0.3056 | H, Na, Nb never need one. E and Z need none half the time. A PLL skip needs none |
| After PLL | 37/144 | 0.2569 | In E and Z the solver picks the angle that avoids a final AUF |

**N perms are the nice example:** each has only 4 states, which differ only by a final AUF, so they never need a pre-AUF.

Convolving with flat 0.25/0.75 AUFs gives mean 32.54 (true value 32.43), a maximum of 50 (true maximum 49), and total variation distance 0.0093 from the truth. The lower tail is off by 15%: P(T ≤ 20) is 0.0130 against the true 0.0152.

**Using the correct AUF probabilities does not fix this.** With exact AUF marginals the mean becomes right, but the distance stays at 0.0093. The error comes from dependence between an AUF and its alg: symmetric cases have their own lengths *and* need fewer AUFs. That makes a stronger point for the paper than "the number was slightly off".

---

## 4. Where convolution is valid (exact integer tests)

Each test compares joint counts with the product of marginals, in exact integers.

| Pair tested | Independent? |
|---|---|
| Last-pair block vs last-layer block | **Yes, exactly** |
| OLL block vs PLL block | **Yes, exactly** |
| AUF before OLL vs OLL alg length | No |
| AUF before PLL vs PLL alg length | No |
| AUF after PLL vs PLL alg length | No |

### Why the blocks are independent (a proof you can include)

Let G be the 62,208 last-layer states, each equally likely.

1. **Last pair is independent of last layer.** For a fixed placement and alg, the map from state before the alg to state after it is a bijection on G. A bijection sends the uniform distribution to the uniform distribution, whichever placement occurred. So the LL state tells you nothing about the last-pair block.
2. **OLL is independent of PLL.** G splits into 216 orientation patterns of 288 states each. The chosen alg depends only on the pattern, and it maps those 288 states one-to-one onto the 288 oriented states. So the PLL state is uniform whichever OLL pattern occurred.

This is the justification the outline asks for: the independence is proved from bijections, not assumed.

### Three computations of the baseline

| Comparison with enumeration | Total variation distance |
|---|---:|
| Direct convolution of blocks | 3.8 × 10⁻¹⁷ |
| DFT of blocks (N = 50) | 3.2 × 10⁻¹⁵ |
| Flat AUFs, 7 separate steps | 9.3 × 10⁻³ |

### Padding

The blocks have maxima 12, 15 and 22, so N ≥ 12 + 15 + 22 + 1 = **50**.

| N | Distance from truth | Mean if read naively |
|---|---:|---:|
| 20 | 0.989 | 11.68 |
| 32 | 0.610 | 12.90 |
| 40 | 0.049 | 30.48 |
| 50 | 3.2 × 10⁻¹⁵ | 32.43 |

N = 40 is the instructive one: it looks almost right and is quietly wrong.

### The 4-point hand example

Input: one flat AUF padded to length 4, a = [0.25, 0.75, 0, 0]. Three AUFs have maximum total 3, so N = 4 is enough. Roots are 1, −i, −1, i.

| Evaluation point | A | A³ |
|---|---|---|
| 1 | 1 | 1 |
| −i | 0.25 − 0.75i | −0.40625 + 0.28125i |
| −1 | −0.5 | −0.125 |
| i | 0.25 + 0.75i | −0.40625 − 0.28125i |

The inverse gives [1/64, 9/64, 27/64, 27/64], which is exactly Binomial(3, 0.75). Note that the value at z = 1 is always 1: it is the sum of the probabilities.

---

## 5. Which techniques break independence

The test: does convolving the technique's own step distributions reproduce its true distribution?

| Technique | Convolution error (TV) | Convolved SD | True SD | Behaviour |
|---|---:|---:|---:|---|
| Blind U2 or sledge insert | ~10⁻¹⁵ | same | same | Still exact |
| Speed-gated angles | ~10⁻¹⁵ | same | same | Still exact |
| Fewest-move angles (each step decided alone) | ~10⁻¹⁵ | same | same | Still exact |
| Cancellations | 0.0031 | 4.973 | 4.967 | Breaks slightly |
| U2 + sledge, two-step | 0.019 | 4.863 | 4.758 | Breaks |
| COLL, two-step | 0.013 | 5.147 | 5.039 | Breaks |
| OLL angles, two-step | 0.030 | 5.000 | 4.776 | Breaks |
| U2 + sledge, full look-ahead | 0.037 | 5.051 | 4.742 | Breaks |
| OLLCP, two-step | 0.073 | 4.595 | 3.905 | Breaks badly (fig 9) |
| ZBLS + ZBLL | 0.072 | 5.493 | 3.506 | Breaks badly |
| ZBLL only | 0.110 | 6.842 | 5.509 | Breaks badly |

**This corrects the project's working assumption.** `CLAUDE.md` says choice rules, cancellations *and* alternate angles all link adjacent steps. Angles chosen from the current step alone do **not**: the OLL alg changes, but the PLL state stays uniform. The exact dividing line is:

> Convolution of blocks stays exact unless a decision or a merge uses information from across a step boundary.

The mean always stays correct, because of linearity. That is exactly why full distributions matter: a means-only analysis could never detect this failure. Convolution overstates the spread when a choice is involved, because a choice picks whichever step is short and so makes the steps negatively correlated.

### Conditioning repairs it (mix, then convolve)

For COLL (used whenever available), split the LL states into three classes and convolve inside each:

| Class | Weight | Mean moves (LL) | Steps independent within class? |
|---|---:|---:|---|
| Edges not oriented: OLL, then PLL | 7/8 | 24.84 | Yes |
| Edges oriented, corners not: COLL, then EPLL | 13/108 | 20.00 | Yes |
| All oriented: PLL only | 1/216 | 14.31 | Yes |

The mixture matches enumeration to **4 × 10⁻¹⁵**, while plain convolution is off by 0.026. So the DFT is not abandoned here: it still works one level down, within each class. That gives the paper a stronger ending than "the DFT failed".

---

## 6. All techniques

| Technique | Mean | SD | Saved | Fewer / same / more | P(T≤25) | P(T≤30) | LS | LS+OLL | OLL+PLL | Net algs | Saved per 100 algs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Baseline | 32.43 | 4.83 | +0.00 | 0.0% / 100.0% / 0.0% | 0.082 | 0.327 | 7.81 | 18.12 | 24.62 | 0 |  |
| U2 insert, always | 32.48 | 4.85 | -0.05 | 15.4% / 69.0% / 15.6% | 0.081 | 0.324 | 7.86 | 18.17 | 24.62 | 0 |  |
| Sledge insert, always | 32.51 | 4.85 | -0.08 | 15.0% / 69.0% / 16.0% | 0.081 | 0.322 | 7.89 | 18.20 | 24.62 | 0 |  |
| U2 insert, two-step | 32.14 | 4.78 | +0.29 | 8.8% / 87.7% / 3.6% | 0.089 | 0.344 | 7.82 | 17.83 | 24.32 | 0 |  |
| Sledge insert, two-step | 32.09 | 4.76 | +0.34 | 8.9% / 87.4% / 3.7% | 0.090 | 0.346 | 7.77 | 17.78 | 24.32 | 0 |  |
| U2 + sledge, two-step | 31.94 | 4.76 | +0.49 | 13.0% / 82.1% / 4.9% | 0.095 | 0.356 | 7.79 | 17.62 | 24.15 | 0 |  |
| U2 + sledge, full look-ahead | 31.23 | 4.74 | +1.20 | 21.2% / 78.8% / 0.0% | 0.113 | 0.421 | 7.82 | 17.91 | 23.41 | 0 |  |
| Cancellations | 31.91 | 4.97 | +0.52 | 24.0% / 76.0% / 0.0% | 0.105 | 0.365 | 7.81 | 17.79 | 24.11 | 0 |  |
| Angles OLL, speed-gated | 32.32 | 4.83 | +0.11 | 10.7% / 82.5% / 6.8% | 0.086 | 0.333 | 7.81 | 18.01 | 24.51 | 42 | 0.26 |
| Angles PLL, speed-gated | 32.82 | 4.85 | -0.39 | 3.8% / 77.8% / 18.4% | 0.072 | 0.297 | 7.81 | 18.12 | 25.02 | 18 | -2.18 |
| Angles both, speed-gated | 32.71 | 4.85 | -0.28 | 12.9% / 64.6% / 22.6% | 0.075 | 0.303 | 7.81 | 18.01 | 24.90 | 60 | -0.47 |
| Angles both, speed-gated, no awkward | 32.74 | 4.85 | -0.31 | 12.0% / 65.5% / 22.5% | 0.074 | 0.301 | 7.81 | 18.04 | 24.94 | 57 | -0.55 |
| Angles both, fewest moves | 31.88 | 4.73 | +0.55 | 20.9% / 75.4% / 3.7% | 0.091 | 0.367 | 7.81 | 17.86 | 24.08 | 48 | 1.14 |
| Angles OLL, two-step | 31.65 | 4.78 | +0.78 | 17.1% / 82.9% / 0.0% | 0.100 | 0.386 | 7.81 | 18.05 | 23.84 | 42 | 1.85 |
| COLL, always when available | 32.02 | 5.02 | +0.41 | 7.9% / 90.5% / 1.6% | 0.101 | 0.368 | 7.81 | 18.46 | 24.21 | 40 | 1.03 |
| COLL, two-step | 31.97 | 5.04 | +0.46 | 7.9% / 92.1% / 0.0% | 0.103 | 0.374 | 7.81 | 18.34 | 24.16 | 40 | 1.16 |
| OLLCP, two-step | 28.50 | 3.90 | +3.93 | 67.5% / 32.5% / 0.0% | 0.190 | 0.691 | 7.81 | 19.83 | 20.69 | 329 | 1.19 |
| ZBLL only (no edge control) | 31.52 | 5.51 | +0.91 | 11.0% / 88.5% / 0.5% | 0.153 | 0.399 | 7.81 | 18.93 | 23.72 | 472 | 0.19 |
| ZBLS + ZBLL | 25.21 | 3.51 | +7.22 | 87.8% / 4.0% / 8.1% | 0.517 | 0.944 | 9.31 | 24.60 | 15.90 | 733 | 0.99 |
| ZBLS (shortest) + ZBLL | 24.36 | 3.22 | +8.07 | 90.7% / 3.4% / 5.9% | 0.621 | 0.983 | 8.46 | 23.75 | 15.90 | 733 | 1.10 |
| ZBLS + ZBLL + cancellations | 24.93 | 3.63 | +7.50 | 88.6% / 3.7% / 7.7% | 0.545 | 0.947 | 9.31 | 24.33 | 15.63 | 733 | 1.02 |
| OLS + PLL | 26.78 | 4.41 | +5.65 | 79.8% / 5.9% / 14.2% | 0.353 | 0.793 | 12.40 | 12.47 | 14.38 | 8,013 | 0.07 |
| 1LLL | 23.73 | 3.03 | +8.70 | 95.6% / 3.0% / 1.5% | 0.719 | 0.990 | 7.81 | 23.73 | 15.93 | 3,837 | 0.23 |
| All CFOP techniques, speed-gated | 30.96 | 5.01 | +1.47 | 46.7% / 36.6% / 16.7% | 0.141 | 0.438 | 7.79 | 17.18 | 23.17 | 100 | 1.47 |
| All CFOP techniques, fewest moves | 30.37 | 4.92 | +2.06 | 51.4% / 43.0% / 5.6% | 0.161 | 0.492 | 7.79 | 17.13 | 22.58 | 88 | 2.34 |

"Fewer / same / more" is the paired comparison with the baseline on the same state. "Net algs" counts what must be learned beyond the baseline's 119, with replaced algs subtracted. U2/sledge inserts and cancellations add no algs, though they do add recognition work. The LS+OLL column is the cumulative mean through the orientation step. For ZB and 1LLL that step is the one-look alg, so that column is not comparable with CFOP.

### Inserts

- Only **14 of 41** F2L cases (50 of 150 placements, so 1/3 of states) have a baseline alg ending in one of your three insert forms. The rest have no swap under your definition. No baseline alg ends in a mirrored `F' U F` form, so mirrors add nothing.
- Blind swapping costs moves: "always U2" adds 0.05, "always sledge" adds 0.08. The technique's value is entirely in the choice.
- Because the variants merge with what comes before, swaps are often cheaper than they look. `R U R' U2 + U R U2 R'` becomes `R U R' U' R U2 R'`.
- U2 and sledge are **complementary**. Separately they save 0.29 and 0.34; together, 0.49.

### The two-step rule is greedy, and it shows

- The two-step rule saves 0.49. Full look-ahead through PLL saves **1.20**, so it leaves 0.71 moves on the table.
- It picks suboptimally in **14.8%** of states.
- It does **worse than making no choice at all in 4.9%** of states, by up to 21 moves. Full look-ahead never does (fig 8).
- The cause: the rule treats PLL as irrelevant, but the insert changes which PLL you get.

### Alternate angles: faster is not fewer

| | Qualifying alternates | Mean STM vs baseline | Longer / same / shorter |
|---|---:|---:|---|
| OLL | 143 (24 cases) | +1.50 | 55% / 31% / 14% |
| PLL | 77 (11 cases) | +1.94 | 60% / 25% / 16% |

Under your rule (use a faster-or-equal alg when you arrive at its angle):

| Rule | OLL | PLL | Both |
|---|---:|---:|---:|
| Speed-gated (saved) | +0.11 | **−0.39** | −0.28 |

If the solver instead takes a qualifying alg only when it is actually shorter, the saving is +0.55 (both steps). Reading OLL's next step as well gives +0.78 (OLL only).

Your speed ranking is doing what you said it does, choosing faster algs, and that is a real result rather than a flaw: **this study measures moves, and the speed ranking is optimising something else.** The pro data agrees (section 8). `results/tables/alternates.csv` lists every alternate.

### Cancellations

- Mean saving 0.5161. In 76% of states nothing cancels; the maximum is 18 moves.
- **Bilinear prediction:** because the blocks are independent, E[saving] = Σ P(a) P(b) J(a, b) over segment pairs, where J(a, b) is the saving when b follows a. That gives 0.33114 (LS/OLL) + 0.18326 (OLL/PLL) + 0.00174 (LS/PLL on an OLL skip) = **0.51614**, matching enumeration to 14 decimal places. The independence result lets you predict the *mean* cancellation exactly without simulating any solves.
- Almost all of it is one pattern: `R'` ending an alg followed by `R` starting the next. That alone is worth 0.20 of the 0.33 at the LS/OLL boundary.
- AUF merges into the alg next to them contribute zero, since the data has already stripped those.
- **One known gap, measured:** the model merges neighbouring segments in pairs. If a whole OLL alg cancels away, the last pair and PLL can then merge through it, and the model misses that. It happens in 1 state in about 61,000 and is worth 0.00004 moves (corrected mean 0.51618). Worst case: 21 counted, 17 true. So cancellation results are exact to about 4 decimal places, not perfectly exact. All other results are unaffected.

### Interactions (`results/tables/factorial.csv`)

Without angles and cancellations together, combinations are nearly additive (interactions between −0.09 and +0.01). Every combination with both has a positive interaction (+0.04 to +0.46), largest here:

| | Saved |
|---|---:|
| Speed-gated angles alone | −0.28 |
| Cancellations alone | +0.52 |
| Sum expected | +0.23 |
| **Actually together** | **+0.70** (interaction +0.46) |

**Mechanism, tested:** a non-zero AUF sits between two algs and blocks them from cancelling. Speed-gated angles raise the chance of no AUF before OLL from 26.9% to 46.3%, and before PLL from 30.6% to 55.6%. Cancellation savings go from 0.52 to 0.98. So the speed-gated rule looks bad alone and much better in a real solve.

### COLL and ZB

- COLL (40 algs) saves 0.46 by two-step, taken in 7.9% of states. Using it whenever available saves 0.41 and sometimes costs moves (1.6% of states).
- Without edge control, ZBLL (472 algs) saves only 0.91, because it is usable in 1/8 of solves. ZBLS + ZBLL together save 7.22. **Most of ZB's value comes from guaranteeing the good state**, at a cost of 1.50 more moves in the last pair (9.31 vs 7.81).
- OLLCP (329 algs) saves 3.93, the best value per alg of any large set: 1.19 moves per 100 algs. COLL gets 1.16, ZB 0.99.
- **ZB is less sure of a win than its mean suggests:** it is worse than baseline in 8.1% of states, and P(T ≤ 25) rises from 8% to 52%.

**Caution on ZBLS and OLLCP:** these sets have almost no votes, so "standard" falls back to the lowest speed score. Picking the shortest listed ZBLS instead gives 8.07 saved, not 7.22 (last pair 8.46, not 9.31). So 7.22 reflects a realistic-speed choice, and 8.07 is the best the listed ZBLS algs can do.

---

## 7. Your own examples, verified

Do these on a real cube. Each scramble is applied to a solved cube; `results/examples.json` has every option in full.

**A longer insert that skips OLL.** Scramble:
`R2 F R U R U' R' F' R U2 R' U2 R U2 F R' F' R U R U' R' R U R' U' R U' R' U`

| Last pair | OLL | PLL | Total |
|---|---|---|---:|
| `U' R U R' U R U' R'` | `R U R' U' R' F R F'` | `U2 R' U2 R U2 R' F R U R' U' R' F' R2` | 30 |
| `U' R U R2 F R F'` (sledge) | skipped | `U2 R' U2 R U2 R' F R U R' U' R' F' R2` | **21** |

This one is even better than your "2 moves longer" example: the sledge merges and comes out a move shorter.

**A U2 insert worth its extra move.** Scramble:
`U' M2 U' M U2 M' U' M2 U F U R U' R' F' R U2 R' U' R U' R' R U R' F R' F' R U`

| Last pair | OLL | PLL | Total |
|---|---|---|---:|
| `U' R' F R F' R U' R'` | `R U R' U R U2 R' F R U R' U' F'` | `U' M2 U M U2 M' U M2 U` | 30 |
| `U' R' F R F' U R U2 R'` | `F R U R' U' F'` | `U' M2 U M U2 M' U M2 U` | **24** |

**The two-step rule losing.** Scramble:
`U' M2 U2 M U M2 U M2 U M R U2 R2 U' F' U F R2 U' R' U2 R U R' F R U R' U' F' U'`

| Last pair | Through OLL | Total |
|---|---:|---:|
| `U F U R U' R' F' R U' R'` (standard) | 21 | **31** |
| `U F U R U' R' F' U R U2 R'` (rule picks this) | 20 | 37 |

The rule saves one move through OLL and pays seven at PLL.

**Largest cancellation.** OLL 37 (`F R' F' R U R U' R'`) directly into a T perm (`R U R' U' R' F R2 U' R' U' R U R' F'`): 22 moves written, 9 performed, 13 cancelled. Chance of that exact pair with no AUFs: 1 in 15,552.

---

## 8. Against real solves

Data: reco.nz reconstructions. 862 solves in `reco_lsll.txt` and 15,795 LL rows in `pro_ll.json`. **The last slot was FR in only 366 of the 862 (42%)**: FL 278, BR 125, BL 93. Other slots are viewed through a y rotation, with FR-only results shown separately.

### How close the model is

| Quantity | Model | Pros |
|---|---:|---:|
| Last pair, AUF included | 7.81 | 7.48 ± 0.09 (n = 582) |
| Same, FR only | 7.81 | 7.40 ± 0.14 (n = 225) |
| OLL alg | 9.62 | 9.26 ± 0.03 (n = 6,040) |
| PLL alg (non-skip) | 13.06 | 13.14 ± 0.04 (n = 7,928) |

The model and the pros are within about 0.4 moves on each step. (OLL and PLL compare alg lengths without AUFs, since the rows do not record them.)

### The uniform state assumption does not hold for pros

| Test | Model | Pros | p |
|---|---:|---:|---:|
| Edges oriented after last pair (labelled CFOP) | 12.5% | 20.1% (n = 582) | 2 × 10⁻⁷ |
| Same, FR only | 12.5% | 22.2% (n = 225) | 5 × 10⁻⁵ |
| F2L case frequencies uniform | uniform | no (χ² = 90.7) | 1 × 10⁻⁵ |
| Same, FR only | uniform | no (χ² = 64.0) | 0.012 |
| OLL case frequencies uniform | uniform | no (χ² = 1429) | < 10⁻²⁰⁰ |
| PLL case frequencies uniform | uniform | no (χ² = 499) | < 10⁻⁹⁰ |

Easy last pairs (F2L 1, 2, 3, 30), Sunes (OLL 26, 27) and U perms are over-represented.

**Two causes, which this data cannot separate:**
1. Real skill: pros choose their last pair, and control edges.
2. Selection: reconstructed solves are disproportionately fast and record solves, which tend to have easy cases.

**Pros use longer algs than the baseline's.** At the pros' actual case frequencies, baseline algs would average 9.01 (OLL) and 12.41 (PLL), but pros actually use 9.26 and 13.14. Only 58% of OLL rows and 42% of PLL rows use the exact baseline alg. This backs up "faster is not fewer": the best solvers in the world knowingly pay 0.73 extra PLL moves.

**What to say in the paper:** the model describes an idealised solver on uniformly random states. It matches pros closely on step lengths but not on which states they face. Comparisons *between techniques* are more trustworthy than the absolute totals, since each technique is measured on identical states.

---

## 9. How the model was checked

| Check | Result |
|---|---|
| Cube engine vs magiccube (a separate library) | All 54 single-move tokens and 300 random 40-move sequences identical |
| Independent brute-force solver, state by state, baseline | **0 mismatches** in 20,000 random states |
| Same, two-step rule | **0 mismatches** in 5,000 |
| Same, full look-ahead | **0 mismatches** in 5,000 |
| Cancellations, separate recount | 5 differences in 20,000: 4 are flaws in the recount (it double-adds a merge inside an alg), 1 is the three-way case in section 6 |
| Sample distribution vs exact | χ² p = 0.33; sample mean 32.47 vs exact 32.43 |
| Random states reach all placements evenly | χ² p = 0.33, all 150 seen |
| Whole solves replayed in magiccube | 1000 of 1000 solved |
| Bilinear cancellation formula vs enumeration | Agree to 14 decimal places |
| Automated tests | `python3 -m pytest tests` |

The brute-force solver recognises cases by trial on raw stickers and shares no case tables with the model. It also found two things during development: a real bug, where a long cancellation was cut off by a window that was too short (OLL 37 into T perm), now fixed; and the three-way gap in section 6, now measured.

---

## 10. Assumptions you can overrule

These were the obvious readings, and each is a small code change.

1. **Two-step scope.** Your answers did not settle where the rule applies (`CLAUDE.md` open question 5). It is applied wherever a real decision exists: insert choice, COLL vs OLL, OLL alternates. Each is reported separately.
2. **Tie-breaks.** If options tie, keep the baseline, then the earliest listed option.
3. **Speed-gated rule.** At the baseline's own angle, keep the baseline. From another angle, take the fastest qualifying alg that needs no AUF.
4. **AWK-flagged algs stay in,** because your rule was rank only. Excluding them barely changes things (−0.31 vs −0.28).
5. **Generated algs count as alternates,** following your "more searched algorithms".
6. **Solved pair, edges not oriented, under ZB:** solve with OLL + PLL, since no ZBLS case exists for it.
7. **Pre-AUF choice when several work:** no AUF if possible, then U or U' before U2 (all cost 1).
8. **Cancellations** merge only at boundaries; written algs count as written. Slice merges such as `R' + r → M'` count, as STM requires. `R + L'` stays 2 moves, even though it equals a slice plus a rotation.

---

## 11. Limitations to reflect on

1. **Movecount is not time.** The clearest evidence is in this study itself: speed-ranked algs increase the count, and pros use longer algs by choice.
2. **Recognition is free in this model.** The two-step rule needs the solver to read the OLL before inserting, and full look-ahead needs them to read the PLL too. The 1.20 figure is a ceiling, not something a person can do in a solve.
3. **Uniform states are false for real solves** (section 8), and that difference would change the numbers.
4. **Only baseline F2L algs define where inserts apply.** A solver could pick a different alg in the other 27 cases to enable a swap, so 0.49 is conservative.
5. **Only the listed algs exist.** All "fewest moves" results are limited to what is in the dataset.
6. **"Most upvoted" is a proxy for "standard".** It reflects a community, and ZBLS and OLLCP barely have votes.
7. **Scope excludes the cross and first three pairs.** So the last pair cannot be chosen here, while real solvers do choose it.
8. **Greedy is not optimal even in the "full" version.** Full look-ahead here covers the insert decision only, with later choices fixed.

---

## 12. Reproduce

```
python3 -m pytest tests                  # about 20 s
python3 scripts/run_all.py               # about 11 min: results/results.json, results/tables/
python3 scripts/real_solves.py           # results/real_solves.json
python3 scripts/make_figures.py          # results/figures/
python3 scripts/examples.py              # results/examples.json
python3 scripts/validate.py 20000        # about 5 min: results/validation.json
python3 scripts/triple_cancellation.py   # about 1 min: results/triple_cancellation.json
```

| Figure | Shows | Paper section |
|---|---|---|
| `fig1_blocks` | The three block distributions | 3 |
| `fig2_flat_auf_error` | Enumeration vs flat-AUF convolution, and the error | 3 |
| `fig3_wraparound` | N = 32 vs N = 50 | 4 |
| `fig9_independence_failure` | Right mean, wrong tails | 5 |
| `fig8_paired_saving` | Two-step vs full look-ahead, state by state | 5 |
| `fig6_moves_saved` | Saving by technique | 6 |
| `fig4_overlay_pmf`, `fig5_overlay_cdf` | Overlay of four systems | 6 |
| `fig7_tradeoff` | Moves saved vs algs learned | 6 |
| `fig10_last_pair_vs_real` | Model vs pros, last pair | 7 |
