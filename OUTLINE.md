# CFOP Last-Slot-Last-Layer Paper: Outline

Rebuilt after the research phase. Every number below is an exact result from `RESEARCH_FINDINGS.md`, which has the full tables, the proofs and the caveats. Anything only you can supply is marked `[NEEDED FROM AUTHOR]`.

## Research question

**How can probability distributions be used to compare how much last-slot and last-layer techniques reduce the movecount of a CFOP solve?**

This holds up well, so I kept your working version. The results would also support a sharper variant if you prefer it: *"...and when can those distributions be combined by convolution?"*

## Central argument

A CFOP last slot and last layer can be treated as a sum of random movecounts, and if the steps are independent the DFT combines them in one operation, which is the tool you set out to learn.

1. **Independence is provable, but in a specific form.** Whole blocks (an AUF together with its alg) are exactly independent, by a bijection argument. Separate AUFs are not, so the natural seven-step model is slightly wrong however good its AUF probabilities are.
2. **The DFT is exact for the baseline** (mean 32.43, SD 4.83), and the interesting techniques then break it, along a precise line: convolution fails when a decision or a merge uses information from across a step boundary.
3. **The failure is invisible in the mean.** Means add regardless. The damage is in the spread and the tails, which is why full distributions are needed. For OLLCP, naive convolution gives SD 4.60 against a true 3.90.
4. **The DFT comes back** when you condition on the state and mix, matching enumeration to about 10⁻¹⁵.
5. **Compared state by state,** the results challenge common intuitions: faster algs can add moves (PLL angles, +0.39), a greedy rule can lose to making no choice (4.9% of states), and ZB's advantage (7.22) comes mainly from guaranteeing the good state.

## Page budget

| Section | Pages | Change from suggested budget |
|---|---:|---|
| 1. Introduction, aim, prediction | 1 | |
| 2. Scope, definitions, model | 2 | |
| 3. Baseline distributions | 2 | |
| 4. Combining distributions: polynomials and DFT | 3.5 | |
| 5. Techniques and the failure of independence | 3.5 | +0.5: this is the analytic core |
| 6. Results and comparison | 1.5 | −0.5: Section 5 already interprets most results, so this is one table and two figures |
| 7. Validation, reflection, conclusion | 2 | |
| **Total** | **15.5** | Appendix tables do not count towards the argument, but check whether your limit includes them |

## Still needed from you

Everything that testing could settle is settled. These are left because they are yours:

| Item | Why it can't come from me |
|---|---|
| `[NEEDED FROM AUTHOR]` Your honest prediction, and its timing | Section 1. Must be what you actually expected. See the note there |
| `[NEEDED FROM AUTHOR]` Your story of chaining convolutions and learning the DFT | Section 1. Your voice and your experience |
| `[NEEDED FROM AUTHOR]` A hand-worked version of the 4-point DFT | Section 4. Answer checks are in the findings; the working must be yours |
| `[NEEDED FROM AUTHOR]` Confidence with the bridge from your table method to the summation formula | Section 4. `CLAUDE.md` flags this as something you must understand before submitting |
| `[NEEDED FROM AUTHOR]` Line spacing, citation style, whether appendices count towards length | Assignment rules I don't have |
| `[NEEDED FROM AUTHOR]` Which of the optional items to cut if over length | Judgement call; suggested cuts are in Section 6 |

Two decisions where I picked a default that you may want to change (`RESEARCH_FINDINGS.md` section 10 lists all assumptions):
- Reporting speed-gated angles *and* fewest-move angles side by side, rather than one as the headline.
- Applying the two-step rule at every real decision: insert, COLL vs OLL, and OLL alternates.

---

## Section 1: Introduction, aim, prediction

**Purpose and claim:** Establishes why the question matters to you and what will be measured. Claims that the comparison needs whole distributions, not just averages.

**Pages:** 1

**Key content:**
- Your dissatisfaction with chaining pairwise convolutions, leading to the DFT.
- The research question.
- Scope in one sentence: FR last slot and last layer, in STM.
- Why these steps: your stated view that the cross and first three pairs are not deterministic enough to model well.
- Movecount is not time, and this will come back.

**Math shown vs. summarized:** None yet, beyond naming T as total moves.

**Visuals:** Optional diagram of the seven-step sequence. Don't spend much space on it.

**Data needed:** `[NEEDED FROM AUTHOR]` your motivation story.

**Personal engagement:** Your real story, and a prediction. Natural things to predict are which technique saves most per alg, and whether "two-step" would ever do worse than no choice. **Only state a prediction you truly made before seeing these results.** If you did not make one, say so, and instead describe what you would have guessed and which result surprised you most. An honest surprise earns more than a retrofitted prediction, and an examiner may probe it.

**Reflection:** Flag that the DFT was chosen for structural reasons (see Section 4), so you are not claiming a speed benefit at these sizes.

**Rubric:** Organization (clear aim), Personal engagement.

**Transition:** To compare anything, a solve must first become something countable.

---

## Section 2: Scope, definitions, model

**Purpose and claim:** Fixes the metric, state space, algorithm rule and the techniques. Claims that the whole problem can be enumerated exactly, so results carry no sampling error.

**Pages:** 2

**Key content:**
- **STM:** every face, wide or slice turn is 1 (U2 = 1); rotations are 0.
- **State space, as a counting argument:** 62,208 last-layer states (4! · 3³ · 4! · 2³ / 2) × 150 FR pair placements (5 · 3 · 5 · 2) = **9,331,200**. Also 216 orientation patterns × 288 permutations = 62,208, which is used later.
- **Baseline sequence:** the seven steps.
- **Alg rule:** most upvoted; alternates are those ranked at least as fast.
- **Technique table:** U2 insert (`U R U2 R'`), sledge (`U' R' F R F'`), both replacing `R U' R'`; angles; cancellations; COLL; ZBLS + ZBLL; plus the two-step rule as the decision method.
- **Uniform assumption,** stated plainly and marked for testing in Section 7.

**Math shown vs. summarized:** Show the state count in full. Summarize the data processing in a sentence each.

**Visuals:** Technique table (one row each: definition, what it changes, how it's modelled). Small example of `R U' R'` vs `U R U2 R'` from one pair state.

**Data needed:** All available. Source: `RESEARCH_FINDINGS.md` section 2.

**Personal engagement:**
- Finding the double-count: SpeedCubeDB writes `U R U' R'`, with the AUF inside the alg. Stripping that was a real modelling decision.
- Finding that only 164 of 622 listed F2L algs are for FR, and that simulation rather than labels is what told them apart.
- Your choice of FR only. Pros end at FR 42% of the time.

**Reflection:** Limits of "most upvoted" as a definition of "standard": it represents a community, and ZBLS and OLLCP barely have votes. Enumeration replaces sampling, but only for the idealized solver.

**Rubric:** Mathematical communication, Organization.

**Transition:** With states and algs fixed, each step has a distribution.

---

## Section 3: Baseline distributions

**Purpose and claim:** Builds the step distributions. Claims that even the "simple" AUF has structure that a flat 0.25/0.75 misses, because of symmetry.

**Pages:** 2

**Key content:**
- Case probabilities come from symmetry: OLL skip 1/216, PLL skip 1/72. N perm and H perm have 4 states, E and Z have 8, others 16.
- Means: F2L alg 7.09, OLL alg 9.58, PLL alg 12.88.
- **Your flat AUF vs exact:**

| AUF | Exact P(0) | Decimal |
|---|---:|---:|
| Before last pair | 7/25 | 0.2800 |
| Before OLL | 29/108 | 0.2685 |
| Before PLL | 11/36 | 0.3056 |
| After PLL | 37/144 | 0.2569 |

- One worked explanation: why an N perm never needs a pre-AUF but a T perm does 3/4 of the time.
- Blocks: last pair 7.81, OLL 10.31, PLL 14.31.
- Means add: 32.43, exactly 350237/10800.

**Math shown vs. summarized:** Show one AUF probability derived fully (7/25 is the simplest: 6 of 150 placements always need none, and 1/4 of the other 144 need none, giving 42/150). Summarize the rest in the table. Full case tables go in the appendix (`results/tables/`).

**Visuals:** `fig1_blocks` (x: moves, y: probability) and the AUF table.

**Data needed:** `results/tables/f2l_baseline.csv`, `oll_baseline.csv`, `pll_baseline.csv`. All available.

**Personal engagement:** Your original 0.25/0.75 figure was a sensible estimate. Testing it and discovering where it breaks, then explaining *why*, is exactly the "testing predictions" the rubric rewards. Also: predict which cases deviate before reading the table, then check it on a real cube.

**Reflection:** Is 0.25/0.75 "wrong"? Its effect on the mean is +0.11 moves out of 32. It's a fine approximation for a mean, and worse for tails (P(T ≤ 20) off by 15%). Whether an error matters depends on what you ask.

**Rubric:** Use of mathematics, Personal engagement, Reflection.

**Transition:** You now have seven distributions to combine, and this is where chaining convolutions felt inadequate.

---

## Section 4: Combining distributions: polynomials and the DFT

**Purpose and claim:** Develops the DFT from your own route (polynomial evaluation and systems of equations), and proves when it applies.

**Pages:** 3.5

**Key content, in this order:**
1. **Distribution as polynomial:** P(z) = Σ pₘ zᵐ. Adding independent variables means multiplying polynomials.
2. **Why values multiply:** E[z^(X+Y)] = E[z^X] E[z^Y] for independent X, Y.
3. **Independence, proved rather than assumed** (the bijection argument). This is your strongest piece of original mathematics.
   - Last pair vs last layer: for a fixed placement and alg, states map one-to-one, so uniform stays uniform.
   - OLL vs PLL: 216 patterns of 288 states; each pattern's alg maps its 288 states one-to-one onto the 288 oriented states.
   - Separate AUFs fail this, and the exact tests confirm both results. So group the AUF with its alg.
4. **Evaluate at roots of unity.** Why those points: they are the ones where the inverse is simple.
5. **Inverse by solving the system through cancellation,** including why dividing by N appears.
6. **The bridge** from your table method to the summation formula, in one short step. `[NEEDED FROM AUTHOR]` confidence with this.
7. **Hand example:** a = [0.25, 0.75, 0, 0], cubed.

| Point | A | A³ |
|---|---|---|
| 1 | 1 | 1 |
| −i | 0.25 − 0.75i | −0.40625 + 0.28125i |
| −1 | −0.5 | −0.125 |
| i | 0.25 + 0.75i | −0.40625 − 0.28125i |

   Inverse: [1/64, 9/64, 27/64, 27/64], which is Binomial(3, 0.75), a neat connection worth stating.
8. **Padding:** N ≥ 12 + 15 + 22 + 1 = 50, shown by failing. N = 40 gives distance 0.049 and a plausible-looking mean of 30.48; N = 32 reads as a mean of 12.90.
9. **Applying it:** DFT of the three blocks gives the baseline. Agreement with full enumeration: 3 × 10⁻¹⁵. Compare with flat AUFs: 0.0093.
10. **Why the error isn't fixed by correct AUF numbers:** with exact marginals, the mean is right, but the distance stays 0.0093. The dependence is the real error.

**Math shown vs. summarized:** Show items 1 to 8 in full. Direct convolution of all blocks only as a one-line check. Mention the FFT in one sentence at most.

**Visuals:** `fig3_wraparound` (x: index read as total moves, y: probability). `fig2_flat_auf_error` (left: both distributions; right: the error in percentage points). DFT table above.

**Data needed:** All available, `results/results.json`. `[NEEDED FROM AUTHOR]` your own hand working.

**Personal engagement:** Your whole route is personal: you came at this through polynomial evaluation rather than the textbook formula. Show that route rather than reproducing a standard derivation. The N = 40 failure is a good place to show a mistake you'd want to catch.

**Reflection:** Honest justification. At N = 50, direct convolution is entirely manageable, and the DFT is not needed for speed. Its value is structural: all distributions combined at once, each evaluation point turning convolution into multiplication. Say this plainly; that is the claim your motivation can support.

**Rubric:** Use of mathematics (sophistication, rigour), Personal engagement, Mathematical communication.

**Transition:** That solves your original problem. But it relied on independence, and the techniques you most care about are designed to link steps.

---

## Section 5: Techniques and the failure of independence

**Purpose and claim:** The analytic turning point. Claims a precise rule: convolution of blocks stays exact unless a decision or a merge uses information from across a step boundary.

**Pages:** 3.5

**Key content:**

**5.1 The dividing line.** Table of techniques with convolution error:

| Stays exact (~10⁻¹⁵) | Breaks |
|---|---|
| Blind insert swap | Two-step inserts (0.019) |
| Speed-gated angles | Full look-ahead (0.037) |
| Fewest-move angles, per step | Cancellations (0.0031) |
| | OLL angles, two-step (0.030) |
| | COLL, two-step (0.013) |
| | OLLCP, two-step (0.073) |

This corrects the project's starting assumption that alternate angles break independence. They change an alg, but they never look across the boundary.

**5.2 Choice rules: mixture over states.** P(T = m) = Σₛ P(s) · 1[min over options of total moves in s = m]. Why the mean is unharmed and the spread is not: choosing whichever is short makes steps negatively correlated, so Var(sum) falls below the sum of Vars. Show that explicitly for OLLCP: SD 4.60 by convolution, 3.90 true.

**5.3 Inserts.** Only 14 of 41 cases (1/3 of states) can swap. Blind swapping costs moves (+0.05, +0.08). With two-step: 0.29, 0.34, together 0.49. Show `R U R' U2 + U R U2 R'` merging.

**5.4 The greedy rule's weakness.** Saves 0.49, while full look-ahead saves 1.20. Suboptimal in 14.8% of states; worse than no choice in 4.9%. Use the worked example: 1 move saved through OLL, 7 lost at PLL.

**5.5 Cancellations: a formula, not a simulation.** Mean saving E[S] = Σ P(a) P(b) J(a, b). This works *because of* the Section 4 proof: 0.33114 + 0.18326 + 0.00174 = 0.51614, matching enumeration to 14 decimal places. So independence gives the mean of something that itself breaks independence. Dominant pattern: `R'` then `R`.

**5.6 Angles: faster is not fewer.** Alternates average +1.5 (OLL) and +1.9 (PLL) moves. Speed-gated PLL costs 0.39 moves. Shorter-only gives +0.55. With cancellations the picture flips: +0.46 interaction, since skipping an AUF unblocks cancellations (P(no AUF before OLL) goes 0.27 → 0.46; cancellation saving 0.52 → 0.98).

**5.7 COLL: conditioning brings the DFT back.** Three classes, convolve within each, mix. Error 4 × 10⁻¹⁵ against 0.026 for plain convolution. This is the moment the DFT returns.

**5.8 ZBLS + ZBLL.** ZBLL alone: 0.91. With ZBLS: 7.22. Why: conditioning on edges oriented has probability 1/8 unless you pay 1.50 moves to guarantee it.

**Math shown vs. summarized:** Show in full: the mixture formula, the variance argument, E[S], the COLL mixture. Summarize individual technique calculations in the table. Cut 5.8 to a paragraph if space is short.

**Visuals:** `fig9_independence_failure` (right mean, wrong tails) and `fig8_paired_saving` (x: moves saved, y: probability). Dividing-line table.

**Data needed:** All available, `RESEARCH_FINDINGS.md` sections 5 to 7.

**Personal engagement:**
- Your own sledge example, made concrete: the worked scramble, done on a real cube, saves 9 moves.
- Your own speed-ranking rule producing a result you may not have expected. Discussing that openly is strong engagement.
- Working out *why* the angle/cancellation interaction happens, then testing it, is exactly the rubric's "making and testing predictions".

**Reflection:** Limits of the two-step rule as a model of a human: it needs the OLL read before inserting, and full look-ahead needs the PLL read too, so 1.20 is a ceiling. Inserts here are limited to baseline algs ending in your forms, so 0.49 is conservative.

**Rubric:** All five, mainly Use of mathematics and Reflection.

**Transition:** With each technique modelled validly, they can all be put side by side.

---

## Section 6: Results and comparison

**Purpose and claim:** One compact comparison. Claims that value depends on what you count: moves, tails, or moves per alg.

**Pages:** 1.5

**Key content:** Summary table (subset of `RESEARCH_FINDINGS.md` section 6):

| Technique | Mean | SD | Saved | Fewer / more (paired) | P(T ≤ 25) | Net algs |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 32.43 | 4.83 | | | 0.082 | 0 |
| U2 + sledge, two-step | 31.94 | 4.76 | 0.49 | 13.0% / 4.9% | 0.095 | 0 |
| Cancellations | 31.91 | 4.97 | 0.52 | 24.0% / 0% | 0.105 | 0 |
| Angles, fewest moves | 31.88 | 4.73 | 0.55 | 20.9% / 3.7% | 0.091 | 48 |
| COLL, two-step | 31.97 | 5.04 | 0.46 | 7.9% / 0% | 0.103 | 40 |
| All CFOP, fewest moves | 30.37 | 4.92 | 2.06 | 51.4% / 5.6% | 0.161 | 88 |
| OLLCP, two-step | 28.50 | 3.91 | 3.93 | 67.5% / 0% | 0.190 | 329 |
| ZBLS + ZBLL | 25.21 | 3.51 | 7.22 | 87.8% / 8.1% | 0.517 | 733 |
| 1LLL | 23.73 | 3.03 | 8.70 | 95.6% / 1.5% | 0.719 | 3,837 |

Points to draw out:
- **Moves per 100 algs:** OLLCP 1.19, COLL 1.16, angles 1.14, ZB 0.99, 1LLL 0.23.
- **Full distributions matter:** ZB raises P(T ≤ 25) sixfold, from 8% to 52%, and narrows the spread too. A mean alone hides both.
- **The stage metrics you asked for:** LS, LS + OLL, OLL + PLL, and LSLL are all in the findings table.
- **Interactions are mostly additive,** with the angle/cancellation exception.

**Math shown vs. summarized:** Table only. Full table and factorial grid in the appendix.

**Visuals:** `fig4_overlay_pmf` (one overlay). `fig7_tradeoff` (x: algs added, log scale; y: moves saved). Summary table. `fig6_moves_saved` is an optional extra if space allows, and `fig5_overlay_cdf` can replace the overlay if you want to emphasise thresholds.

**Data needed:** All available.

**Personal engagement:** Which trade-off would you actually take as a solver? That's a real decision only you can make, and the numbers now support it. `[NEEDED FROM AUTHOR]`

**Reflection:** "Saved per alg" treats all algs as equally hard to learn and use. Is it fair to compare a technique you can adopt tomorrow (cancellations) with 733 algs?

**Rubric:** Organization (concise), Mathematical communication, Reflection.

**Suggested cuts if over length:** OLS, ZB with shortest ZBLS, the factorial grid (keep the one interaction), and `fig6`.

**Transition:** These numbers describe an idealized solver. Does that solver resemble a real one?

---

## Section 7: Validation, reflection, conclusion

**Purpose and claim:** Tests the model two ways. Claims the model is internally certain, close to pros on step lengths, wrong on which states pros face, and that comparisons between techniques are more trustworthy than absolute totals.

**Pages:** 2

**Key content:**

**Internal validation (brief):** independent brute-force solver: 0 mismatches over 20,000 states; magiccube replay: 1000 of 1000; DFT vs enumeration agree to 10⁻¹⁵. Mention the bug the checks caught (OLL 37 into T perm cut off by a short window). It shows the validation did real work. Also state the one known gap: if a whole OLL cancels away, the pairwise model misses further merging (1 state in 61,000, worth 0.00004 moves).

**Against pros:**

| | Model | Pros |
|---|---:|---:|
| Last pair | 7.81 | 7.48 ± 0.09 |
| OLL alg | 9.62 | 9.26 ± 0.03 |
| PLL alg | 13.06 | 13.14 ± 0.04 |
| Edges oriented after last pair | 12.5% | 20.1% |

Uniform states are rejected: F2L p ≈ 10⁻⁵, OLL and PLL far beyond that. Two causes that can't be separated: real skill (choosing pairs, edge control), and selection (reconstructed solves skew fast).

**Reflection:**
- **Movecount vs time, finally.** Pros pay 0.73 extra PLL moves compared with baseline algs at their own case frequencies. That's the most compelling evidence available that they are optimising something else, and your own speed-gated result says the same.
- Recognition is free in this model.
- What you'd do next: real-solve state distribution in place of uniform, and an estimate of recognition cost.
- Was the DFT worth learning for this? `[NEEDED FROM AUTHOR]`. Your honest view belongs here.

**Conclusion:** Answer the question directly: distributions compare techniques validly when each is combined in the way its dependence allows. Return to the DFT: it solved the problem you set it, its limits then told you more about the techniques than the DFT itself did.

**Math shown vs. summarized:** χ² and binomial results as stated values only; don't develop those tests. Keep them secondary so they don't distract from the argument.

**Visuals:** `fig10_last_pair_vs_real` (x: moves in last pair; y: probability). Model vs pros table.

**Data needed:** All available. `results/real_solves.json`, `results/validation.json`.

**Personal engagement:** Try the worked examples on your own cube; your own experience with whether a two-step decision is recognisable in real time is authentic evidence no dataset here gives. `[NEEDED FROM AUTHOR]`

**Reflection:** This whole section.

**Rubric:** Reflection, Organization, Personal engagement.

---

## Appendices

- **A.** Full F2L, OLL and PLL tables (`results/tables/*_baseline.csv`).
- **B.** Full technique table and factorial grid (`techniques.csv`, `factorial.csv`).
- **C.** Code structure and how to reproduce, briefly, noting technology was used for enumeration while all key mathematics is shown by hand in the body.

## Quality checks

| Check | Status |
|---|---|
| Clear central argument, not a list of techniques | Yes: the dividing line, and its consequences |
| Every technique has a valid mathematical model | Yes: all exact, section 5 of the findings gives the type |
| Each use of convolution justified by independence | Yes: proved in Section 4, tested in 5.1 |
| Reflection and engagement spread throughout | Yes: each section has its own |
| Repetitive material moved to tables or appendix | Yes: case data is in the appendix |
| Reader sees why full distributions matter | Yes: fig 9, and ZB's P(T ≤ 25) |
| Page budget at most 16 | Yes: 15.5 |
