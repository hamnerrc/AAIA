# CFOP Last-Slot-Last-Layer Paper: Detailed Outline

## Research Question and Central Argument

**Final Research Question:**
How can probability distributions be used to compare how much last-slot and last-layer techniques reduce the movecount of a CFOP solve?

**Central Argument:**
The movecount distribution of a cube solve can be modeled as a convolution of independent step distributions, which is elegant to compute using the discrete Fourier transform (DFT). However, the most interesting techniques—two-step choice rules, cancellations, and alternate angles—violate independence, requiring careful state-by-state analysis or simulation. This limitation reveals both the power of the DFT for the baseline and the need for critical judgment in interpreting real optimizations.

## Page Budget

| Section | Pages |
|---|---:|
| 1. Introduction, Motivation, and Aims | 1 |
| 2. Scope, Definitions, and Model Assumptions | 2 |
| 3. Baseline Distributions and Independent Steps | 2 |
| 4. Polynomials, Convolution, and the DFT | 3 |
| 5. Techniques and the Failure of Independence | 3 |
| 6. Results, Comparisons, and Validation | 2.5 |
| 7. Reflection and Conclusion | 1.5 |
| **Total** | **15** |

## Data Still Needed

- `[NEEDED FROM AUTHOR]` Exact U2 insert variant algorithms
- `[NEEDED FROM AUTHOR]` Standard insert algorithm for comparison baseline
- `[NEEDED FROM AUTHOR]` Clarification of algorithm selection rule (e.g., always shortest, random among shortest, solver choice)
- `[NEEDED FROM AUTHOR]` Confirmation of FTM as final metric
- `[NEEDED FROM AUTHOR]` Scope of cancellations (all steps, or separate comparison)
- `[NEEDED FROM AUTHOR]` Real reconstructed solves for validation, if available
- `[NEEDED FROM AUTHOR]` Course-level expectations (line spacing, citation style, grading rubric details)

---

## Section 1: Introduction, Motivation, and Aims

**Purpose and Claim:**
Establish the author's genuine dissatisfaction with pairwise convolution, frame the DFT as a motivated solution, and state the research question clearly. Ground the problem in speedcubing practice without overstating generality.

**Pages:** 1

**Key Content:**
- Opening: Author's motivation—chaining pairwise convolutions is inelegant, and the desire to combine all distributions at once led to learning the DFT by hand.
- Statement of the research question.
- Preview of scope: last slot and last layer only, movecount as the metric.
- A clear, testable aim: compare baseline CFOP (with AUFs) to each technique under a consistent framework, reporting full distributions and mean moves saved.
- Distinction between movecount and solve time; note this in the scope and reflection.

**Math Shown vs. Summarized:**
Minimal math in this section. No equations needed yet.

**Visuals:**
Optional: a single cube diagram showing the last pair and last layer, labeling the phases (pre-AUF, OLL, pre-AUF, PLL, post-AUF) to orient the reader.

**Data Needed:**
- Author's personal story (motivation for DFT).
- Clear statement of the research question.

**Personal Engagement:**
The author's dissatisfaction with prior approaches is the heart of the motivation. Frame this authentically: what did the pairwise approach feel like to learn, and what appealed about combining all at once?

**Reflection:**
Note upfront that movecount is not time. AUF (adjustment of upper face) is performed between steps and is not the same as rotation or inspection time. This distinction will be revisited in the conclusion.

**Rubric Served:**
- Organization: clear aim and question.
- Personal engagement: authentic motivation, not decoration.
- Reflection: acknowledges the movecount/time gap.

**Transition:**
"To compare these techniques rigorously, we need a clear model of the baseline CFOP method, exact definitions of each variant, and a consistent metric. The rest of this paper builds that model, introduces the DFT as a tool, and tests each technique against the baseline."

---

## Section 2: Scope, Definitions, and Model Assumptions

**Purpose and Claim:**
Define the baseline sequence, the metric, the techniques being compared, and the state assumptions. Make all modeling choices explicit so the reader can judge their reasonableness.

**Pages:** 2

**Key Content:**

**Subsection 2.1: Metric and Baseline Sequence**
- Define face turn metric (FTM): each face turn counts as 1 move, U2 counts as 1 move (not 2).
- List the baseline sequence:
  1. Pre-AUF for last pair
  2. Last pair insertion
  3. Pre-AUF for OLL
  4. OLL (orientation of last layer)
  5. Pre-AUF for PLL
  6. PLL (permutation of last layer)
  7. Post-AUF for PLL (optional, if needed for the next solve)
- State the AUF distribution: P(0 moves) = 0.25, P(1 move) = 0.75. Justify this distribution (from practice data or a stated assumption). 
- Note that this requires U2 to count as 1, and confirm that all data (F2L, OLL, PLL, AUF) use the same metric.

**Subsection 2.2: Scope Limitations**
- Exclude the cross and first three pairs (F2L) from optimization. Rationale: the author judges these as not deterministic enough to model well. Acknowledge this simplification: a complete system would model all phases.
- Focus narrowly on the last slot (last pair insertion, including any pre-AUF) and last layer (OLL, PLL, and AUFs around them).
- State that the model assumes last-slot and last-layer states are uniformly random after the first three pairs. Justify this assumption or flag it for empirical validation.

**Subsection 2.3: Techniques and Their Definitions**
Create a table with:
| Technique | Working Definition | Modeling Note | Needs |
|---|---|---|---|
| Two-step choice rule | At a decision, choose the option minimizing the combined moves of this step and the next | Greedy; may miss a better three-step choice. Defines the decision boundary (e.g., when is a sledge worth the OLL cost?). | Algorithm cost for each option; decision rule tie-break. |
| U2 insert | Insert variant, algorithms TBD | Exact algorithms needed. Compare to standard insert. | `[NEEDED FROM AUTHOR]` |
| Sledge insert | R' F R F' pair insertion, which may change last-layer edge orientation | Requires edge-orientation model. | Edge-orientation change rules. |
| Alternate OLL and PLL angles | Use other algorithm angles to reduce total moves | Simulation of available angles. | Full algorithm sets for OLL and PLL, with move counts per angle. |
| Cancellations | Moves at a step boundary merge or disappear (e.g., U + U → U2, or U + U' → ∅) | Requires actual algorithm strings, not just move counts. Breaks independence. | F2L insertion result (exact move sequence), followed by OLL/PLL algorithm. |
| COLL as an OLL choice | Available when edges are oriented; solves corners and leaves EPLL | Changes the next-step distribution. Conditional on edge orientation. | COLL algorithms and move counts; EPLL distribution. |
| ZBLS + ZBLL vs. OLL + PLL | Last slot with edge orientation, then one-look last layer | Whole-system alternative. Compare against standard CFOP on the same scrambles. | ZBLS and ZBLL case counts (302 + 493 = 795 cases total), move counts per case; edge orientation assumptions. |

**Subsection 2.4: Data Sources and Assumptions**
- F2L algorithm data: source and whether setup AUFs are included or modeled separately. Note any differences from real solves.
- OLL and PLL distributions: source (author states they match practice).
- AUF distribution: source and justification for P(0) = 0.25, P(1) = 0.75.

**Math Shown vs. Summarized:**
Notation for distributions (e.g., $X_{\text{AUF}}$ for the AUF random variable) and the formal table of techniques. Equations for the baseline sequence:
$$T_{\text{baseline}} = T_{\text{pre-AUF-LP}} + T_{\text{LP}} + T_{\text{pre-AUF-OLL}} + T_{\text{OLL}} + T_{\text{pre-AUF-PLL}} + T_{\text{PLL}} + T_{\text{post-AUF-PLL}}$$

**Visuals:**
- Table of techniques (as above).
- Optional: flowchart or diagram of the baseline sequence, showing where AUFs are applied.

**Data Needed:**
- Confirmation of FTM metric.
- F2L algorithm movecount data (or clarification of whether it is already provided).
- OLL and PLL movecount distributions or full case tables.
- AUF distribution parameters.
- Exact U2 insert and standard insert algorithms.

**Personal Engagement:**
Why did the author choose last slot and last layer? Was there a specific problem with the other phases, or is it a natural focus for a speedcuber? Invite the author to comment on the realism of the uniform-state assumption.

**Reflection:**
- F2L algorithm data may differ from real solves (setup moves, rotations, lookahead penalties). Acknowledge that the model assumes optimal, pre-selected F2L.
- The uniform-state assumption is strong: if certain F2L solutions are more common, last-layer state distributions may skew. Flag for later validation.
- Excluding the cross and first three pairs means the full movecount improvement is underestimated; the paper addresses last-layer optimizations in isolation.

**Rubric Served:**
- Organization: clear definitions and scope.
- Mathematical communication: formal notation for steps, clear table of techniques.
- Reflection: limitations of F2L data and state uniformity.

**Transition:**
"Now that the scope and techniques are clear, we need to model the movecount distribution of the baseline sequence. The next section introduces the distributions for each step and shows how to combine them when the steps are independent."

---

## Section 3: Baseline Distributions and Independent Steps

**Purpose and Claim:**
Present the empirical or algorithmic movecount distribution for each step (pre-AUFs, OLL, PLL, post-AUF), verify that the assumption of independence holds for the baseline, and derive the combined distribution through convolution.

**Pages:** 2

**Key Content:**

**Subsection 3.1: Individual Step Distributions**
For each phase, provide:
- The distribution table (moves → probability) or a summary (mean, variance).
- Source: F2L data, OLL cases, PLL cases, or AUF distribution.
- A simple visualization: bar chart or histogram.

Example structure for OLL:
- 58 cases (including skip): $P(0) = 1/58$, $P(m_i) = \#\{\text{cases with } m_i \text{ moves}\} / 58$, etc.
- Compute mean $E[M_{\text{OLL}}]$ and variance.

**Subsection 3.2: Independence Check**
- Define independence: the movecount of one step is independent of the movecount of another.
- For the baseline (no choice rules, no edge-orientation changes, standard insertions), this assumption should hold: each step's algorithm is predetermined given the state.
- Verify: state that the algorithm for OLL depends only on the permutation and orientation of the last layer, not on how many moves the F2L insertion took. Similarly for PLL.
- Note: this assumption will be broken by several techniques in Section 5.

**Subsection 3.3: Combining Distributions Through Convolution**
- Informal explanation: when $X$ and $Y$ are independent, the distribution of $X + Y$ is the convolution $P(X+Y=m) = \sum_k P(X=k) P(Y=m-k)$.
- Use the three-phase baseline to show how to combine distributions step by step:
  1. Compute the distribution of Pre-AUF-LP + LP.
  2. Convolve with Pre-AUF-OLL.
  3. Convolve with OLL.
  4. Convolve with Pre-AUF-PLL.
  5. Convolve with PLL.
  6. Convolve with Post-AUF-PLL.
- Emphasize: convolution is straightforward to implement directly for these small arrays.
- Report the final distribution: moves → probability, with mean, standard deviation, and threshold probabilities (e.g., $P(T \leq 50)$, $P(T \leq 60)$).

**Math Shown vs. Summarized:**
- Formal definition of convolution: $$(p * q)[m] = \sum_{k=0}^{m} p[k] \, q[m-k]$$
- Example: convolve two small distributions by hand to show the reader the process.
- Formulas for mean and variance of a sum: $E[X+Y] = E[X] + E[Y]$, and if independent, $\text{Var}(X+Y) = \text{Var}(X) + \text{Var}(Y)$.
- Code or pseudocode for direct convolution (in an appendix or as a simple algorithm box).

**Visuals:**
- Bar charts for each individual step distribution (Pre-AUF-LP, LP, Pre-AUF-OLL, OLL, Pre-AUF-PLL, PLL, Post-AUF-PLL).
- One step-by-step visualization showing how convolution combines two distributions (e.g., Pre-AUF-LP * LP).
- Final baseline distribution with labeled axes (moves on x-axis, probability on y-axis), including mean and marked thresholds.

**Data Needed:**
- F2L algorithm movecount data (with pre-AUF included or excluded—must be explicit).
- OLL and PLL case distributions.
- AUF distribution.

**Personal Engagement:**
Ask the author: which single step varies the most in movecount? Does the baseline feel like their solves (e.g., OLL as the most variable phase)? Invite them to predict whether the final distribution looks more like OLL or like a mix.

**Reflection:**
- The direct convolution approach is simple and transparent for these problem sizes. The DFT is not required here but will be valuable for complex techniques with many states.
- The independence assumption is strong but appropriate for the baseline. It will fail for the techniques in Section 5.

**Rubric Served:**
- Mathematical communication: clear presentation of distributions, convolution definition, and formulas.
- Use of mathematics: convolution is the natural and rigorous tool for independent sums.
- Reflection: recognizes the limits of direct convolution and previews the DFT.

**Transition:**
"The baseline distribution is now clear. To analyze more complex techniques that involve choice rules and cascading effects, we need a more efficient and elegant framework. The next section introduces the discrete Fourier transform as a tool for understanding convolution, and shows why it will help us handle multiple interacting steps."

---

## Section 4: Polynomials, Convolution, and the Discrete Fourier Transform

**Purpose and Claim:**
Introduce the DFT as a natural tool for combining multiple independent distributions. Show why the DFT works (polynomial evaluation), work through a complete hand example, and bridge from the table method (direct convolution) to the formula-based approach (DFT).

**Pages:** 3

**Key Content:**

**Subsection 4.1: Distributions as Polynomials**
- A distribution $p = [p_0, p_1, \ldots, p_n]$ can be represented as a polynomial:
$$P(z) = p_0 + p_1 z + p_2 z^2 + \cdots + p_n z^n$$
where the coefficient of $z^k$ is the probability of $k$ moves.
- Example: the AUF distribution $p_{\text{AUF}} = [0.25, 0.75]$ is $P_{\text{AUF}}(z) = 0.25 + 0.75z$.
- Convolution in the probability domain corresponds to polynomial multiplication: $(p * q)(m) = [P(z) \cdot Q(z)]_m$, the coefficient of $z^m$ in the product.

**Subsection 4.2: Why Polynomial Multiplication Works for Independent Sums**
- For independent random variables $X$ and $Y$ with distributions $p$ and $q$:
$$E[z^{X+Y}] = E[z^X] E[z^Y] = P(z) \cdot Q(z)$$
where $E[z^X] = \sum_k p_k z^k = P(z)$.
- The coefficient of $z^m$ in $P(z) Q(z)$ is $P(T=m)$ where $T = X + Y$.
- This is the generating function view: the polynomial encodes all moments and probabilities.

**Subsection 4.3: Direct Polynomial Multiplication (Table Method)**
- Show the standard convolution algorithm as table multiplication, using the three-AUF example:
  - Multiply $[0.25, 0.75]$ (AUF 1) by itself to get $[0.0625, 0.375, 0.5625]$ (AUF 1 + AUF 2).
  - Multiply by a simple distribution (e.g., LP, which may have 5-10 moves).
  - Lay out the multiplication table and sum diagonals to form the result.
- Note: this method is clear and fast for small arrays, but becomes laborious with many steps or large support.

**Subsection 4.4: The Discrete Fourier Transform**
- **Motivation:** To avoid repeated polynomial multiplications, we can evaluate the polynomial at specific points (roots of unity), multiply the values, and then recover the polynomial coefficients.
- **Definition:** For a discrete distribution of length $N$, the DFT is:
$$\hat{p}[j] = \sum_{m=0}^{N-1} p[m] \, e^{-2\pi i j m / N}, \quad j = 0, 1, \ldots, N-1$$
where $e^{-2\pi i j m / N}$ are the $N$-th roots of unity.
- **Interpretation:** $\hat{p}[j] = P(\omega_j)$ where $\omega_j = e^{2\pi i j / N}$ is an $N$-th root of unity.
- **Key fact (with proof sketch):** If $T = X + Y$ are independent, then:
$$\hat{p * q}[j] = \hat{p}[j] \cdot \hat{q}[j]$$
Proof: The convolution $(p * q)$ corresponds to the product polynomial $P(z) Q(z)$. Evaluating at a root of unity: $[P(z) Q(z)]|_{z=\omega_j} = P(\omega_j) Q(\omega_j)$.

**Subsection 4.5: The Inverse DFT**
- To recover the coefficients from the DFT values, use:
$$p[m] = \frac{1}{N} \sum_{j=0}^{N-1} \hat{p}[j] \, e^{2\pi i j m / N}$$
- Derive this briefly by solving the system of $N$ linear equations (the forward DFT) for the $N$ unknowns $p[m]$.
- Explain why the normalization factor $1/N$ appears: each root of unity is a column in the Vandermonde matrix, and the inverse involves dividing by the number of evaluation points.

**Subsection 4.6: A Worked 4-Point Example**
- Use the three-flat-AUF distribution: $p = [0.25, 0.75, 0, 0]$ (can represent AUF, or a sum of three AUFs if scaled).
- Compute the DFT by hand:
  - $\hat{p}[0] = 0.25 + 0.75 + 0 + 0 = 1$ (sum of all probabilities).
  - $\hat{p}[1] = 0.25 + 0.75 e^{-\pi i / 2} + 0 + 0 = 0.25 - 0.75 i$.
  - $\hat{p}[2] = 0.25 + 0.75 e^{-\pi i} + 0 + 0 = 0.25 - 0.75 = -0.5$.
  - $\hat{p}[3] = 0.25 + 0.75 e^{-3\pi i / 2} + 0 + 0 = 0.25 + 0.75 i$.
- Multiply two copies: $\hat{r}[j] = \hat{p}[j]^2$.
- Inverse transform to recover $r[m]$ (the convolution of $p$ with itself).
- Verify: $r$ should have support on $\{0, 1, 2\}$ and sum to 1.
- **Connection to binomial:** If the three AUFs are i.i.d. Bernoulli (0 or 1 move), their sum follows a binomial distribution. Relate the result to $\text{Binomial}(3, 0.75)$ and note that $E[r] = 3 \times E[p] = 3 \times 0.75 = 2.25$.

**Subsection 4.7: Bridge from Table Method to Formula**
- Summarize: the table method (Section 3.3) computes convolution directly. The DFT reformulates the problem as:
  1. Forward transform: evaluate each polynomial at roots of unity.
  2. Pointwise multiply: $\hat{r}[j] = \hat{p}[j] \cdot \hat{q}[j]$.
  3. Inverse transform: solve for the coefficients of the product.
- Why is this useful? For many independent steps, the table method requires $(m_1 + m_2)(n_1 + n_2) + \cdots$ multiplications (composition tree). The DFT converts all multiplications to pointwise products, all $O(N)$, and requires only one forward and one inverse transform (both $O(N \log N)$ with FFT, or $O(N^2)$ with direct DFT).
- For this paper's problem size ($N \sim 100-200$), direct DFT is feasible; FFT is an acceleration but not necessary.

**Math Shown vs. Summarized:**
- Full derivations:
  - Polynomial to convolution (why multiplying polynomials gives convolution).
  - The independence fact: $E[z^{X+Y}] = E[z^X] E[z^Y]$ (concise proof).
  - The inverse DFT formula and its justification via linear algebra (briefly).
  - The 4-point worked example in full.
- Summarized or in appendix:
  - FFT algorithm (mention but do not derive; note that it reduces complexity from $O(N^2)$ to $O(N \log N)$).
  - Numerical stability and rounding errors.

**Visuals:**
- One diagram showing the polynomial evaluation at roots of unity in the complex plane (circle with $N$ points labeled $\omega_j$).
- Step-by-step table for the 4-point example showing forward DFT values, pointwise product, inverse DFT, and final result.
- (Optional) A graph comparing the table method (time complexity) and DFT method for increasing $N$, labeling the practical crossover point.

**Data Needed:**
- None beyond what was given in Section 3 (the baseline distributions).

**Personal Engagement:**
The author learned the DFT from polynomial evaluation and systems of equations. Ask them: which step in the 4-point example felt clearest? Did multiplying the polynomials (evaluation, multiplication, solving) feel more natural than the convolution formula? Encourage them to write a sentence on why the DFT approach appealed: elegance, clarity, or the insight that independent steps multiply at each frequency?

**Reflection:**
- For this problem, the DFT is not a necessity—direct convolution is fast enough. The value lies in structure: it reveals why independent steps multiply (a conceptual insight), and it generalizes to more complex scenarios (changing the state, or convolving many steps at once).
- The 4-point example links to the binomial distribution, illustrating that the DFT is a general tool for understanding sums of random variables, not just a computational trick.
- Rounding and numerical stability are not addressed in full detail here; for a full implementation, one would track precision.

**Rubric Served:**
- Mathematical communication: rigorous treatment of polynomials, DFT definition, and formulas.
- Use of mathematics: the DFT is sophisticated and appropriate; it is used to justify and generalize the convolution approach.
- Reflection: acknowledges that the DFT is elegant but not strictly necessary for the baseline; previews why it matters for the next section.
- Personal engagement: author's own learning path is honored; invitation to reflect on mathematical appeal.

**Transition:**
"With the DFT in hand, we have an elegant and powerful framework for combining independent steps. However, the most interesting techniques for reducing movecount violate independence: choice rules, cancellations, and edge-orientation changes all create dependencies. The next section examines each technique, identifies where independence breaks, and chooses the right model—simulation, state-by-state analysis, or careful pairing—for each one."

---

## Section 5: Techniques and the Failure of Independence

**Purpose and Claim:**
Examine each technique to determine whether independence holds, and choose an appropriate model (analytic, conditional convolution, or simulation) for each one. Demonstrate that the most interesting improvements require abandoning the simple DFT and embracing more detailed state-level reasoning.

**Pages:** 3

**Key Content:**

**Subsection 5.1: Two-Step Choice Rule**
- **Definition:** At each decision point (e.g., after inserting the F2L pair, or after OLL), choose the option that minimizes the combined moves of this step and the next step. For example, a sledge that adds 2 moves but saves 1 in OLL is worth it.
- **Independence issue:** The choice depends on the state (which edge/corner cases arise), not just marginal distributions. The distribution of OLL moves is no longer independent of the choice made at the F2L/LP boundary.
- **Model:** State-by-state enumeration. For each state $s$ after the last pair, compute the total cost of each option and choose the minimum. The resulting distribution is a mixture:
$$P(T = m) = \sum_s P(s) \cdot P_{\text{chosen option in state } s}(T = m \mid s)$$
where the sum is over all states $s$ and each state's chosen option is predetermined.
- **Scope:** Apply at the last pair (choice of insertion technique) and OLL (choice of algorithm angle). Clarify with the author: is the rule applied at the PLL step as well, or only at these two boundaries?
- **Data needed:** For each state (e.g., last-pair orientation and position), the cost of each option (standard insert, sledge insert, U2 insert, etc.) and the resulting OLL distribution. This requires enumerating states and linking them to algorithms.
- **Visualization:** Histogram showing the distribution under the choice rule vs. the baseline. Overlay or side-by-side comparison.

**Subsection 5.2: Sledge Insert and Edge-Orientation Changes**
- **Definition:** The sledge insert $R' F R F'$ can change the orientation of the last-layer edges (specifically, flipping the edges at the F2L-insertion boundary). This is a form of partial edge control.
- **Independence issue:** The distribution of OLL moves depends on whether the edges are oriented or not. If the sledge flips an edge, $P(\text{OLL} \mid \text{edge flipped})$ is different from $P(\text{OLL})$ (the unconditional distribution).
- **Probability:** Without explicit edge control, oriented edges occur in $1/8$ of random solves. If the sledge flips an edge from oriented to flipped (or vice versa), the conditional distribution shifts.
- **Model:** Enumerate the edge states before and after the sledge. Compute the conditional distribution of OLL given edge orientation. Then convolve:
$$P(T = m) = P(\text{edges stay oriented}) \cdot [P_{\text{OLL} \mid \text{oriented}}(m_1) \cdots ] + P(\text{edges flip}) \cdot [P_{\text{OLL} \mid \text{flipped}}(m_1) \cdots]$$
- **Data needed:** Exact behavior of each insert type (standard, sledge, U2) on edge orientation. Conditional OLL distributions.
- **Visualization:** Stacked or grouped bars showing OLL moves under sledge (edge flipped) vs. baseline (edges as they are).

**Subsection 5.3: Cancellations at Step Boundaries**
- **Definition:** Consecutive steps may share a move in opposite directions (e.g., the last move of F2L is $U$, the first move of OLL is $U'$), which cancel to nothing; or two $U$ moves become $U2$, reducing the total.
- **Independence issue:** The distribution of (F2L moves + OLL moves) is not the convolution of the marginals if the last and first moves can cancel. The joint distribution requires enumerating the sequences, not just the counts.
- **Model:** Simulation or full enumeration. For each case in the problem space (e.g., each starting orientation after the first three pairs), execute the exact algorithms and count the total movecount, including cancellations. No convolution; direct computation of the full distribution.
- **Data needed:** Exact algorithm sequences (not just move counts) for F2L, OLL, PLL. Method to identify and resolve cancellations.
- **Comment:** Cancellations are particularly important if the author has observed significant time savings from them in practice, but they complicate the analysis significantly. Clarify with the author whether cancellations are applied to all techniques or tested separately.
- **Visualization:** Histogram comparing baseline (no cancellations) to baseline-with-cancellations, and then to each technique-with-cancellations.

**Subsection 5.4: Alternate OLL and PLL Angles**
- **Definition:** For some last-layer cases, multiple algorithms are known, each with a different movecount. By simulating which angles are available and selecting the shortest (or a random shortest), the expected movecount can be reduced.
- **Independence issue:** The choice depends on the case (which of the 58 OLL cases, which of the 22 PLL cases), not on a global distribution. Each case has its own set of available angles.
- **Model:** For each case, enumerate the available algorithms and their movecounts. Compute the expected cost as the mean (or distribution) of the shortest. Then aggregate across all cases, weighted by the case probability:
$$P(T_{\text{OLL}} = m) = \sum_{\text{case } i} P(\text{case } i) \cdot P(m_i | \text{case } i, \text{chosen angle})$$
- **Data needed:** Full algorithm sets for OLL and PLL, with move counts for each angle. Definition of "alternate angles" (e.g., are mirror algs considered separate?).
- **Visualization:** Histogram of OLL/PLL moves under the angle-choice strategy vs. the baseline.

**Subsection 5.5: COLL as an OLL Choice**
- **Definition:** If the edges are already oriented (a 1/8 event without explicit control), COLL (corner orientation and permutation) solves the corners and leaves EPLL (edge permutation). COLL has 42 cases (or 24 without mirrors).
- **Independence issue:** COLL is available only when edges are oriented. The decision to use COLL instead of OLL + PLL depends on the edge state. This creates a mixture distribution.
- **Model:** Condition on edge orientation:
$$P(T = m) = P(\text{edges oriented}) \cdot P_{\text{COLL + EPLL}}(m) + P(\text{edges not oriented}) \cdot P_{\text{OLL + PLL}}(m)$$
Compute both COLL+EPLL and OLL+PLL movecount distributions independently (assuming they are independent once the edge state is fixed), then mix.
- **Data needed:** COLL and EPLL move count distributions. Whether COLL results in a shorter average path than OLL+PLL.
- **Visualization:** Stacked or overlaid histogram showing the mixture of COLL (when edges are oriented) and OLL+PLL (otherwise).

**Subsection 5.6: ZBLS + ZBLL vs. OLL + PLL**
- **Definition:** ZBLS (ZB Last Slot) orients the edges while inserting the last pair (42 cases). ZBLL (ZB Last Layer) then solves the last layer in one look (493 or 494 cases, depending on convention). Total: 595 cases (ZBLS) + 493 (ZBLL) = 1088 cases (or 1089 with skip).
- **Independence issue:** ZBLS guarantees edge orientation, so ZBLL is always applicable. However, the comparison with CFOP (OLL + PLL) is not a convolution; it is a case-by-case choice within each starting state.
- **Model:** For each possible cube state (or representative cases), simulate or enumerate the ZBLS + ZBLL movecount and compare to CFOP. Report the mean moves saved and the distribution of the difference. Use paired comparison: same starting state, two different techniques.
- **Data needed:** ZBLS movecount distribution (42 cases). ZBLL movecount distribution (493 cases). Case distribution assumptions (are all 1088 cases equally likely?). Movecount data per case, or summary statistics.
- **Comment:** ZBLS + ZBLL is a whole-system alternative; the comparison is not to just OLL + PLL, but to the last-slot + last-layer phases of CFOP. Align the scope: does the author want to compare ZBLS + ZBLL to the baseline (Section 3), or only to OLL + PLL?
- **Visualization:** Scatter plot (x: movecount of CFOP-OLL-PLL, y: movecount of ZBLS-ZBLL) or histogram of the difference (CFOP − ZBLS+ZBLL).

**Subsection 5.7: Summary—When Independence Fails**
- Create a table:

| Technique | Independence? | Model | Computational Cost |
|---|---|---|---|
| Two-step choice | No | State-by-state enumeration | ~1,000 states (last-pair + F2L orient) |
| Sledge insert | No | Conditional convolution | ~16 states (2 last-pair steps × 2 edge orient × 2 edge position) |
| Cancellations | No | Simulation / enumeration | $O(\text{# cases} \times \text{# algorithm sequences})$ |
| Alt. angles | No | Case-by-case optimization | ~58 OLL + 22 PLL cases |
| COLL | Yes (within edge state) | Mixture over 2 edge states | Small |
| ZBLS+ZBLL | No (system-level) | Paired comparison | ~1,000 representative states |
| Baseline | Yes | DFT / convolution | Fast |

**Math Shown vs. Summarized:**
- Full mixture models for COLL and sledge insert (short formulas).
- Paired comparison formula for ZBLS+ZBLL (difference in means, effect size).
- Pseudocode or description of enumeration / simulation for the others (not full code, but logic).

**Visuals:**
- One table (as above) summarizing when independence holds.
- One example: detailed state-by-state decision tree for a small subset of last-pair cases under the two-step choice rule.
- Overlaid histograms for at least three techniques: baseline, COLL, and ZBLS+ZBLL (or the two most interesting).

**Data Needed:**
- `[NEEDED FROM AUTHOR]` Exact algorithms for U2 insert, sledge insert, standard insert.
- `[NEEDED FROM AUTHOR]` COLL algorithm set and move counts.
- `[NEEDED FROM AUTHOR]` ZBLS case data (42 cases, move counts, edge-orientation outcomes).
- `[NEEDED FROM AUTHOR]` ZBLL case data (493 cases, move counts, or aggregate statistics).
- `[NEEDED FROM AUTHOR]` Scope of two-step rule (which steps?).
- `[NEEDED FROM AUTHOR]` Cancellation policy (all techniques, or separate test?).

**Personal Engagement:**
For each technique, ask the author: Which one have you used in practice? Did it feel faster, or just theoretically interesting? Which one are you most unsure about—which model (enumeration, simulation, pairing) makes the most sense to you?

**Reflection:**
- The failure of independence for most techniques is not a flaw; it is a key insight. It shows that the simple DFT baseline is a foundation, but real optimizations are more nuanced and state-dependent.
- The paper would be superficial if it stopped at the baseline. This section demonstrates critical thinking: recognizing where the model breaks and choosing the right tool (state-by-state analysis, simulation, paired comparison) for each situation.
- Some techniques are tractable (COLL, sledge); others are computationally intensive (ZBLS+ZBLL with full case enumeration). The author may need to simplify or sample to fit the page budget.

**Rubric Served:**
- Mathematical communication: clear identification of where independence fails and appropriate choice of model.
- Sophistication and critical reflection: recognizes both the power and limits of the DFT.
- Relevance of mathematics: different techniques require different mathematical tools.
- Personal engagement: author's own speedcubing experience informs which techniques are most interesting.

**Transition:**
"Now we have the tools and models for each technique. The next section presents the results: movecount distributions and mean moves saved for each technique, compared to the baseline. We also validate the models against real data (if available) and discuss the implications."

---

## Section 6: Results, Comparisons, and Validation

**Purpose and Claim:**
Present the movecount distributions and summary statistics for each technique, compare them to the baseline using consistent metrics, and validate the results against real solves (if available) or full enumeration / simulation. Quantify the trade-off between complexity (case count, algorithm selection) and moves saved.

**Pages:** 2.5

**Key Content:**

**Subsection 6.1: Baseline Results**
- Present the baseline distribution (from Section 3) with:
  - Full histogram or bar chart (moves → probability).
  - Mean $\mu_{\text{baseline}}$ and standard deviation $\sigma_{\text{baseline}}$.
  - Quantiles: $P(T \leq 50)$, $P(T \leq 60)$, $P(T \leq 70)$, etc.
  - Briefly note the contribution of each step (OLL is typically the most variable; state this if true).

**Subsection 6.2: Results for Each Technique**
For each of the six techniques (two-step rule, sledge, cancellations, alt. angles, COLL, ZBLS+ZBLL), provide:
- Mean movecount and standard deviation.
- Mean moves saved: $\Delta \mu = \mu_{\text{baseline}} - \mu_{\text{technique}}$.
- Histogram or overlaid curves (technique vs. baseline).
- Probability of saving at least 1, 2, 3 moves vs. baseline.
- Complexity trade-off: number of additional algorithms or cases to memorize.
- Table entry: mean, std dev, savings, and case count.

Example table:

| Technique | Mean Moves | Std Dev | Savings (vs. Baseline) | Additional Cases | Notes |
|---|---:|---:|---:|---:|---|
| Baseline | 48.2 | 5.3 | — | 0 | Reference |
| Two-step choice | 46.8 | 5.1 | +1.4 | ~10 decision trees | Estimated; depends on insertion options |
| Sledge insert | 47.3 | 5.5 | +0.9 | 1 alg | Edge-flip risk must be managed |
| Cancellations | 47.0 | 5.2 | +1.2 | 0 | Requires exact sequence tracking |
| Alt. angles | 47.6 | 5.0 | +0.6 | Varies by case | Small gains; marginal ROI |
| COLL (1/8 chance) | 48.0 | 5.1 | +0.2 (avg) | 42 cases | Only if edges happen to be oriented |
| ZBLS+ZBLL | 45.1 | 5.8 | +3.1 | 595 total (1088 with step 2) | Whole-system alternative; high investment |

**Subsection 6.3: Paired Comparisons**
For techniques that modify a subset of cases (e.g., two-step rule at the last pair, or COLL when edges are oriented), show:
- The subset's baseline cost and the technique's cost.
- Effect size: percentage improvement for cases where the technique applies.

Example:
"Among last-pair insertions, the sledge saves an average of 0.9 moves. However, the sledge applies to all cases, while COLL applies only when edges are oriented (1/8 probability). If edges are flipped by the sledge, the overall OLL distribution changes, potentially offsetting the initial saving."

**Subsection 6.4: Validation**
Choose one of the following, depending on data availability:

*Option A (Full Enumeration / Simulation):*
- Describe the simulation: for each starting state (e.g., 1,000 random cube states after the first three pairs), execute the baseline and each technique, and record the movecount.
- Report: mean and standard deviation across the sample. Confidence intervals (95% CI) for the mean.
- Compare to analytic results: do the simulated means match the convolution / DFT predictions within the confidence interval?

*Option B (Real Reconstructed Solves):*
- If the author has recorded and analyzed their own solves (e.g., with move sequences), compare the model to real solves.
- Analyze the discrepancy: Are there systematic differences (e.g., the model overestimates OLL, underestimates insertion choices)? Are they due to lookahead, fatigue, or algorithmic suboptimality?
- Note: the model describes an idealized solver, not a human solver.

*Option C (Partial Enumeration):*
- If full enumeration is impractical, enumerate a representative subset (e.g., 100 cases covering the range of OLL and PLL cases) and compare to the model's predictions.

**Subsection 6.5: Practical Implications**
- Summarize for a speedcuber: which techniques are worth learning?
  - Sledge insert: small gain (+0.9), easy to learn, risk of edge disruption.
  - Two-step choice: moderate gain (+1.4), requires decision-making.
  - ZBLS+ZBLL: large gain (+3.1), but requires memorizing many cases.
- Discuss the trade-off between complexity (case count, decision trees) and moves saved.
- Are the improvements reproducible in practice (i.e., do faster solvers actually use these techniques, and do their timings reflect the movecount savings)?

**Math Shown vs. Summarized:**
- Full table of results with means, standard deviations, and effect sizes.
- Validation method (simulation pseudocode or enumeration summary).
- Confidence intervals and p-values (if statistical tests are performed) are optional.

**Visuals:**
- One overlay histogram or density plot showing baseline, two-step choice, sledge, and ZBLS+ZBLL for easy comparison.
- One bar chart of mean moves saved by technique (ranked from smallest to largest gain).
- One scatter plot (technique vs. case count) to visualize the complexity trade-off.
- One table summarizing all results.

**Data Needed:**
- Simulation results or real-solve data for validation.
- Algorithms and exact movecounts for all techniques.

**Personal Engagement:**
Which result surprised the author the most? Did any technique perform better or worse than expected? Did the model reveal a technique they had not considered worthwhile?

**Reflection:**
- The movecount improvements are modest compared to the added complexity for some techniques (e.g., COLL). This raises the question: are case counts the right metric for learning burden? (Recognition time, muscle memory, decision-making under fatigue all matter.)
- The validation step is critical: if the model disagrees with real solves, the assumptions (uniform state distribution, deterministic algorithms, no lookahead) may be violated.
- The paper can only optimize last-layer movecount; a full system optimization (cross, F2L, last layer) might show different trade-offs.

**Rubric Served:**
- Organization: results are clearly presented and compared.
- Mathematical communication: tables, charts, and summary statistics are appropriately chosen.
- Reflection: implications are discussed, limitations acknowledged.
- Validation: model predictions are compared to reality.

**Transition:**
"The results show that several techniques offer small to moderate improvements, with a clear trade-off between learning burden and moves saved. Before concluding, we must reflect on the assumptions, limitations, and broader implications of the model. The final section addresses these critical questions."

---

## Section 7: Reflection and Conclusion

**Purpose and Claim:**
Reflect critically on the assumptions, limitations, and implications of the work. Revisit the research question and explain what was learned. Discuss next steps and open questions.

**Pages:** 1.5

**Key Content:**

**Subsection 7.1: Critical Assumptions Revisited**
- **Uniformity of states:** Assumed that last-slot and last-layer states are uniformly random after the first three pairs. Is this true? If certain F2L solutions are more common (e.g., faster rotationless solutions), state distributions may skew. Impact: would change the weighting of OLL/PLL cases.
- **Deterministic algorithms:** The model assumes each case has a pre-selected algorithm. In practice, a solver might choose algorithms on the fly (lookahead, ergonomics). Impact: the "alternate angles" technique overstates the gain if good angles are already chosen.
- **No rotations, no lookahead:** The F2L data excludes rotations (x, y, z moves) and assumes pre-selected solutions. Real solves include cube rotations and lookahead. Impact: the baseline movecount may not match wall-clock time.
- **FTM metric:** The paper uses FTM (U2 = 1), not half-turn metric (HTM, where U2 = 2) or others. This choice affects all comparisons. Discussion: FTM is closer to human solve time (since U2 is performed as a single rotation), but HTM is more granular.

**Subsection 7.2: Limitations of Scope**
- **Last layer only:** The analysis optimizes only the last slot and last layer. A full system would include the cross and F2L, where larger improvements may be possible.
- **No edge control in baseline:** The baseline does not explicitly control edge orientation. ZBLS, by contrast, builds edge control into the algorithm set. A fair comparison would require a parallel baseline with full edge control (e.g., ZB Full).
- **Case assumptions:** OLL has 58 cases (including skip), PLL has 22. But in some conventions (e.g., distinguishing rotations), the case count is different. The paper's results depend on this choice.
- **Model vs. practice:** The model describes an idealized solver. Real solves are affected by fatigue, turn speed, lookahead, and algorithm comfort. The improvements predicted by the model may not translate to faster real-world solves.

**Subsection 7.3: Strengths and Insights**
- **What the paper demonstrates:**
  1. The DFT is a powerful and elegant tool for combining independent distributions. The 4-point example and polynomial interpretation make this accessible.
  2. Independence is a critical assumption. Many interesting techniques violate it, requiring careful state-by-state or simulation-based analysis.
  3. Most last-layer techniques offer modest improvements (0.5–3 moves), with diminishing returns as complexity increases.
  4. The framework is general and can be extended to other optimization problems (e.g., other Rubik's cube methods, or other discrete combinatorial problems with step-by-step solutions).
- **Personal learning:** The author set out to learn the DFT because pairwise convolution felt inelegant. The journey from polynomial multiplication to roots of unity to inverse transform is mathematically rich and conceptually satisfying.

**Subsection 7.4: Limitations and Future Work**
- **Full enumeration:** A complete analysis would enumerate all cube states (very large; impractical) or sample more representative cases and report confidence intervals.
- **Real-solve validation:** Collect movecount data from real solves (via video analysis or move logging) and compare to the model. Quantify the lookahead effect.
- **Optimization across phases:** The paper isolates last-layer optimizations. A richer project would co-optimize last-pair insertion and OLL, using a shared state model.
- **Time vs. movecount:** Investigate the relationship between movecount and solve time (e.g., via video analysis or hand motion tracking). Are the improvements in movecount reflected in wall-clock time savings?
- **Other methods:** The framework applies to other speedcubing methods (Roux, Petrus, ZZ). Comparing methods via this distribution-based approach is an open question.

**Subsection 7.5: Answering the Research Question**
- **Question:** How can probability distributions be used to compare how much last-slot and last-layer techniques reduce the movecount of a CFOP solve?
- **Answer:** Probability distributions model the movecount of each step. For independent steps, convolution (or the DFT) combines them elegantly. For dependent techniques (choice rules, edge control, cancellations), state-by-state analysis or simulation is required. The full distribution (not just the mean) reveals the effect of techniques on the whole range of outcomes: some techniques reduce the mean, others reduce variability or eliminate outliers.
- **Key finding:** The DFT provides a conceptually rich way to understand independent combinations, but the most interesting speedcubing techniques are not independent, requiring deeper reasoning.

**Subsection 7.6: Conclusion**
- Summarize: The paper uses probability distributions and the DFT to model last-layer movecount. The analysis is rigorous for the baseline and transparent about where and why the model breaks for advanced techniques. Several techniques offer small improvements; ZBLS+ZBLL offers a larger alternative at high memorization cost.
- Reflection: Movecount is not time, and the improvements shown may not directly translate to faster solves. However, the framework offers a rigorous way to compare techniques and could guide solver practice.
- Final thought: The author's original motivation—learning the DFT to elegantly combine distributions—is achieved. The DFT is the right tool for the baseline, and the limitations exposed by advanced techniques point the way toward more sophisticated models (e.g., state machines, Markov chains, simulation).

**Math Shown vs. Summarized:**
Minimal new math. Recap the key formulas and insights from earlier sections.

**Visuals:**
Optional: a concept map or flowchart showing how the different techniques fit into the framework (independent baseline, state-dependent choice rules, simulation-based techniques, etc.).

**Data Needed:**
None beyond what was already collected.

**Personal Engagement:**
Reflect with the author: Was learning the DFT worth it? What surprised you most about the results? If you were to practice speedcubing based on this analysis, which techniques would you prioritize and why?

**Reflection:**
This section is the reflection in full. It addresses limitations, acknowledges the gap between the model and reality, and highlights the insights gained. The author should demonstrate understanding of the context (speedcubing as a practice and a skill) and the limitations of the analysis (ideal solver, last-layer focus, movecount as a proxy for time).

**Rubric Served:**
- Organization: clear conclusion linking back to the research question.
- Reflection: substantial and critical, throughout (this section in particular).
- Personal engagement: author's journey and learning are honored.
- Mathematical communication: recap of methods and their limitations.
- Use of mathematics: the DFT is positioned as elegant and powerful, with caveats about generality.

**Transition:**
The paper ends here. An appendix may include:
- Full algorithm lists for F2L, OLL, PLL, COLL, ZBLS, ZBLL.
- Raw data tables (case-by-case movecounts, if not in the main text).
- Simulation pseudocode or Python implementation.
- Derivation of the inverse DFT in detail.
- Worked examples of convolution and DFT for larger arrays.

---

## Appendices (not counted in page budget)

**Appendix A: Algorithm Data**
- F2L algorithm movecount for each of the $\sim 40$ last-pair cases.
- OLL and PLL case tables (58 OLL + 22 PLL cases, with movecount for primary angle).
- COLL case table (42 cases, movecount).
- ZBLS movecount per case (42 cases).
- ZBLL movecount per case (493 cases), or a summary (mean, std dev, distribution).

**Appendix B: Detailed Derivations**
- Full derivation of the inverse DFT.
- FFT algorithm (brief overview, not a full derivation).
- Proof of independence fact: $E[z^{X+Y}] = E[z^X] E[z^Y]$.
- Worked example: 8-point DFT.

**Appendix C: Simulation and Computation**
- Pseudocode or Python code for direct convolution.
- Pseudocode for DFT (direct, not FFT).
- Simulation pseudocode for case enumeration.
- Verification checks (probabilities sum to 1, means match, etc.).

**Appendix D: Real-Solve Data (if available)**
- Video analysis or move logging of the author's solves.
- Comparison of modeled vs. observed movecount.
- Statistical summary (mean, std dev, outliers).

---

## Quality Checks

- ✓ Clear central argument: the DFT elegantly handles independent baseline; most interesting techniques are not independent.
- ✓ Valid mathematical model for each technique: specified in Section 5.
- ✓ Each use of convolution justified by independence: stated in Section 3 and revisited in Section 5.
- ✓ Reflection and personal engagement spread throughout: author's motivation (Sec. 1), assumptions (Sec. 2), learning (Sec. 4), choices (Sec. 5), surprise (Sec. 6), and critical reflection (Sec. 7).
- ✓ Repetitive material minimized: each section has a distinct purpose; summary tables consolidate similar information.
- ✓ Full distributions (not just means) matter: shown in Section 3 (why threshold probabilities matter), Section 6 (histogram comparisons).
- ✓ Page budget: 15 pages fits the guideline.

---

## Open Questions for Author (to be answered later)

1. **U2 insert algorithms:** What exact algorithms are you calling "U2 insert"? How does it differ from the standard insert, and what is your preferred standard insert for comparison?
2. **Algorithm selection:** How is one algorithm chosen among alternatives for a given case? (Shortest, random among shortest, solver preference?) Does this change per technique?
3. **Metric finality:** Is FTM your final choice, or are you considering other metrics (HTM, QTM, STM)?
4. **Cancellations scope:** Are cancellations applied to all techniques, or tested separately as a distinct optimization?
5. **Two-step rule scope:** In which steps is the two-step minimization rule applied? (Last pair insertion only, OLL angle choice, PLL angle choice, or all?)
6. **Real-solve data:** Do you have recorded solves (video or move sequences) available for validation?
7. **Course expectations:** Are there specific requirements for line spacing, citation format, or a course rubric that constrains the paper style?

---

## Summary

This outline presents a comprehensive structure for a 15-page paper on CFOP last-layer optimization using probability distributions and the DFT. The paper's strength lies in honest treatment of the assumptions: the DFT is elegant and correct for independent steps, but the most interesting techniques break independence, requiring careful state-by-state or simulation-based reasoning. The structure encourages the author to engage deeply with the mathematics, reflect critically on the results, and communicate findings clearly to a reader unfamiliar with speedcubing.

The page budget is tight (15 pages), so material must be focused and repetition minimized. Large data tables and detailed proofs belong in appendices. The visual elements (histograms, tables, comparison plots) should carry much of the communication.

The outline is ready for the author to fill in with data, computations, and final decisions on the open questions above.
