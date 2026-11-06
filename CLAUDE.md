# CFOP Last-Slot-Last-Layer Paper: Outline Brief

## Project Status (Updated 2026-11-05)

**Completed:**
- ✅ Detailed 15-page outline (`OUTLINE.md`) with section-by-section guidance, page budgets, data requirements, and quality checks.
- ✅ Identified current codebase: spreadsheet (LSLL distributions.xlsx) with inelegant AUF handling, and structured JSON algorithm data (cube-simulator-and-algs/).
- ✅ Generated `questions.txt` with 15 clarification questions for author.

**Next phase (awaiting author answers):**
1. Confirm data sources (JSON vs. spreadsheet preference).
2. Define movecount metric precisely (FTM counting rules).
3. Clarify algorithm angle handling (single best, or distribution across angles).
4. Define AUF distribution parameters and application scope.
5. Explain magiccube integration (what is it, how to use it for simulation/testing).
6. Provide exact algorithms for U2 insert, sledge insert, standard baseline.
7. **Then:** Build clean Python codebase with separated Distribution classes, convolution pipeline, and simulation hooks.

**Key decisions made (for continuity):**
- The outline is comprehensive, following the CLAUDE.md rubric closely. It identifies where independence breaks and appropriate models for each technique.
- The existing spreadsheet combines AUFs inelegantly; the new codebase will track them separately.
- Algorithm JSON data in cube-simulator-and-algs/ is structured (cases, angles, scores, votes); need to decide on aggregation strategy.

---

## Your job

Create a detailed outline of the paper and save it as `OUTLINE.md`. Do not write the paper.

The author writes in their own voice and has a graded-exam mindset. Every section must earn marks under the rubric below. Do not invent data, results, algorithms, or personal anecdotes. Mark anything unavailable as `[NEEDED FROM AUTHOR]`.

Ask the author before finalizing if a blocking item in "Open questions" would change the structure. Otherwise, note your assumption and continue.

## Paper at a glance

**Working research question (the author may revise it):**
How can probability distributions be used to compare how much last-slot and last-layer techniques reduce the movecount of a CFOP solve?

**Scope:** Last slot and last layer only. The cross and first three pairs are excluded because the author considers them not deterministic enough to model well.

**Primary outcome:** Movecount, using one consistent metric. Current working assumption: face turn metric (FTM), where U2 counts as 1 move. Movecount is not solve time; address this difference in the reflection.

**Length:** Target about 15 pages, maximum 16. The guideline is 12-20 pages, but repetitive, tangential, or wordy material loses marks. Large tables and raw data go in appendices.

**Personal hook:** The author was dissatisfied with chaining pairwise convolutions and wanted to combine all distributions at once. That led them to learn the DFT by hand. Treat this as a genuine motivation, not decoration.

## Model

### Baseline sequence

1. Pre-AUF for last pair
2. Last pair
3. Pre-AUF for OLL
4. OLL
5. Pre-AUF for PLL
6. PLL
7. Post-AUF for PLL

AUF distribution from the author: P(0 moves) = 0.25, P(1 move) = 0.75. This requires U2 to count as 1 move.

### Data supplied by the author

- F2L algorithm data, which may differ from real solves. Note this in the reflection.
- Exact OLL, PLL, and AUF distributions. The author states that these match practice, so no general disclaimer is required.
- Alternate OLL and PLL angles will be simulated and aggregated into a distribution.

### Techniques to compare

| Technique | Working definition | Modeling note |
|---|---|---|
| Two-step choice rule | At a decision, choose the option minimizing the combined moves of this step and the next | The author's example: a sledge 2 moves longer that skips OLL is worth it |
| U2 insert | Insert variant, exact algorithms `[NEEDED FROM AUTHOR]` | No reliable public definition was confirmed in research |
| Sledge insert | Pair insertion using `R' F R F'`, which can change last-layer edge orientation | Compare against the standard insert |
| Alternate OLL and PLL angles | Using other algorithm angles to reduce total moves | Simulate |
| Cancellations | Moves at a step boundary merge or disappear, such as `U` + `U` -> `U2`, or `U` + `U'` -> nothing | Requires actual algorithm strings |
| COLL as an OLL choice | Available when edges are oriented; solves corners and leaves EPLL | Changes the next-step distribution |
| ZBLS + ZBLL vs OLL + PLL | Last slot plus edge orientation, then one-look last layer | Compare as a whole-system alternative |

### Research notes

- ZBLS has 302 cases, and with ZBLL totals 795 cases: https://speedsolving.com/wiki/index.php?title=ZBLS
- ZBLS solves the last F2L pair and orients last-layer edges. ZBLL then finishes the last layer: https://www.speedsolving.com/wiki/index.php/ZB_Method
- ZBLL case counts vary by convention. The wiki lists 493, and one derivation gets 494 including the skip. State the convention used: https://www.speedsolving.com/wiki/index.php/ZBLL and https://cuberoot.me/math/lsll
- Without edge control, oriented edges arise in 1 of 8 solves: https://www.speedsolving.com/wiki/index.php/ZBLL
- COLL orients and permutes corners while preserving edge orientation, leaving EPLL. It has 42 cases including mirrors, or 24 without: https://speedsolving.com/wiki/index.php?title=COLL
- Sledge insert and its edge-orientation effect: https://www.speedsolving.com/wiki/index.php/Partial_Edge_Control and https://www.speedsolving.com/wiki/index.php?title=Sledgehammer
- OLL has 58 cases and PLL 22 cases when each skip is included: https://cuberoot.me/math/lsll

## Mathematical requirements

### Critical constraint: convolution requires independence

Convolving steps is valid only when their movecounts are independent. Choice rules, cancellations, and alternate angles can link adjacent steps, so a simple convolution of marginal distributions can be wrong.

The outline must state where each tool is valid:

- **Independent baseline:** Convolution and the DFT are appropriate after verifying the assumption.
- **Choice rules:** Model each state's selected option, for example:
  `P(T = m) = Σ_s P(s) · 1[min_option total moves in state s = m]`.
  This is a mixture over states, not a product of unrelated distributions.
- **Cancellations and angles:** Enumerate or simulate the actual sequences.
- **Conditioned steps:** Where appropriate, convolve within a fixed state, then mix across states.
- **Technique comparisons:** Do not treat two techniques on the same scramble as independent. For `P(A beats B)` or moves saved per state, use paired state-by-state comparison.

Treat this limitation as a major analytical turning point: the DFT elegantly solves the author's original problem, but the most interesting techniques break the assumption it relies on. This is a key opportunity for sophistication and critical reflection.

### Justify the DFT honestly

Do not claim the DFT is required for speed at these small array sizes, where direct convolution is manageable. Justify it through structure: independent steps multiply at each evaluation point. Compare it fairly with direct convolution and simulation. Mention the FFT only briefly unless it becomes relevant.

### Required development

Include these only if they serve the argument:

1. Distribution as polynomial; convolution as polynomial multiplication.
2. Why values multiply: for independent `X, Y`, `E[z^(X+Y)] = E[z^X] E[z^Y]`.
3. Padding: `N >= (sum of maximum movecounts) + 1`, plus the wraparound error if `N` is too small.
4. DFT as evaluation at roots of unity.
5. Inverse DFT by solving the evaluation equations through cancellation, including why dividing by `N` appears.
6. One compact hand-worked 4-point example. Preferred: three flat AUFs, `[0.25, 0.75, 0, 0]`, which also connect to the binomial distribution.
7. Means add under linearity. This also shows why full distributions matter: tails and thresholds such as `P(T <= m)` need more than the mean.

The author learned this from polynomial evaluation and systems of equations. They are comfortable with the process but have not fully connected it to the summation formula. Plan one short bridge from the table method to the formula, and flag it as something the author must understand before submission.

### Things to verify

- **Symmetry and skips:** Flat AUFs may conflict with skips and symmetric cases. An OLL skip needs no pre-OLL AUF, and a PLL skip needs no pre-PLL AUF. Check the exact data for this, including whether it already embeds AUFs.
- **Double counting:** Confirm that F2L algorithm data excludes setup AUFs if they are modeled separately.
- **Metric:** Every data source and algorithm must use the same metric.
- **State assumption:** Justify or test the assumption that last-slot and last-layer states are uniformly random after three pairs.
- **Decision rule:** Define a tie-break rule and state that a two-step rule is greedy, so it may miss a better three-step choice.
- **Checks:** Probabilities sum to 1, the distribution mean matches a separately computed mean, and simulated results converge or are fully enumerated.

## Comparison and validation plan

For each technique, plan to report:

- Full movecount distribution with labeled axes
- Mean, standard deviation, and threshold probabilities
- Mean moves saved against the baseline
- Paired comparison where possible
- Algorithms required, plus moves saved per additional algorithm

Keep the comparison concise with one overlay plot and one summary table.

Validation options:

- Compare analytic results to full enumeration or simulation, with stated precision.
- If the author has reconstructed solves, compare the model with real solves. Otherwise, state that the model describes an idealized solver.

Ask the author before adding any real-solve data collection.

## Rubric (condensed)

**1. Organization: coherent, well organized, concise**
Include an introduction, clear aim, body, and conclusion, all logically linked. Put relevant graphs and tables where they are discussed; use appendices only for large data or extras. Avoid irrelevant material and repetitive calculations.

**2. Mathematical communication: relevant, appropriate, consistent**
Use correct notation and terminology, avoiding calculator or computer notation unless software-generated. Define key terms and variables. Use appropriate representations such as formulas, diagrams, tables, and graphs. Present proofs logically. All graphs must be labeled.

**3. Personal engagement: outstanding**
Show authentic, high-quality engagement in numerous places, driving the work forward creatively. Examples include independent thinking, original presentation, multiple perspectives, and making and testing predictions. Avoid textbook-style reproduction. The work must leave the impression of complete understanding of the context. Effort alone does not count.

**4. Reflection: substantial and critical, throughout**
Go beyond describing results. Link to aims, consider limitations, compare approaches, discuss strengths, weaknesses, and implications, consider other perspectives, and identify next steps. Reflection should develop the exploration, not appear only at the end.

**5. Use of mathematics: relevant, precise, sophisticated, rigorous, thoroughly understood**
Mathematics must be appropriate to the course and support the aim. Avoid overcomplication where simpler mathematics works. Be error-free with appropriate accuracy. Justify or prove relevant claims, and demonstrate understanding through reasoning and examples rather than correct answers alone. Sophistication can include challenging concepts, multiple perspectives, and links between areas of mathematics. Technology is encouraged, but understanding must be shown. Doing a few things well beats doing many things poorly.

## Required outline format

Start `OUTLINE.md` with:

- The final research question
- A short statement of the central argument
- A page-budget table totaling 15-16 pages
- A list of data still needed

Then provide this for each section and subsection:

- **Purpose and claim:** What this section establishes
- **Pages:** Budget
- **Key content:** Specific definitions, derivations, calculations, or simulations
- **Math shown vs. summarized:** What appears in full and what moves to an appendix
- **Visuals:** Specific graph, table, or diagram, including axes
- **Data needed:** Exact source or `[NEEDED FROM AUTHOR]`
- **Personal engagement:** A real opportunity for the author to add their own decision or prediction, without invented details
- **Reflection:** Limitation or decision raised at that point
- **Rubric:** Criteria served
- **Transition:** How it leads to the next section

Suggested page budget, which you may adjust with a reason:

| Section | Pages |
|---|---:|
| Introduction, aim, prediction | 1 |
| Scope, definitions, model assumptions | 2 |
| Baseline distributions | 2 |
| Combining distributions: polynomials and DFT | 3.5 |
| Techniques and the failure of independence | 3 |
| Results and comparison | 2 |
| Validation, reflection, conclusion | 2 |

## Quality checks before finishing

- Is there a clear central argument, not just a list of techniques?
- Does every technique have a valid mathematical model?
- Is each use of convolution justified by independence?
- Are reflection and personal engagement spread throughout?
- Is there anything repetitive that could go in a table or appendix?
- Can a reader see why full distributions, rather than means alone, matter?
- Does the page budget stay within 16 pages?

## Open questions for the author (answered later after first pass)

1. What exact U2 insert algorithms are meant, and which standard insert is the comparison?
2. Which algorithm set defines each step's movecount, and how is one chosen among alternatives?
3. Is FTM the final metric?
4. Are cancellations applied to every technique, or tested as a separate factor?
5. Is the "next two steps" rule applied at every decision, including the last-layer steps?
6. Are real reconstructed solves available for validation?
7. Does the assignment specify line spacing, citation style, or course-level expectations?

---

## Instructions for Next Agent

**When the author provides answers to questions.txt:**

1. **Read and organize clarifications:**
   - Update this CLAUDE.md with any new information (metric finality, algorithm definitions, magiccube details).
   - Note any scope changes or new priorities.

2. **Design and build the Python codebase:**
   - Structure: `src/distributions/`, `src/convolution/`, `src/simulation/`, `src/utils/`
   - Core components:
     - `Distribution` class: movecount → probability mapping, with operations (convolve, DFT, mean/std/quantiles).
     - `Metric` class: handles FTM parsing from algorithm strings.
     - `AlgorithmSet` class: loads and aggregates OLL/PLL/COLL/etc. from JSON and/or spreadsheet.
     - `Convolution` pipeline: direct convolution + DFT, configurable by technique.
     - `Simulator` module: integrate magiccube (or alternative) for testing sledges, U2 inserts, edge-orientation effects.
   - Data layer: parse JSON files cleanly, extract movecounts, handle multiple angles (by score/vote weighting or enumeration).
   - Config file: AUF parameters, metric choice, data sources (JSON vs. spreadsheet), algorithm selection strategy.

3. **Verify and validate:**
   - Parse movecount from algorithm strings and compare against magiccube (once integrated).
   - Test sledge insert edge-orientation changes.
   - Test U2 insert algorithms (once defined).
   - Verify convolution against direct enumeration on a few cases.

4. **Generate baseline distributions:**
   - Compute independent baseline (pre-AUF, LS, OLL, PLL, post-AUF) using the chosen method (JSON data).
   - Convolve to final baseline distribution.
   - Report mean, std dev, quantiles.

5. **Implement each technique:**
   - Two-step choice: state-by-state enumeration.
   - Sledge + COLL + ZBLS+ZBLL: conditional/mixture models.
   - Cancellations: simulation (sequence-level analysis).
   - Alt. angles: distribution aggregation.
   - Create comparison tables and visualizations.

6. **Output:**
   - `data_distributions.py`: loaded and cached distributions (baseline for each step, combined baseline, per-technique results).
   - `convolution.py`: convolution and DFT logic.
   - `simulation.py`: magiccube integration and test harness.
   - `results.json` or pickle: final distributions for each technique (for paper figures).
   - `validation_report.txt`: comparison vs. enumeration/simulation, confidence intervals.
   - Update CLAUDE.md with any changes.

**Blocking issues:**
- If magiccube is unavailable or incompatible, use alternative (e.g., `rubiks-cube` Python package, or pure algorithm string parsing + manual verification).
- If JSON and spreadsheet data conflict, consult author; prioritize JSON unless told otherwise.
- If movecount metric is ambiguous, implement both FTM and HTM, with config switch.

**Timeline estimate:**
- Data loading + parsing: 2–3 hours.
- Distribution + convolution module: 2–3 hours.
- Baseline computation: 1 hour.
- Technique implementations: 4–6 hours (state enumeration, simulation, aggregation).
- Testing + validation: 2–3 hours.
- **Total: ~15–20 hours** (parallel work welcome).
