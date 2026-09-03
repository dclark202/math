# Numerical checks for Appendix A and Section 4 of `dac-der-id-spaces.tex`

These scripts re-run, on explicit small cases, every "direct check" that the proofs in Appendix A (and Lemma 4.4 / Proposition 4.6) rely on. They are evidence, not proofs: each script enumerates or samples chains of partitions, weighted trees, flags and cells for arities `t ≤ 5`, evaluates both sides of the identity in question with the formulas of the paper, and asserts equality (to 1e-9) or the stated property. Pure Python 3 (no third-party packages); run `python <script>.py` from this directory. Every script ends with an `ALL ... PASSED` line, or a traceback pointing at the first failing case.

## Models (imported by the checks)

| file | implements |
|---|---|
| `partitions.py` | chains of partitions `P(n)` (Def 2.4), faces/degeneracies, profiles, the decomposition `Ψ` and composition `γ ∘ (β_i)` of Lemma 4.4 |
| `trees.py` | laminar trees (Def A.6), weightings (Def A.10), reduction (Lemma A.16), alive partitions (Def A.19), the levelling map `Θ` (Constr. A.22) and its inverse (proof of Thm A.24) |
| `cuts.py` | flags, positions, pieces (§A.46), framed cut systems and their cut maps with the thresholds of Def A.48, tree frames and labelled frames (Example A.50), Ching's cut `c^tr` (Def A.31) |
| `wcs.py` | windowed cut systems (Def A.54): tree cells, constant cells (Example A.55), substitution (A.25), anchors (Def A.52), the three legs of the contraction (Prop A.57), the composition of labellings (Def A.62) |

## Checks

| script | statements checked | what is verified | size |
|---|---|---|---|
| `checks_theta.py` | Lemma A.20, (A.20), Lemma A.23, Thm A.24 | `Θ` respects all face/degeneracy relations, is a bijection (round trips in both directions), and Lemma A.23(1)–(3) | all nondegenerate chains of `P(n)` and all trees, `n ≤ 5` |
| `checks_cuts.py` | Def A.31, Example A.50(1)(2), coassociativity (Thm A.34) | (A.21)/(A.22); levelled cut non-basepoint iff `alive⟨a⟩ = α`; cut-then-cut = cut along the 3-flag | all trees `n ≤ 5`, random weightings |
| `checks_continuity.py` | Prop A.49 / A.58 as originally stated | the counterexamples: admissible systems whose cut maps are discontinuous (this is what motivated axioms (W5), (W6)); the paper's own cells behave at the same limits | arities 2, 3 |
| `checks_wcs.py` | Def A.56, Prop A.57, Thm A.60 (A.26), Prop A.57 (`o`), Lemma A.63 (`σ`) | (W2)–(W4) of composites; associativity of substitution (two bracketings); cut of a composite = cut then cut; `o` and `σ` respect composition | all trees `n ≤ 5`, random flags and cells |
| `checks_legs.py` | continuity of the paper's cells | tree, constant, leg-1/2/3 cells and composites tend to the basepoint along sequences approaching every flag-absence wall and the basepoint of `B(t)` | 2944 sequences per cell type |
| `probe_w5.py` | axioms (W5), (W6) (repaired Def A.54) | the ratios `w_ν / min(v−u, E_ν−h(p(B)), w_anc)` stay bounded below and `(E_anc−E_ν)/w_ν` bounded above on tree, constant, leg cells and all their pairwise composites, including near-degenerate heights | all trees `n ≤ 5` |
| `checks_w56.py` | closure of (W5), (W6) (repaired Prop A.57) | (W5), (W6) hold with the *explicit constants* derived in the proof: composites (`λλ'`, `(K+1)/λ'`), Legs 1–2 (`min(λ, (1/λ+K)^{-1}, 1)`), Leg 3 (`t`) | 7370 positions × 25 cell types |
| `checks_transport.py` | Prop A.71, Lemma A.64, Constr. A.69 | chain-side levelled action (iterated `Ψ` + rescalings, peeled for general labellings) transported by `Θ` equals `κ_{F^e} ∘ Θ`, including degenerate simplices, unrelated flags and tied heights; `κ_{W^e} = κ_{F^e}`; coface compatibility of the multi-cut map | `t ≤ 5`, 1700 labellings, 11993 coface cases |
| `checks_sec4.py` | Lemma 4.4, Prop 4.6, line "trivial unless λ_p = λ_{p+1}" in the proof of Thm 1.2 | uniqueness of the decomposition, `compose ∘ decompose = id`, the set-level squares (4.14), associativity of decompositions | all `k`-simplices of `P(t)`, `t ≤ 5`, `k ≤ 4` (74374 squares) |

The full audit that produced these scripts is in `../AUDIT.md`; the repair of Definition A.54 / Proposition A.58 is described in `../REPAIR-NOTES.md`.
