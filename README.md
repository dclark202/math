# Math

Mathematical writings of Duncan Clark. Two groups here: research and course notes. 
- Research is two papers I wrote as a graduate student *On the Goodwillie derivatives of the identity in structured ring spectra* and *The partition poset complex and the Goodwillie derivatives of the identity in spaces*, along with my PhD thesis. Both papers were revised in August/September 2026 with the help of Claude Fable 5/5.1, resulting in an errata list (with genuine mathematical corrections) to the former, and a largely expanded version of the second. *The derivatives of the identity and the spectral Lie opera* is an additional memo which asserts equivalence of the operad I produce in the second paper with the well-known spectral Lie operad.

- Course notes consists of three textbooks I compiled for courses I frequently taught at Milwaukee School of Engineering

## Research

**[On the Goodwillie derivatives of the identity in structured ring spectra](derivatives-id-algO/dac-der-id-algO.pdf)** · 48 pp. · [arXiv:2004.02812](https://arxiv.org/abs/2004.02812) · Tbilisi Math. J., special issue on homotopy theory, spectra, and structured ring spectra.
Constructs a natural highly homotopy coherent operad structure on the derivatives of the identity on algebras over an operad `O` in spectra, shows every connected `O`-algebra carries a left action of those derivatives, and proves the resulting operad is equivalent to `O` itself. Introduces **N**-colored operads with levels as the framework for the comparison.

**[Errata and corrections to the above](derivatives-id-algO/errata.pdf)** · 21 pp. · September 2026.
Lists the corrections to the first paper by section and supplies a corrected framework for its Sections 6–8. The main results, Theorem 1.1(a) and (b), stand, but their published proofs do not: the composition product of **N**-colored objects with levels is not associative (Proposition 6.17 is false), the comparison map θ of Section 5.5 is not invertible even for Σ-copowered objects (Remark 7.9), and the proof of Proposition 7.16 uses maps that do not exist. The note replaces profiles by leveled trees, defines operads with levels through their component maps, brackets iterated box-circle products from the right so that θ itself supplies the compositions, and re-proves the main theorems in that setting. Prepared with the assistance of Claude; the finite combinatorial claims involved were checked by computer.

**[The partition poset complex and the Goodwillie derivatives of the identity in spaces](derivatives-id-spaces/dac-der-id-spaces.pdf)** · 47 pp. · [arXiv:2007.05440](https://arxiv.org/abs/2007.05440)
Gives a new construction of an operad structure on the derivatives of the identity in spaces via a pairing of cosimplicial objects, and proves it agrees with the structure originally described by Ching — both are restrictions of a single algebra over an operad of windowed cut systems on weighted trees. The version here is the one uploaded to the arXiv in September 2026; it adds the comparison appendix (Appendix A) to the 2020 version. As disclosed in its acknowledgements, the proof of the comparison theorem was developed with substantial assistance from Anthropic's Claude models. The identities that proof relies on are re-run in small arities by the pure-Python scripts in [`derivatives-id-spaces/audit-scripts/`](derivatives-id-spaces/audit-scripts/README.md); each script ends with an `ALL ... PASSED` line.

**[The derivatives of the identity and the spectral Lie operad](derivatives-id-spaces/memo-infty-homology/infty-homology-memo.pdf)** · 14 pp. · memo, September 2026.
A companion note to the second paper, written with substantial assistance from Claude and not intended for publication. It reduces the statement that both point-set operad structures on the derivatives present the Blans–Heuts spectral Lie operad, in the ∞-categorical sense, to a single residual lemma, with a ledger of what is cited, proved, and open; and it verifies the comparison theorem on homology through arity 6 by an independent computation, identifying the homology operad as the operadic suspension of the Lie operad. The scripts behind that computation are in [`derivatives-id-spaces/audit-scripts/homology/`](derivatives-id-spaces/audit-scripts/homology/README.md).

**[PhD thesis](phd-thesis/dac-phd-thesis.pdf)** · 143 pp. · Ohio State University, 2021.
Develops the same material at length: operads and their algebras, an overview of functor calculus, the box product, **N**-colored operads with levels, and the derivatives of the identity, with two appendices. The two papers above are the later and more polished treatments; the thesis is included as the fuller account of the background.


## Course notes

Three self-contained undergraduate texts, each with worked examples and exercises at a range of difficulty (challenge problems are marked with a `*`). Written to accompany a lecture course rather than replace one.

**[Abstract Algebra and Number Theory](course-notes/Abstract-algebra/Abstract-algebra.pdf)** · 74 pp.
Divisibility and the Euclidean algorithm, modular arithmetic, RSA encryption, groups, symmetric and dihedral groups, normal subgroups and quotients, group actions and counting, matrix groups.

**[Linear Algebra](course-notes/Linear-algebra/Linear-algebra.pdf)** · 93 pp.
Linear systems and row reduction, vector spaces, determinants, eigenvalues and diagonalization, Markov chains, inner products and orthogonality, least squares, singular value decomposition.

**[Multivariable Calculus](course-notes/Multivariable-calculus/Multivariable-calc.pdf)** · 99 pp.
Vectors, curves and surfaces, curvature and torsion, partial derivatives, optimization and Lagrange multipliers, multiple integrals in several coordinate systems, vector fields, and Green's, Stokes' and the divergence theorems.

Alongside the algebra notes, [`course-notes/Abstract-algebra/Python script/`](course-notes/Abstract-algebra/Python%20script/) holds fourteen Jupyter notebooks implementing the algorithms as they come up: the Euclidean algorithm, the Chinese remainder theorem, Euler's totient, primitive roots, Cayley tables, and working RSA, ElGamal and Caesar ciphers.

## License

The three sets of course notes are licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) — use and adapt them for non-commercial purposes with attribution, sharing any derivative under the same terms.
