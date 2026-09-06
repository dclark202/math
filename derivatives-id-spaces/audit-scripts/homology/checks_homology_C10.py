"""C10: homology of the partition poset complex P(n) (Definition 2.4) and its Sigma_n-character.

For n = 1, ..., 6 (`python checks_homology_C10.py [nmax]`, default 6, about 15 s) this script asserts
  (a) H~_i(P(n); Z) = 0 for i != n-1 (zero rank AND no torsion) and H~_{n-1}(P(n); Z) = Z^{(n-1)!};
  (b) the character of Sigma_n on H~_{n-1}(P(n); Q), obtained as traces of the matrices of the action on an
      explicit basis of cycles (chains.PChainComplex), equals the character of Lie(n) (x) sgn, where
      Lie(n) = Ind_{C_n}^{Sigma_n}(omega) is computed independently by brute force from the induced-character
      formula (chains.induced_lie_character; this is Klyachko's description of Lie(n)).
Along the way: boundary^2 = 0; equivariance of the boundary; the cycle basis is integral (hence a Z-basis of
the integer cycles); M(e) = I and M(s t) = M(s) M(t) on homology; traces are class functions; the Hopf trace
formula (trace on H~_{n-1} = (-1)^{n-1} sum_k (-1)^k #{basis chains of C_k fixed by s}, valid once (a) is
known) reproduces the character without any matrices; and for n <= 5 a dense Smith normal form of every
boundary matrix agrees with the sparse computation.  For n <= 5 every permutation-indexed check runs over
all of Sigma_n; for n = 6 they run over all class representatives, the generators of Sigma_6 and a random
sample (the character comparison itself is complete: every conjugacy class).
Ends with `ALL C10 CHECKS PASSED` or raises at the first failure.
"""
import random
import sys
import time
from math import factorial

from chains import (PChainComplex, conjugacy_classes, cycle_type, format_cycle_type, cycle_notation, generators,
                    identity, compose, sign, induced_lie_character, klyachko_closed_form, mat_mul, trace)

random.seed(1)


def check_n(n):
    t0 = time.time()
    full = n <= 5                                     # run permutation-indexed checks over all of Sigma_n
    cc = PChainComplex(n)
    cc.check_dd()
    classes = conjugacy_classes(n)
    perms = [s for cl in classes.values() for s in cl]
    reps = [cl[0] for cl in classes.values()]
    assert len(perms) == factorial(n)
    e = identity(n)

    def dedupe(seq):
        seen, out = set(), []
        for s in seq:
            if s not in seen:
                seen.add(s)
                out.append(s)
        return out

    subset = perms if full else dedupe([e] + reps + generators(n) + random.sample(perms, 30))
    cc.check_equivariance(subset)                    # contains the generators, so this covers all of Sigma_n

    # (a) homology over Z
    hom = cc.homology(dense_snf=full)
    print("n = %d: nondegenerate non-basepoint simplices per degree: %s"
          % (n, ", ".join("C_%d = %d" % (k, cc.dim[k]) for k in cc.degrees)))
    for k in cc.degrees:
        h = hom[k]
        assert h["rank_Q"] == h["rank_Z"]
        assert not h["torsion"], ("torsion in H~_%d" % k, h["torsion"])
        if k != n - 1:
            assert h["betti"] == 0, ("H~_%d has rank %d" % (k, h["betti"]))
    assert hom[n - 1]["betti"] == factorial(n - 1), (hom[n - 1]["betti"], factorial(n - 1))
    print("  H~_i(P(%d); Z) = 0 for i != %d, H~_%d = Z^%d, no torsion in any degree (all Z-echelon pivots are units: %s%s)"
          % (n, n - 1, n - 1, hom[n - 1]["betti"], all(hom[k]["unit_pivots"] for k in cc.degrees),
             "; dense Smith normal forms agree" if full else ""))

    # (b) action on the top homology
    free, K, integral = cc.top_cycle_basis()
    assert len(K) == factorial(n - 1)
    assert integral, "cycle basis is not integral"
    maxentry = max(abs(c) for v in K for c in v.values())
    mats = {s: cc.top_action_matrix(s, verify=True) for s in (perms if full else dedupe([e] + reps))}
    assert mats[e] == [[1 if i == j else 0 for j in range(len(K))] for i in range(len(K))]
    if n <= 4:
        pairs = [(s, t) for s in perms for t in perms]
    elif n == 5:
        pairs = [(random.choice(perms), random.choice(perms)) for _ in range(300)]
    else:
        pairs = [(random.choice(perms), random.choice(perms)) for _ in range(6)]
    for s, t in pairs:
        Ms = mats[s] if s in mats else cc.top_action_matrix(s, verify=False)
        Mt = mats[t] if t in mats else cc.top_action_matrix(t, verify=False)
        Mst = mats[compose(s, t)] if compose(s, t) in mats else cc.top_action_matrix(compose(s, t), verify=False)
        assert mat_mul(Ms, Mt) == Mst, "not a left action on homology"
    chi_hom = {}
    for ct, cl in classes.items():
        chi_hom[ct] = trace(mats[cl[0]])
        assert chi_hom[ct] == cc.top_trace(cl[0])
    for s in subset:                                 # the trace is a class function
        assert cc.top_trace(s) == chi_hom[cycle_type(s)], ("trace not constant on the class", s)
    # Hopf trace formula cross-check (valid because H~_k = 0 for k != n-1 was verified above)
    for ct, cl in classes.items():
        for s in (cl[:3] if full else cl[:1]):
            fc = cc.fixed_counts(s)
            assert (-1) ** (n - 1) * sum((-1) ** k * fc[k] for k in cc.degrees) == chi_hom[ct], ("Hopf trace", ct)

    # Lie(n) by brute-force induction, then (x) sgn
    chi_lie = induced_lie_character(n)
    lie = {}
    for ct, cl in classes.items():
        vals = {chi_lie[s] for s in cl}
        assert len(vals) == 1
        lie[ct] = vals.pop()
    assert lie[tuple([1] * n)] == factorial(n - 1)

    print("  character table: H~_%d(P(%d); Q) vs Lie(%d) (x) sgn; cycle basis: %d integral cycles, max |entry| %d, "
          "%d matrices verified, %d products checked" % (n - 1, n, n, len(K), maxentry, len(mats), len(pairs)))
    print("  %-12s %-16s %6s %4s %10s %8s %14s %6s %12s"
          % ("cycle type", "representative", "size", "sgn", "chi(H~top)", "chi(Lie)", "chi(Lie x sgn)", "equal", "closed form"))
    closed_ok = True
    for ct, cl in classes.items():
        sg = sign(cl[0])
        expected = sg * lie[ct]
        cf = klyachko_closed_form(ct)
        closed_ok &= (cf == lie[ct])
        print("  %-12s %-16s %6d %4d %10d %8d %14d %6s %12d"
              % (format_cycle_type(ct), cycle_notation(cl[0]), len(cl), sg, chi_hom[ct], lie[ct], expected,
                 "yes" if chi_hom[ct] == expected else "NO", cf))
    print("  closed-form Klyachko values mu(m) m^(d-1) (d-1)! agree with the brute-force induction: %s"
          % ("yes" if closed_ok else "NO"))
    for ct in classes:
        assert chi_hom[ct] == sign(classes[ct][0]) * lie[ct], ("character mismatch", n, ct, chi_hom[ct], lie[ct])
    twisted = [ct for ct in classes if sign(classes[ct][0]) * lie[ct] != lie[ct]]
    print("  classes on which Lie(%d) and Lie(%d) (x) sgn differ: %s"
          % (n, n, ", ".join(format_cycle_type(c) for c in twisted) or "none"))
    print("  time %.2fs" % (time.time() - t0))


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    T0 = time.time()
    for n in range(1, nmax + 1):
        check_n(n)
    print("total time %.1fs" % (time.time() - T0))
    print("ALL C10 CHECKS PASSED")
