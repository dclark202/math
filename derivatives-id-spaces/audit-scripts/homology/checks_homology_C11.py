"""Check C11: the two operad structures compared in Theorem 1.3 induce the SAME maps on the top homology of the
partition complexes.

For every partition alpha of {1..k}, k <= KMAX (default 5; `python checks_homology_C11.py [kmax]`), the based maps
    c^lev_alpha, c^tr_alpha : |P(k)| --> |P(n)| ^ |P(k_1)| ^ ... ^ |P(k_n)|
(the levelled cut = half-cut cell of the box-product structure of Theorem 1.2, and Ching's cut transported through
the levelling homeomorphism Theta; see degree.py for the formulas and the ORIENTATION CONVENTIONS (O1)-(O5)) are
pushed forward on H~_{k-1} by degree counting over the preimages of random generic points, in the cycle bases of
chains.py (source) and the Kunneth bases (target).  Before C11 is asserted the following guards must pass:
  (a) unit profiles (1;k) and (k;1..1): both induced maps are s * identity, s a global sign, which is reported;
  (b) Sigma_k-equivariance, pointwise (forward maps on random points) and on homology:
      M_{sigma.alpha} A(sigma) = T(alpha, sigma) M_alpha, with the Koszul sign of the factor reordering;
  (c) arity 3 reproduces Appendix B: the coordinate formulas (B.1) c^tr(x,y) = (x, (y-x)/(1-x)) on x < y and
      (B.2) c^lev(x,y) = (2x, 2y-1) on x < 1/2 <= y, the basepoint elsewhere, the worked input (1/4, 3/4)
      |-> (1/4, 2/3) resp. (1/2, 1/2), and the induced maps H~_2(P(3)) -> Z read off from these formulas;
  (d) the degree sums are integers independent of the generic point (3 random points per target cell, inside
      degree.induced_matrix), and the image of every source cycle is certified to be a combination of the
      Kunneth product cycles (so it is a cycle of the target);
  (e) the shuffle sign recorded for every affine piece agrees with the sign of the Jacobian determinant of the
      forward map recomputed by exact central differences (all pieces for k <= 4, a sample for k >= 5).
Ends with ALL C11 CHECKS PASSED or raises at the first failure."""
import os
import sys
import time
import random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
from degree import (Complexes, all_partitions, induced_matrix, jacobian_sign_fd, target_action_matrix,   # noqa: E402
                    act_target_point, matrices_equal, fmt_matrix, fmt_partition, target_basis, blocks_of,
                    FORWARD, POINT, random_heights, lev_forward, tr_forward, HALF)
from chains import all_perms, generators, act_chain, act_partition, mat_mul, cycle_notation  # noqa: E402
from partitions import min_partition, max_partition  # noqa: E402

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
rng = random.Random(2026)
cxs = Complexes()
T0 = time.time()
STRUCTS = ("lev", "tr")


def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def scaled(M, s):
    return [[s * x for x in row] for row in M]


def check_unit(M, k, s_name, which):
    """(a): M = s * I for a global sign s; returns s."""
    m = len(M)
    assert m == len(M[0]) == len(cxs(k)._top[1]), ("unit map has the wrong shape", k, s_name, which)
    s = M[0][0]
    assert s in (1, -1) and matrices_equal(M, scaled(identity(m), s)), \
        ("unit profile %s, structure %s, k = %d: not +-identity:\n%s" % (which, s_name, k, fmt_matrix(M)))
    return s


def random_top_point(k):
    cc = cxs(k)
    return rng.choice(cc.basis[k - 1]), random_heights(k - 1, rng)


def perm_sample(k):
    """All permutations for k <= 4; the generators plus a random sample for k >= 5."""
    if k <= 4:
        return all_perms(k)
    ps = list(generators(k))
    extra = 12 if k == 5 else 4
    while extra > 0:
        p = list(range(1, k + 1))
        rng.shuffle(p)
        ps.append(tuple(p))
        extra -= 1
    return ps


# ---------------------------------------------------------------------------------------------------------
# (c) Appendix B, arity 3, checked before anything else
# ---------------------------------------------------------------------------------------------------------
def check_appendix_B():
    T3 = frozenset([1, 2, 3])
    mn, mx = min_partition(T3), max_partition(T3)
    two = (frozenset([frozenset([1, 2])]), frozenset([frozenset([1]), frozenset([2])]))   # the top cell of P(2)
    cc3 = cxs(3)
    free, K, integral = cc3._top
    assert integral
    n_pts = 0
    for i in (1, 2, 3):
        ibar = T3 - {i}
        alpha_i = frozenset([ibar, frozenset([i])])
        gamma_i = (mn, alpha_i, mx)                         # nabla^2_i, (u_1, u_2) = (x, y)
        others = [(mn, frozenset([T3 - {j}, frozenset([j])]), mx) for j in (1, 2, 3) if j != i]
        slot = blocks_of(alpha_i).index(ibar)               # which smash slot carries B(ibar)

        def expect(bottom_height, top_height):
            pieces = [(POINT, ())] * 2
            pieces[slot] = (two, (top_height,))
            return (two, (bottom_height,)), tuple(pieces)
        samples = [(Fraction(1, 4), Fraction(3, 4))]
        while len(samples) < 60:
            samples.append(random_heights(2, rng))
        for _ in range(60):                                 # coarser samples, more often near the walls of (B.2)
            x = Fraction(rng.randint(1, 999), 1000)
            y = Fraction(rng.randint(1, 999), 1000)
            if x < y:
                samples.append((x, y))
        for x, y in samples:
            assert 0 < x < y < 1
            u = (x, y)
            tr = tr_forward(gamma_i, u, alpha_i)
            assert tr == expect(x, (y - x) / (1 - x)), ("(B.1) fails", i, x, y, tr)
            lev = lev_forward(gamma_i, u, alpha_i)
            if x < HALF < y:                                # generic part of the validity square of (B.2)
                assert lev == expect(2 * x, 2 * y - 1), ("(B.2) fails", i, x, y, lev)
            elif x > HALF or y < HALF:
                assert lev is None, ("(B.2) should be the basepoint off the square", i, x, y, lev)
            for g in others:                                # the other strata carry no vertex ibar
                assert tr_forward(g, u, alpha_i) is None and lev_forward(g, u, alpha_i) is None
            n_pts += 1
        # the worked input of Appendix B
        x, y = Fraction(1, 4), Fraction(3, 4)
        assert tr_forward(gamma_i, (x, y), alpha_i) == expect(Fraction(1, 4), Fraction(2, 3))
        assert lev_forward(gamma_i, (x, y), alpha_i) == expect(HALF, HALF)
        # induced map on H~_2(P(3)) = Z^2 read off from (B.1)/(B.2): both formulas have positive Jacobian
        # (1/(1-x) resp. 4) and cover the single target cell exactly once, so c_*(z) = z[gamma_i].
        expected = [[K[a].get(cc3.index[2][gamma_i], 0) for a in range(len(K))]]
        for s in STRUCTS:
            M, info = induced_matrix(s, 3, alpha_i, cxs, rng)
            assert matrices_equal(M, expected), ("Appendix B induced map", s, i, M, expected)
    print("(c) Appendix B: (B.1), (B.2), basepoint loci, worked input (1/4,3/4) -> (1/4,2/3) / (1/2,1/2), and the "
          "induced maps z -> z[nabla^2_i] reproduced at %d points x 3 cells" % n_pts)


check_appendix_B()

# ---------------------------------------------------------------------------------------------------------
# main loop
# ---------------------------------------------------------------------------------------------------------
results = {}
unit_signs = {}
stats = {}
for k in range(2, KMAX + 1):
    tk = time.time()
    cc = cxs(k)
    parts = all_partitions(k)
    rank = len(cc._top[1])
    print("\n== k = %d: %d partitions, %d top cells, rank H~_%d(P(%d)) = %d" % (k, len(parts), cc.dim[k - 1], k - 1, k, rank))
    n_fd = 0
    pre_range = {s: [None, None] for s in STRUCTS}
    for alpha in parts:
        Ms = {}
        for s in STRUCTS:
            M, info = induced_matrix(s, k, alpha, cxs, rng, npoints=3, stats=stats)
            Ms[s] = M
            S = info["structure"]
            # (d) preimage counts: exactly one for the levelled cut (one shuffle fits the generic heights)
            for e, counts in info["npre"].items():
                for c in counts:
                    if s == "lev":
                        assert c == 1, ("levelled cut: generic point with %d preimages" % c, alpha, e)
                    lo, hi = pre_range[s]
                    pre_range[s] = [c if lo is None else min(lo, c), c if hi is None else max(hi, c)]
            # (e) Jacobian sign cross-check
            cms = S.cellmaps if k <= 4 else rng.sample(S.cellmaps, min(25, len(S.cellmaps)))
            for cm in cms:
                fd = jacobian_sign_fd(cm, alpha)
                assert fd == cm.sign, ("Jacobian sign mismatch", s, fmt_partition(alpha), cm.chain, cm.levels, cm.sign, fd)
                n_fd += 1
        results[(k, alpha)] = Ms
        # ---- C11: the comparison ----
        if not matrices_equal(Ms["lev"], Ms["tr"]):
            print("\nDISAGREEMENT at k = %d, alpha = %s (profile %s)" % (k, fmt_partition(alpha), info["profile"]))
            print("levelled (box-product) structure:\n" + fmt_matrix(Ms["lev"]))
            print("Ching's structure (transported through Theta):\n" + fmt_matrix(Ms["tr"]))
            raise AssertionError("C11 FAILS: the two structures induce different maps on homology")
        if k <= 4:
            print("  alpha = %-8s profile %-12s  c^lev_* = c^tr_* =  (rows: Kunneth basis %s; columns: cycle basis of P(%d))"
                  % (fmt_partition(alpha), str(info["profile"]), str(target_basis(cxs, info["profile"])), k))
            print(fmt_matrix(Ms["lev"], indent="      "))
    # ---- (a) unit profiles ----
    for which, alpha in (("(1;%d)" % k, min_partition(cc.T)), ("(%d;1..1)" % k, max_partition(cc.T))):
        for s in STRUCTS:
            unit_signs[(k, which, s)] = check_unit(results[(k, alpha)][s], k, s, which)
    print("  (a) unit profiles: c_* = s * identity with s = %s"
          % {which: [unit_signs[(k, which, s)] for s in STRUCTS] for which in ("(1;%d)" % k, "(%d;1..1)" % k)})
    # ---- (b) equivariance ----
    perms = perm_sample(k)
    A_of = {sigma: cc.top_action_matrix(sigma) for sigma in perms}
    n_hom, n_pt = 0, 0
    for alpha in parts:
        for sigma in perms:
            A = A_of[sigma]
            T, salpha = target_action_matrix(alpha, sigma, cxs)
            for s in STRUCTS:
                lhs = mat_mul(results[(k, salpha)][s], A)
                rhs = mat_mul(T, results[(k, alpha)][s])
                assert matrices_equal(lhs, rhs), ("equivariance on homology fails", k, s, fmt_partition(alpha), cycle_notation(sigma), lhs, rhs)
                n_hom += 1
        for _ in range(6):
            sigma = rng.choice(perms)
            chain, u = random_top_point(k)
            salpha = act_partition(sigma, alpha)
            for s in STRUCTS:
                lhs = FORWARD[s](act_chain(sigma, chain), u, salpha)
                rhs = act_target_point(FORWARD[s](chain, u, alpha), alpha, sigma)
                assert lhs == rhs, ("pointwise equivariance fails", k, s, fmt_partition(alpha), cycle_notation(sigma), chain, u, lhs, rhs)
                n_pt += 1
    print("  (b) equivariance: %d homology identities (%d permutations), %d pointwise identities" % (n_hom, len(perms), n_pt))
    print("  (d) preimage counts per generic point: lev %s, tr %s;  (e) %d Jacobian signs cross-checked" % (pre_range["lev"], pre_range["tr"], n_fd))
    print("  C11 at k = %d: c^lev_* = c^tr_* for all %d partitions   [%.1f s]" % (k, len(parts), time.time() - tk))

print("\nunit signs:", sorted(set(unit_signs.values())))
print("totals: %d target cells, %d affine pieces; runtime %.1f s" % (stats.get("cells", 0), stats.get("cellmaps", 0), time.time() - T0))
print("ALL C11 CHECKS PASSED")
