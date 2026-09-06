"""Check C12: the maps induced on top homology by the cocompositions of C11 assemble into an operad on the dual
groups O(k) := L(k)^vee = H~^{k-1}(P(k); Z), L(k) := H~_{k-1}(P(k); Z) (free of rank (k-1)!, chains.py), and this
operad is the sign-twisted, degree-shifted Lie operad in every arity <= KMAX (default 5, about 1 s;
`python checks_homology_C12.py [kmax]`; kmax = 6 takes about 16 s, almost all of it recomputing the C11 matrices).

CONVENTIONS (fixed before any check; the script reports what it finds under them and adjusts nothing)
  (D1) L(k) has the integral cycle basis K_0 .. K_{(k-1)!-1} of chains.py and the LEFT action sigma.K_j =
       sum_i A(sigma)[i][j] K_i (A = PChainComplex.top_action_matrix).  O(k) = L(k)^vee carries the dual basis
       K_0^vee, ...; an element of O(k) is its column of coordinates.  O(k) is a RIGHT Sigma_k-module through
       (phi.sigma)(x) := phi(sigma.x), i.e. phi.sigma = A(sigma)^T phi  ((phi.sigma).tau = phi.(sigma tau)).
  (D2) O(k) sits in cohomological degree k-1 (it is H~^{k-1}(P(k)); in homological terms degree 1-k, the degree
       of H_*(partial_k Id) = H~^{-*}(P(k))).  A tensor product O(n) (x) O(k_1) (x) ... (x) O(k_n) is indexed by the
       dual Kunneth basis in the lexicographic order of degree.target_basis, the dual of a product basis vector
       being the product of the dual basis vectors (no sign in (L(n) (x) ...)^vee = O(n) (x) ...).  The ONLY source
       of signs is the reordering of tensor factors: moving a factor of arity p past a factor of arity q costs the
       Koszul sign (-1)^{(p-1)(q-1)} of the degrees in (D2); so two factors of even arity anticommute and every
       other pair commutes.  This is the sign (O5) built into C11 (b), and it is recomputed here independently.
  (D3) For a partition alpha = {T_1 < ... < T_n} of {1..k} (blocks ordered by least element, |T_i| = k_i) the
       composition map is gamma_alpha := M_alpha^T : O(n) (x) O(k_1) (x) ... (x) O(k_n) -> O(k), the transpose of
       the C11 matrix M_alpha of (c_alpha)_* (degree.induced_matrix; both structures of C11 are recomputed and
       re-asserted equal, so either may be used).  A block of size 1 contributes O(1) = Z.1, 1 = the dual of the
       point class of P(1); 1 is the unit of the operad.  The composition along a tree is the iteration of these
       maps; the blocks of every partition are always taken in their least-element order.

CHECKS (each ends in an assertion; the actual matrices are printed at a failure)
  (3) units: gamma_{(1;k)}(1 (x) x) = x and gamma_{(k;1..1)}(x (x) 1 (x) ... (x) 1) = x for every basis vector x
      of O(k), k <= KMAX (the transposes of C11 (a), re-read on the dual side, and gamma(1 (x) 1) = 1 for k = 1).
  (1) equivariance: for every partition alpha of {1..k} and every sigma in Sigma_k (all of Sigma_k for k <= 5,
      generators and a random sample for k = 6):  gamma_{sigma.alpha}(y).sigma = gamma_alpha(y.sigma)  for every
      dual Kunneth basis vector y of the target of sigma.alpha, where y.sigma := T(alpha, sigma)^T y is
      koszul * (y_0.sbar) (x) (the factor at the block sigma(T_i) moved to the slot of T_i and acted on by the
      induced permutation sigma_i of {1..k_i}), sbar the block permutation and koszul the sign (D2) of the
      reordering (recomputed here from sbar and the arities and compared with degree.induced_block_data).
      In matrices: A(sigma)^T M_{sigma.alpha}^T = M_alpha^T T(alpha, sigma)^T.
  (2) associativity (coassociativity of the cocompositions): for every chain alpha <= beta of partitions of
      {1..k} (alpha refines beta; beta has blocks B_1 < .. < B_m, alpha has T_1 < .. < T_n, C_j = {i : T_i in B_j},
      rho = beta/alpha the induced partition of {1..n}) the two composites
        route 1:  L(k) --c_beta--> L(m) (x) (x)_j L(|B_j|) --(x)_j c_{alpha|B_j}--> L(m) (x) (x)_j (L(n_j) (x) (x)_{i in C_j} L(k_i))
        route 2:  L(k) --c_alpha--> L(n) (x) (x)_i L(k_i) --c_rho (x) id--> L(m) (x) (x)_j L(n_j) (x) (x)_i L(k_i)
      agree after the factors of route 1 are brought into the order of route 2 with the Koszul sign (D2); alpha|B_j
      is alpha restricted to B_j and relabelled through the order-preserving bijection B_j = {1..|B_j|}.  The
      transpose of this identity is the associativity axiom gamma(gamma(a; b_1..b_m); c_1..c_n) =
      +-gamma(a; gamma(b_1; c_{C_1}), ..., gamma(b_m; c_{C_m})) of the operad O in its partition-indexed form.
  (4) the Lie operad.  e := K_0^vee, the generator of O(2) = Z.
      (i)   e.(12) = ? e (the Sigma_2-action on L(2)); reported.
      (ii)  the three bracketings A = [[x1,x2],x3] = gamma_{12|3}(e (x) e (x) 1), B = [[x1,x3],x2] =
            gamma_{13|2}(e (x) e (x) 1), C = [x1,[x2,x3]] = gamma_{1|23}(e (x) 1 (x) e) in O(3) = Z^2 satisfy
            exactly one primitive integer relation, which is computed and reported (the Jacobi identity).
      (iii) the left-normed brackets [[..[x_1, x_{s(2)}], x_{s(3)}] .., x_{s(k)}], s(1) = 1, built by iterating
            gamma_{T|{s(k)}}(e (x) - (x) 1) (T = {1, s(2), .., s(k-1)}) from the root down, form a Z-basis of O(k):
            the determinant of their (k-1)! x (k-1)! coordinate matrix is +-1 (reported).  The same brackets are
            rebuilt leaf-first, as gamma_{s(1)s(2)|singletons}(z (x) e (x) 1 .. (x) 1) with z a left-normed bracket
            of arity k-1 at the root, and must agree (an explicit instance of (2)).
Consequence (arities <= KMAX): O is a graded operad; e is symmetric of degree 1 and satisfies the Jacobi
relation of (ii); the free operad on such an e modulo that relation is the operadic suspension of Lie, whose
arity-k part is free of rank (k-1)! with the left-normed brackets as a basis, so its map to O is a surjection
between free Z-modules of equal rank, hence an isomorphism of operads in arities <= KMAX.
Ends with ALL C12 CHECKS PASSED or raises at the first failure."""
import os
import sys
import time
import random
from fractions import Fraction
from itertools import permutations, product
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
from degree import (Complexes, all_partitions, induced_matrix, target_basis, target_action_matrix,     # noqa: E402
                    induced_block_data, blocks_of, order_iso, profile_of, matrices_equal, fmt_matrix, fmt_partition)
from chains import all_perms, generators, cycle_notation  # noqa: E402
from partitions import min_partition, max_partition  # noqa: E402

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 5
STRUCTS = ("lev", "tr")
rng = random.Random(2026)
cxs = Complexes()
T0 = time.time()


# ---------------------------------------------------------------------------------------------------------
# small exact linear algebra: sparse matrices as lists of row dicts, dense matrices as lists of lists
# ---------------------------------------------------------------------------------------------------------

def sp_of(M):
    return [{j: v for j, v in enumerate(row) if v} for row in M]


def sp_transpose(A, ncols):
    T = [dict() for _ in range(ncols)]
    for i, row in enumerate(A):
        for j, v in row.items():
            T[j][i] = v
    return T


def sp_mul(A, B):
    out = []
    for row in A:
        acc = {}
        for j, a in row.items():
            for l, b in B[j].items():
                acc[l] = acc.get(l, 0) + a * b
        out.append({l: v for l, v in acc.items() if v})
    return out


def sp_dense(A, ncols):
    return [[row.get(j, 0) for j in range(ncols)] for row in A]


def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def det_bareiss(M):
    """Exact integer determinant (fraction-free Gaussian elimination)."""
    A = [list(r) for r in M]
    n = len(A)
    if n == 0:
        return 1
    sgn, prev = 1, 1
    for c in range(n - 1):
        if A[c][c] == 0:
            p = next((i for i in range(c + 1, n) if A[i][c] != 0), None)
            if p is None:
                return 0
            A[c], A[p] = A[p], A[c]
            sgn = -sgn
        for i in range(c + 1, n):
            for j in range(c + 1, n):
                A[i][j] = (A[i][j] * A[c][c] - A[i][c] * A[c][j]) // prev
        prev = A[c][c]
    return sgn * A[n - 1][n - 1]


def det_fraction(M):
    A = [[Fraction(x) for x in r] for r in M]
    n = len(A)
    d = Fraction(1)
    for c in range(n):
        p = next((i for i in range(c, n) if A[i][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]
            if f:
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return d


def kernel_Q(cols):
    """Basis of the kernel of the matrix with the given dense columns (Fractions), as coefficient vectors."""
    m = len(cols)
    n = len(cols[0])
    A = [[Fraction(cols[j][i]) for j in range(m)] for i in range(n)]
    piv, r = [], 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        A[r] = [x / A[r][c] for x in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(m) if c not in piv]
    ker = []
    for f in free:
        v = [Fraction(0)] * m
        v[f] = Fraction(1)
        for i, c in enumerate(piv):
            v[c] = -A[i][f]
        ker.append(v)
    return ker


def koszul_sign(degrees, positions):
    """Sign of reordering graded factors: factor s (degree degrees[s]) goes to position positions[s];
    prod over pairs s < t with positions[s] > positions[t] of (-1)^{degrees[s] degrees[t]}."""
    sgn = 1
    for s in range(len(degrees)):
        for t in range(s + 1, len(degrees)):
            if positions[s] > positions[t] and degrees[s] % 2 and degrees[t] % 2:
                sgn = -sgn
    return sgn


# ---------------------------------------------------------------------------------------------------------
# the C11 matrices (both structures, re-asserted equal), the dual bases and the composition maps
# ---------------------------------------------------------------------------------------------------------
M = {}            # partition -> dense matrix of (c_alpha)_* (rows: Kunneth basis, columns: cycle basis)
ROWS = {}         # partition -> list of Kunneth index tuples (row order)
ROWPOS = {}       # partition -> dict tuple -> row
MCOL = {}         # partition -> list over source columns of dict {row tuple: coefficient}
PARTS = {}        # k -> partitions of {1..k}
rank = {}         # k -> (k-1)!

print("computing the C11 matrices for k <= %d (both structures) ..." % KMAX)
for k in range(1, KMAX + 1):
    tk = time.time()
    cc = cxs(k)
    rank[k] = len(cc._top[1])
    PARTS[k] = all_partitions(k)
    for alpha in PARTS[k]:
        Ms = {}
        for s in STRUCTS:
            Ms[s], info = induced_matrix(s, k, alpha, cxs, rng)
        if not matrices_equal(Ms["lev"], Ms["tr"]):
            print("levelled:\n" + fmt_matrix(Ms["lev"]) + "\nChing:\n" + fmt_matrix(Ms["tr"]))
            raise AssertionError("C11 input fails at alpha = %s" % fmt_partition(alpha))
        M[alpha] = Ms["lev"]
        ROWS[alpha] = target_basis(cxs, info["profile"])
        ROWPOS[alpha] = {ix: r for r, ix in enumerate(ROWS[alpha])}
        assert len(ROWS[alpha]) == len(M[alpha]) and all(len(row) == rank[k] for row in M[alpha])
        MCOL[alpha] = [{ROWS[alpha][r]: M[alpha][r][j] for r in range(len(M[alpha])) if M[alpha][r][j]}
                       for j in range(rank[k])]
    print("  k = %d: %d partitions, rank O(%d) = %d   [%.1f s]" % (k, len(PARTS[k]), k, rank[k], time.time() - tk))


def unit_partition(k):
    return min_partition(frozenset(range(1, k + 1)))          # profile (1; k): one block


def discrete_partition(k):
    return max_partition(frozenset(range(1, k + 1)))          # profile (k; 1..1): all singletons


ONE = [1]                                                     # the unit 1 in O(1)


def gamma(alpha, ys):
    """gamma_alpha(y_0 (x) y_1 (x) ... (x) y_n) in O(k) (dense coordinate lists) = M_alpha^T applied to the
    coefficient vector of the tensor in the dual Kunneth basis (no signs, (D2))."""
    blocks = blocks_of(alpha)
    assert len(ys) == len(blocks) + 1
    assert len(ys[0]) == rank[len(blocks)] and all(len(y) == rank[len(b)] for y, b in zip(ys[1:], blocks))
    k = sum(len(b) for b in blocks)
    out = [0] * rank[k]
    Ma = M[alpha]
    for r, ix in enumerate(ROWS[alpha]):
        c = ys[0][ix[0]]
        for i in range(1, len(ix)):
            if not c:
                break
            c *= ys[i][ix[i]]
        if c:
            row = Ma[r]
            for j in range(rank[k]):
                if row[j]:
                    out[j] += c * row[j]
    return out


# cached action matrices A(sigma) (left action on L(k)) in sparse column form and A(sigma)^T in row form
_A = {}


def A_of(k, sigma):
    if (k, sigma) not in _A:
        Ad = cxs(k).top_action_matrix(sigma, verify=(k <= 5))
        _A[(k, sigma)] = (Ad, sp_of(Ad), sp_transpose(sp_of(Ad), rank[k]))
    return _A[(k, sigma)]


def act_dual(k, sigma, phi):
    """phi.sigma = A(sigma)^T phi for a dense coordinate list phi of O(k)."""
    Ad = A_of(k, sigma)[0]
    return [sum(Ad[i][j] * phi[i] for i in range(rank[k])) for j in range(rank[k])]


def dual_action_T(alpha, sigma):
    """T(alpha, sigma)^T as a sparse matrix (rows: Kunneth basis of the alpha-target, columns: of the
    sigma.alpha-target) from cached action matrices; identical in content to degree.target_action_matrix.
    Also returns the block permutation, the induced permutations and the Koszul sign."""
    salpha, sbar, sigmas, koszul = induced_block_data(alpha, sigma)
    blocks = blocks_of(alpha)
    n = len(blocks)
    prof = profile_of(alpha)
    # the Koszul sign of (D2), recomputed from sbar and the arities, must agree with degree.py
    assert koszul == koszul_sign([len(b) - 1 for b in blocks], sbar), (fmt_partition(alpha), sigma, koszul)
    An = A_of(n, tuple(s + 1 for s in sbar))[1]
    Ai = [A_of(prof[i], sigmas[i])[1] for i in range(n)]
    TT = []
    spos = ROWPOS[salpha]
    for ix in ROWS[alpha]:
        row = {}
        # nonzero entries of column ix[0] of A_n(sbar): rows i0 with An[i0][ix[0]] != 0
        fac0 = [(i0, r.get(ix[0], 0)) for i0, r in enumerate(An) if r.get(ix[0], 0)]
        facs = [[(ii, r.get(ix[i + 1], 0)) for ii, r in enumerate(Ai[i]) if r.get(ix[i + 1], 0)] for i in range(n)]
        for choice in product(fac0, *facs):
            ix2 = [None] * (n + 1)
            ix2[0] = choice[0][0]
            val = koszul * choice[0][1]
            for i in range(n):
                ix2[sbar[i] + 1] = choice[i + 1][0]
                val *= choice[i + 1][1]
            row[spos[tuple(ix2)]] = row.get(spos[tuple(ix2)], 0) + val
        TT.append({c: v for c, v in row.items() if v})
    return TT, salpha, sbar, sigmas, koszul


def perm_sample(k):
    if k <= 5:
        return all_perms(k)
    ps = list(generators(k))
    for _ in range(10):
        p = list(range(1, k + 1))
        rng.shuffle(p)
        ps.append(tuple(p))
    return ps


# ---------------------------------------------------------------------------------------------------------
# (3) units
# ---------------------------------------------------------------------------------------------------------
print("\n(3) units")
n_unit = 0
for k in range(1, KMAX + 1):
    for j in range(rank[k]):
        x = [1 if i == j else 0 for i in range(rank[k])]
        lhs = gamma(unit_partition(k), [ONE, x])
        assert lhs == x, ("unit axiom (1;k) fails", k, j, lhs)
        rhs = gamma(discrete_partition(k), [x] + [ONE] * k)
        assert rhs == x, ("unit axiom (k;1..1) fails", k, j, rhs)
        n_unit += 2
    assert matrices_equal(M[unit_partition(k)], identity(rank[k])) and matrices_equal(M[discrete_partition(k)], identity(rank[k]))
print("  gamma_(1;k)(1 (x) x) = x and gamma_(k;1..1)(x (x) 1 .. (x) 1) = x on every basis vector, k <= %d: %d identities"
      % (KMAX, n_unit))

# ---------------------------------------------------------------------------------------------------------
# (1) equivariance on the dual side
# ---------------------------------------------------------------------------------------------------------
print("\n(1) equivariance: gamma_{sigma.alpha}(y).sigma = gamma_alpha(y.sigma), phi.sigma = A(sigma)^T phi, Koszul sign (D2)")
MT = {alpha: sp_transpose(sp_of(M[alpha]), rank[sum(len(b) for b in alpha)]) for alpha in M}
n_eq, n_neg, n_cross, n_forced = 0, 0, 0, 0
example = None
for k in range(2, KMAX + 1):
    tk = time.time()
    perms = perm_sample(k)
    neg_k = 0
    for alpha in PARTS[k]:
        for sigma in perms:
            TT, salpha, sbar, sigmas, koszul = dual_action_T(alpha, sigma)
            AT = A_of(k, sigma)[2]
            lhs = sp_mul(AT, MT[salpha])            # A(sigma)^T M_{sigma.alpha}^T
            rhs = sp_mul(MT[alpha], TT)             # M_alpha^T T(alpha, sigma)^T
            if lhs != rhs:
                R = len(ROWS[salpha])
                print("alpha = %s, sigma = %s, sigma.alpha = %s, koszul = %d" % (fmt_partition(alpha), cycle_notation(sigma), fmt_partition(salpha), koszul))
                print("A(sigma)^T M_{sigma.alpha}^T =\n" + fmt_matrix(sp_dense(lhs, R)))
                print("M_alpha^T T(alpha,sigma)^T =\n" + fmt_matrix(sp_dense(rhs, R)))
                raise AssertionError("C12 (1) FAILS: equivariance on the dual side")
            n_eq += 1
            if koszul == -1:
                neg_k += 1
                if example is None:
                    example = (k, alpha, sigma, salpha, sbar, sigmas)
                # the sign is forced, not vacuous: with the opposite sign the identity must fail whenever
                # the two sides are nonzero (rhs is linear in the global sign)
                if any(row for row in lhs):
                    assert lhs != [{c: -v for c, v in row.items()} for row in rhs], "Koszul sign not forced"
                    n_forced += 1
    n_neg += neg_k
    # cross-check the cached T against degree.target_action_matrix on a few pairs
    for _ in range(6 if k <= 5 else 3):
        alpha, sigma = rng.choice(PARTS[k]), rng.choice(perms)
        Tdeg, sa = target_action_matrix(alpha, sigma, cxs)
        TT, salpha, _, _, _ = dual_action_T(alpha, sigma)
        assert sa == salpha and matrices_equal(sp_dense(sp_transpose(TT, len(ROWS[salpha])), len(ROWS[alpha])), Tdeg)
        n_cross += 1
    print("  k = %d: %d identities (%d partitions x %d permutations%s), Koszul sign -1 in %d of them   [%.1f s]"
          % (k, len(PARTS[k]) * len(perms), len(PARTS[k]), len(perms),
             "" if k <= 5 else " = generators + sample", neg_k, time.time() - tk))
k, alpha, sigma, salpha, sbar, sigmas = example
print("  first identity with Koszul sign -1: alpha = %s, sigma = %s, sigma.alpha = %s, block permutation %s, "
      "within-block permutations %s: y.sigma = -(y_0.sbar) (x) (reordered factors)"
      % (fmt_partition(alpha), cycle_notation(sigma), fmt_partition(salpha), tuple(s + 1 for s in sbar), sigmas))
print("  %d identities in all, %d with Koszul sign -1 (in %d of these both sides are nonzero and the identity fails with the sign +1);"
      " %d T matrices cross-checked against degree.target_action_matrix" % (n_eq, n_neg, n_forced, n_cross))

# ---------------------------------------------------------------------------------------------------------
# (2) associativity = coassociativity of the cocompositions along chains alpha <= beta
# ---------------------------------------------------------------------------------------------------------
print("\n(2) associativity: route 1 (c_beta, then c_{alpha|B_j}) = Koszul reordering of route 2 (c_alpha, then c_{beta/alpha})")


def restrict_partition(alpha, B):
    iota = order_iso(B)
    return frozenset(frozenset(iota[x] for x in T) for T in alpha if T <= B)


def check_assoc(k, alpha, rho):
    """The identity for the chain alpha <= beta, beta = union of the blocks of alpha along rho (a partition of
    {1..n}, n = number of blocks of alpha).  Returns (beta, koszul sign of the reordering, nontrivial?)."""
    blocks = blocks_of(alpha)
    n = len(blocks)
    Cs = blocks_of(rho)                                   # C_j, ordered by least block index
    m = len(Cs)
    Bs = [frozenset().union(*[blocks[i - 1] for i in C]) for C in Cs]
    beta = frozenset(Bs)
    assert blocks_of(beta) == Bs, "block orders of beta and beta/alpha disagree"
    restr = [restrict_partition(alpha, B) for B in Bs]
    for j, B in enumerate(Bs):                            # blocks of alpha|B_j <-> T_i, i in C_j, in order
        assert profile_of(restr[j]) == tuple(len(blocks[i - 1]) for i in sorted(Cs[j]))
    n_j = [len(C) for C in Cs]
    # factor labels in the two orders and the Koszul sign of the reordering
    labels1 = [("m",)]
    for j in range(m):
        labels1.append(("n", j))
        labels1.extend(("k", i) for i in sorted(Cs[j]))
    labels2 = [("m",)] + [("n", j) for j in range(m)] + [("k", i) for i in range(1, n + 1)]
    deg = {("m",): m - 1}
    deg.update({("n", j): n_j[j] - 1 for j in range(m)})
    deg.update({("k", i): len(blocks[i - 1]) - 1 for i in range(1, n + 1)})
    pos2 = {lab: p for p, lab in enumerate(labels2)}
    sgn = koszul_sign([deg[lab] for lab in labels1], [pos2[lab] for lab in labels1])
    pos1 = {lab: p for p, lab in enumerate(labels1)}
    reorder = [pos1[lab] for lab in labels2]              # key2[t] = key1[reorder[t]]
    for a in range(rank[k]):
        # route 2
        out2 = {}
        for ix, x in MCOL[alpha][a].items():
            for ix2, y in MCOL[rho][ix[0]].items():
                key = ix2 + ix[1:]
                out2[key] = out2.get(key, 0) + x * y
        out2 = {key: v for key, v in out2.items() if v}
        # route 1
        out1 = {}
        for ix, x in MCOL[beta][a].items():
            terms = [((ix[0],), x)]
            for j in range(m):
                terms = [(key + ixj, c * y) for key, c in terms for ixj, y in MCOL[restr[j]][ix[j + 1]].items()]
            for key, c in terms:
                key2 = tuple(key[t] for t in reorder)
                out1[key2] = out1.get(key2, 0) + sgn * c
        out1 = {key: v for key, v in out1.items() if v}
        if out1 != out2:
            print("alpha = %s <= beta = %s (rho = %s), source cycle %d, Koszul sign %d" % (fmt_partition(alpha), fmt_partition(beta), fmt_partition(rho), a, sgn))
            print("route 1 (reordered):", sorted(out1.items()))
            print("route 2:", sorted(out2.items()))
            raise AssertionError("C12 (2) FAILS: coassociativity")
        if sgn == -1 and out1:                              # the reordering sign is forced by the data
            assert {key: -v for key, v in out1.items()} != out2, "Koszul sign of the reordering not forced"
            forced[0] = True
    nontrivial = (len(Cs) not in (1, n)) and n != k
    return beta, sgn, nontrivial


forced = [False]


n_assoc, n_nontriv, n_signed, n_sforced = 0, 0, 0, 0
for k in range(1, KMAX + 1):
    tk = time.time()
    ck, cn, cs, cf = 0, 0, 0, 0
    for alpha in PARTS[k]:
        n = len(alpha)
        for rho in PARTS[n]:
            forced[0] = False
            beta, sgn, nontrivial = check_assoc(k, alpha, rho)
            ck += 1
            cn += nontrivial
            cs += (sgn == -1)
            cf += forced[0]
    n_assoc += ck
    n_nontriv += cn
    n_signed += cs
    n_sforced += cf
    print("  k = %d: %d chains alpha <= beta (%d with alpha, beta both proper and distinct), reordering sign -1 in %d (forced in %d)   [%.1f s]"
          % (k, ck, cn, cs, cf, time.time() - tk))
print("  %d coassociativity identities (%d nontrivial, %d with Koszul sign -1, of which %d nonzero and failing with the sign +1)"
      % (n_assoc, n_nontriv, n_signed, n_sforced))

# ---------------------------------------------------------------------------------------------------------
# (4) the Lie operad
# ---------------------------------------------------------------------------------------------------------
print("\n(4) the Lie operad")
e = [1]
A2 = A_of(2, (2, 1))[0]
e_swapped = act_dual(2, (2, 1), e)
print("  (i) A((12)) on L(2) = %s, so e.(12) = %s e: e is Sigma_2-%s" % (A2, "+" if e_swapped == e else "-", "invariant" if e_swapped == e else "anti-invariant"))
assert e_swapped in (e, [-1])
sym_sign = 1 if e_swapped == e else -1


def part(*blocks):
    return frozenset(frozenset(b) for b in blocks)


A3 = gamma(part((1, 2), (3,)), [e, e, ONE])
B3 = gamma(part((1, 3), (2,)), [e, e, ONE])
C3 = gamma(part((1,), (2, 3)), [e, ONE, e])
print("  (ii) in the dual basis of O(3): [[x1,x2],x3] = %s, [[x1,x3],x2] = %s, [x1,[x2,x3]] = %s" % (A3, B3, C3))
ker = kernel_Q([A3, B3, C3])
assert len(ker) == 1, ("expected exactly one relation among the three bracketings", ker)
v = ker[0]
den = 1
for x in v:
    den = den * x.denominator // gcd(den, x.denominator)
rel = [int(x * den) for x in v]
g = 0
for x in rel:
    g = gcd(g, abs(x))
rel = [x // g for x in rel]
if rel[0] < 0:
    rel = [-x for x in rel]
assert all(abs(x) == 1 for x in rel), ("Jacobi relation has a coefficient != +-1", rel)
assert all(sum(rel[t] * vec[j] for t, vec in enumerate((A3, B3, C3))) == 0 for j in range(2))


def signed(c, name):
    return ("+ " if c > 0 else "- ") + name


print("  Jacobi identity in O(3):  %s %s %s = 0" % (signed(rel[0], "[[x1,x2],x3]"), signed(rel[1], "[[x1,x3],x2]"), signed(rel[2], "[x1,[x2,x3]]")))
jacobi = tuple(rel)


def model_jacobi_patterns():
    """Independent model for the expected signs: g = gl(1|1) as a Z-graded Lie algebra (E11, E22 in degree 0, E12,
    E21 in degree 1, [X,Y] = XY - (-1)^{|X||Y|} YX), V = g[1] (|sx| = |x| - 1, cohomologically) with the
    graded-symmetric degree-1 operation e(sx, sy) = (-1)^{|x|} s[x, y].  The three trees act on V^{(x)3} through the
    Koszul conventions of End_V: A(v1,v2,v3) = e(e(v1,v2),v3), B = (-1)^{|v2||v3|} e(e(v1,v3),v2) (shuffle of the
    inputs into the block order of 13|2), C = (-1)^{|v1|} e(v1, e(v2,v3)) (the operation e of 1 (x) e passes v1).
    Returns the sign patterns (sA, sB, sC) with sA A + sB B + sC C = 0 on every basis triple; the operadic
    suspension of Lie predicts exactly one."""
    basis = ["E11", "E22", "E12", "E21"]
    deg = {"E11": 0, "E22": 0, "E12": 1, "E21": 1}

    def mat(b):
        i, j = int(b[1]) - 1, int(b[2]) - 1
        Mm = [[0, 0], [0, 0]]
        Mm[i][j] = 1
        return Mm

    def mul(X, Y):
        return [[sum(X[i][l] * Y[l][j] for l in range(2)) for j in range(2)] for i in range(2)]

    def coords(Mm):
        return {b: Mm[int(b[1]) - 1][int(b[2]) - 1] for b in basis if Mm[int(b[1]) - 1][int(b[2]) - 1]}

    def gdeg(x):
        ds = {deg[b] for b in x}
        assert len(ds) == 1
        return ds.pop()

    def bracket(x, y):
        out = {}
        for bx, cx in x.items():
            for by, cy in y.items():
                s = (-1) ** (deg[bx] * deg[by])
                XY, YX = mul(mat(bx), mat(by)), mul(mat(by), mat(bx))
                for b, c in coords([[XY[i][j] - s * YX[i][j] for j in range(2)] for i in range(2)]).items():
                    out[b] = out.get(b, 0) + cx * cy * c
        return {b: c for b, c in out.items() if c}

    def ee(x, y):                                   # e(sx, sy) = (-1)^{|x|} s[x, y], as an element of g
        if not x or not y:
            return {}
        return {b: (-1) ** gdeg(x) * c for b, c in bracket(x, y).items()}

    def vdeg(x):
        return gdeg(x) - 1

    def scale(v, s):
        return {b: s * c for b, c in v.items()}

    def add(*vs):
        out = {}
        for v in vs:
            for b, c in v.items():
                out[b] = out.get(b, 0) + c
        return {b: c for b, c in out.items() if c}

    assert all(ee({by: 1}, {bx: 1}) == scale(ee({bx: 1}, {by: 1}), (-1) ** (vdeg({bx: 1}) * vdeg({by: 1})))
               for bx in basis for by in basis), "e is not graded-symmetric on gl(1|1)[1]"
    triples = []
    for bx, by, bz in product(basis, repeat=3):
        x, y, z = {bx: 1}, {by: 1}, {bz: 1}
        A_ = ee(ee(x, y), z)
        B_ = scale(ee(ee(x, z), y), (-1) ** (vdeg(y) * vdeg(z)))
        C_ = scale(ee(x, ee(y, z)), (-1) ** vdeg(x))
        triples.append((A_, B_, C_))
    assert any(A_ or B_ or C_ for A_, B_, C_ in triples)
    patterns = []
    for sB, sC in product((1, -1), repeat=2):
        if all(not add(A_, scale(B_, sB), scale(C_, sC)) for A_, B_, C_ in triples):
            patterns.append((1, sB, sC))
    return patterns


patterns = model_jacobi_patterns()
assert patterns == [jacobi], ("the Jacobi signs of O(3) are not those of End(g[1]), g = gl(1|1)", jacobi, patterns)
print("  model check: in End(g[1]) for the graded Lie algebra g = gl(1|1), with the same Koszul conventions, exactly the sign pattern %s holds" % (jacobi,))


def lnb(seq):
    """[[..[x_{seq[0]}, x_{seq[1]}], ..], x_{seq[-1]}] in O(len(seq)), seq a permutation of {1..k} with seq[0] = 1,
    built from the root down: gamma_{T | {seq[-1]}}(e (x) lnb(seq[:-1] relabelled) (x) 1), T = {seq[0..k-2]}."""
    k = len(seq)
    assert seq[0] == 1 and sorted(seq) == list(range(1, k + 1))
    if k == 1:
        return ONE
    if k == 2:
        return e
    T = frozenset(seq[:-1])
    iota = order_iso(T)
    y = lnb(tuple(iota[x] for x in seq[:-1]))
    return gamma(frozenset([T, frozenset([seq[-1]])]), [e, y, ONE])


def lnb_leaf_first(seq):
    """The same bracket built from the leaves up: z (x) e (x) 1 ... along alpha = {seq[0] seq[1] | singletons},
    z the left-normed bracket of arity k-1 whose first input is the inner bracket and whose input number
    iota({x}) is the leaf x."""
    k = len(seq)
    if k <= 2:
        return lnb(seq)
    blocks = [frozenset(seq[:2])] + [frozenset([x]) for x in seq[2:]]
    alpha = frozenset(blocks)
    ordered = blocks_of(alpha)
    assert ordered[0] == blocks[0]
    slot = {b: i + 1 for i, b in enumerate(ordered)}      # block -> input number of z
    z = lnb((1,) + tuple(slot[frozenset([x])] for x in seq[2:]))
    ys = [z] + [e if b == blocks[0] else ONE for b in ordered]
    return gamma(alpha, ys)


dets = {}
for k in range(2, KMAX + 1):
    tk = time.time()
    seqs = [(1,) + tail for tail in permutations(range(2, k + 1))]
    rows = [lnb(s) for s in seqs]
    for s, r in zip(seqs, rows):
        r2 = lnb_leaf_first(s)
        assert r2 == r, ("root-first and leaf-first constructions of the left-normed bracket disagree", s, r, r2)
    d = det_bareiss(rows)
    assert det_fraction(rows) == d                       # second opinion on the determinant
    dets[k] = d
    if k <= 4:
        print("  (iii) k = %d: left-normed brackets (rows, in the order %s) in the dual basis of O(%d):" % (k, seqs, k))
        print(fmt_matrix(rows, indent="        "))
    assert abs(d) == 1, ("the left-normed brackets are not a Z-basis of O(%d): determinant %d" % (k, d), rows)
    print("  (iii) k = %d: %d left-normed brackets [[..[x_1,x_s(2)],..],x_s(k)] (root-first = leaf-first), determinant %+d: a Z-basis of O(%d) = Z^%d   [%.1f s]"
          % (k, len(seqs), d, k, rank[k], time.time() - tk))

# ---------------------------------------------------------------------------------------------------------
print("\nsummary")
print("  conventions: right action phi.sigma = A(sigma)^T phi; O(k) in cohomological degree k-1; Koszul sign (-1)^{(p-1)(q-1)} for")
print("  moving a factor of arity p past one of arity q; gamma_alpha = M_alpha^T, blocks in least-element order, no other signs.")
print("  units: %d identities; equivariance: %d identities (%d with Koszul sign -1); associativity: %d chains alpha <= beta (%d nontrivial)."
      % (n_unit, n_eq, n_neg, n_assoc, n_nontriv))
print("  e.(12) = %s e;  Jacobi:  %s [[x1,x2],x3] %s [[x1,x3],x2] %s [x1,[x2,x3]] = 0;  determinants of the left-normed bases: %s"
      % ("+" if sym_sign == 1 else "-", *["+" if c > 0 else "-" for c in jacobi], {k: "%+d" % d for k, d in dets.items()}))
print("  => in arities <= %d the operad on H~^{k-1}(P(k)) is generated by the symmetric degree-1 element e with the Jacobi relation above,"
      "\n     i.e. it is the operadic suspension of Lie (arity k: Lie(k) (x) sgn_k in cohomological degree k-1)." % KMAX)
print("runtime %.1f s" % (time.time() - T0))
print("ALL C12 CHECKS PASSED")
