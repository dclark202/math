"""Reduced normalized chain complex of the pointed simplicial set P(n) (Definition 2.4) over Z, its homology
(ranks over Q by exact Fraction elimination; ranks and torsion over Z by a sparse column-echelon form over Z,
with a dense Smith normal form as cross-check/fallback), an explicit basis of cycles for the top homology
H~_{n-1}(P(n); Q), and the matrices of the Sigma_n-action on it.

Pure Python 3, no third-party packages.  Reuses the chain model of ../partitions.py (chains of partitions
listed coarse -> fine, face d_i = drop entry i, None = basepoint).  Run `python chains.py [nmax]` for a
summary table (default nmax = 6, about 10 s).

Conventions (see README.md):

* C_k has one basis element per nondegenerate non-basepoint k-simplex of P(n), i.e. per strict chain
  min = l_0 < l_1 < ... < l_k = max (only k = 0 for n = 1, where min = max; k = 1..n-1 for n >= 2).
  Basis elements are listed in a fixed canonical order (`ckey`) so that all matrices are reproducible.
* The simplex [l_0 < ... < l_k] is oriented by the order of the chain (the vertex order of the nerve).
  boundary(x) = sum_{i=0}^{k} (-1)^i d_i x, where the face d_i x is replaced by 0 when it is the basepoint
  of P(n) (this happens exactly for i = 0 and i = k when k >= 1: the face no longer starts at min, resp.
  no longer ends at max) or when it is a degenerate simplex (cannot happen for a strict chain; the code
  still tests for it and counts how often each kind of face was dropped).  Concretely
      boundary[l_0 < ... < l_k] = sum_{i=1}^{k-1} (-1)^i [l_0 < ... < l_{i-1} < l_{i+1} < ... < l_k].
  This is the normalized chain complex of P(n) divided by that of the basepoint, so its homology is the
  reduced homology of |P(n)| = Par(n).  There is no degree -1: "reduced" means modulo the basepoint,
  not augmented; in particular H~_0(P(1)) = Z (Par(1) = S^0).
* Sigma_n acts on the LEFT.  A permutation is a tuple p with p[x-1] = sigma(x), (sigma tau)(x) =
  sigma(tau(x)); sigma acts on a partition by sigma.lambda = {sigma(B) : B in lambda}, on a chain
  componentwise, and on C_k by permuting basis elements with coefficient +1 (sigma preserves the order of
  a chain, so no orientation sign appears).  Hence P(sigma tau) = P(sigma) P(tau) on chains, and the same
  on homology (checked numerically in checks_homology_C10.py).
* H~_{n-1} = ker(boundary_{n-1}) (there is no C_n).  A basis is read off a column reduction of
  boundary_{n-1} over Q that records the column operations (`kernel_basis_Q`): the columns that reduce to
  zero ("free" columns f_1 < ... < f_m) give cycles K_1, ..., K_m with K_i[f_i] = 1 and K_i[f_j] = 0 for
  j != i (all other support lies in pivot columns), so the coordinates of any cycle w in this basis are
  simply (w[f_1], ..., w[f_m]) and the matrix of sigma is M(sigma)[i][j] = (sigma.K_j)[f_i].  When every
  pivot of the reduction is +-1 (the case for all n <= 6) the K_i are integer vectors and hence also a
  Z-basis of ker(boundary_{n-1}) (every integer cycle has integer coordinates), and the M(sigma) are
  integer matrices.
* Torsion.  H~_k = ker d_k / im d_{k+1} with C_k free, so tors H~_k = tors(C_k / im d_{k+1}) is read off
  the Smith normal form of d_{k+1}.  `z_echelon` computes a sparse column-echelon form H = d_{k+1} U over Z
  (U unimodular, so im H = im d_{k+1}): its nonzero columns have distinct lowest rows and their number is
  the rank r (over Z and Q).  If every pivot entry is +-1, the r pivot rows and the r nonzero columns form
  a triangular r x r minor of determinant +-1, so the gcd of the r x r minors of H is 1, every invariant
  factor of d_{k+1} is 1, and C_k / im d_{k+1} is free of rank dim C_k - r.  Otherwise (never happens for
  n <= 6) the dense Smith normal form of the nonzero columns of H is computed (`smith_invariants`, which
  the checks also run on the full boundary matrices for n <= 5 as an independent second opinion).
"""
import os
import sys
from fractions import Fraction
from itertools import permutations
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
from partitions import nondegenerate_chains, face, min_partition, max_partition  # noqa: E402


# ----------------------------------------------------------------------------------------------
# chains and permutations
# ----------------------------------------------------------------------------------------------

def ckey(chain):
    """Canonical sort key of a chain (tuple of partitions = frozensets of frozensets)."""
    return tuple(tuple(sorted(tuple(sorted(B)) for B in lam)) for lam in chain)


def is_degenerate(chain):
    return any(chain[i] == chain[i + 1] for i in range(len(chain) - 1))


def all_perms(n):
    return list(permutations(range(1, n + 1)))


def identity(n):
    return tuple(range(1, n + 1))


def compose(s, t):
    """(s t)(x) = s(t(x))."""
    return tuple(s[t[x - 1] - 1] for x in range(1, len(s) + 1))


def inverse(s):
    inv = [0] * len(s)
    for x, y in enumerate(s, 1):
        inv[y - 1] = x
    return tuple(inv)


def cycles(s):
    """Cycles of s (including fixed points), each starting at its least element, ordered by that element."""
    n = len(s)
    seen = [False] * n
    out = []
    for x in range(1, n + 1):
        if not seen[x - 1]:
            cyc, y = [], x
            while not seen[y - 1]:
                seen[y - 1] = True
                cyc.append(y)
                y = s[y - 1]
            out.append(tuple(cyc))
    return out


def cycle_type(s):
    return tuple(sorted((len(c) for c in cycles(s)), reverse=True))


def cycle_notation(s):
    cs = [c for c in cycles(s) if len(c) > 1]
    return "".join("(" + " ".join(map(str, c)) + ")" for c in cs) if cs else "e"


def sign(s):
    return (-1) ** sum(len(c) - 1 for c in cycles(s))


def generators(n):
    """The adjacent transpositions (i i+1) and the n-cycle (1 2 ... n); they generate Sigma_n."""
    gens = []
    for i in range(1, n):
        s = list(range(1, n + 1))
        s[i - 1], s[i] = s[i], s[i - 1]
        gens.append(tuple(s))
    if n > 1:
        gens.append(tuple(list(range(2, n + 1)) + [1]))
    return gens


def act_partition(s, lam):
    return frozenset(frozenset(s[x - 1] for x in B) for B in lam)


def act_chain(s, chain):
    return tuple(act_partition(s, lam) for lam in chain)


def conjugacy_classes(n):
    """dict cycle type -> list of permutations, cycle types sorted with the identity first."""
    classes = {}
    for s in all_perms(n):
        classes.setdefault(cycle_type(s), []).append(s)
    order = sorted(classes, key=lambda ct: (len(ct) * -1, ct))  # (1,..,1) first, (n) last
    return {ct: classes[ct] for ct in order}


def format_cycle_type(ct):
    out, i = [], 0
    while i < len(ct):
        j = i
        while j < len(ct) and ct[j] == ct[i]:
            j += 1
        out.append(str(ct[i]) + ("^%d" % (j - i) if j - i > 1 else ""))
        i = j
    return " ".join(out)


# ----------------------------------------------------------------------------------------------
# exact linear algebra on sparse integer columns (dict row -> int)
# ----------------------------------------------------------------------------------------------

def rank_Q(cols):
    """Rank over Q by exact Gaussian column elimination with fractions.Fraction: each column is reduced
    against the previously found pivot columns, keyed by their lowest nonzero row."""
    pivots = {}
    for col in cols:
        v = {r: Fraction(c) for r, c in col.items()}
        while v:
            r = min(v)
            p = pivots.get(r)
            if p is None:
                a = v[r]
                pivots[r] = {i: c / a for i, c in v.items()}
                break
            a = v[r]
            for i, c in p.items():
                x = v.get(i, 0) - a * c
                if x:
                    v[i] = x
                else:
                    del v[i]
    return len(pivots)


def egcd(a, b):
    """(g, s, t) with s*a + t*b = g = gcd(a, b) >= 0."""
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def z_echelon(cols):
    """Sparse column-echelon form over Z by unimodular column operations (subtracting integer multiples of
    a pivot column, or replacing a pair of columns (p, v) by (s p + t v, -(b/g) p + (a/g) v) with
    s a + t b = g = gcd(a, b), a transformation of determinant 1).  Returns dict lowest_row -> reduced
    column; the reduced columns have distinct lowest rows and their number is the rank of the matrix."""
    pivots = {}
    for col in cols:
        v = dict(col)
        while v:
            r = min(v)
            p = pivots.get(r)
            if p is None:
                pivots[r] = v
                break
            a, b = p[r], v[r]
            if b % a == 0:
                q = b // a
                for i, c in p.items():
                    x = v.get(i, 0) - q * c
                    if x:
                        v[i] = x
                    else:
                        del v[i]
            else:
                g, s, t = egcd(a, b)
                bg, ag = b // g, a // g
                pn, vn = {}, {}
                for i in set(p) | set(v):
                    pi, vi = p.get(i, 0), v.get(i, 0)
                    x, y = s * pi + t * vi, -bg * pi + ag * vi
                    if x:
                        pn[i] = x
                    if y:
                        vn[i] = y
                pivots[r] = pn
                v = vn
    return pivots


def dense(cols, nrows):
    """Dense row-major integer matrix (list of lists) from sparse columns."""
    M = [[0] * len(cols) for _ in range(nrows)]
    for j, col in enumerate(cols):
        for r, c in col.items():
            M[r][j] = c
    return M


def smith_invariants(M):
    """Nonzero diagonal entries d_1 | d_2 | ... of the Smith normal form of an integer matrix M (list of rows),
    computed by unimodular row and column operations.  The rank is len(result); the cokernel is
    Z^{rows - rank} (+) sum Z/d_i."""
    A = [list(r) for r in M]
    m = len(A)
    n = len(A[0]) if m else 0
    diag = []
    t = 0
    while t < min(m, n):
        best = None                                   # nonzero entry of least |value| in A[t:, t:]
        for i in range(t, m):
            row = A[i]
            for j in range(t, n):
                x = row[j]
                if x and (best is None or abs(x) < best[0]):
                    best = (abs(x), i, j)
                    if best[0] == 1:
                        break
            if best and best[0] == 1:
                break
        if best is None:
            break
        _, i, j = best
        A[t], A[i] = A[i], A[t]
        if j != t:
            for row in A:
                row[t], row[j] = row[j], row[t]
        while True:
            if A[t][t] < 0:
                A[t] = [-x for x in A[t]]
            p = A[t][t]
            clean = True
            for i in range(t + 1, m):                 # clear column t
                a = A[i][t]
                if a:
                    q = a // p
                    if q:
                        Ai, At = A[i], A[t]
                        for jj in range(t, n):
                            if At[jj]:
                                Ai[jj] -= q * At[jj]
                    if A[i][t]:
                        clean = False
            for j in range(t + 1, n):                 # clear row t
                a = A[t][j]
                if a:
                    q = a // p
                    if q:
                        for row in A:
                            if row[t]:
                                row[j] -= q * row[t]
                    if A[t][j]:
                        clean = False
            if clean:
                bad = None                            # enforce p | every remaining entry
                for i in range(t + 1, m):
                    if any(x % p for x in A[i][t + 1:]):
                        bad = i
                        break
                if bad is None:
                    break
                At, Ab = A[t], A[bad]
                for jj in range(t, n):
                    At[jj] += Ab[jj]
                continue
            best = None                               # a remainder survived: move the smallest to (t,t)
            for i in range(t, m):
                x = A[i][t]
                if x and (best is None or abs(x) < best[0]):
                    best = (abs(x), i, t)
            for j in range(t, n):
                x = A[t][j]
                if x and (best is None or abs(x) < best[0]):
                    best = (abs(x), t, j)
            _, i, j = best
            A[t], A[i] = A[i], A[t]
            if j != t:
                for row in A:
                    row[t], row[j] = row[j], row[t]
        diag.append(abs(A[t][t]))
        t += 1
    for i in range(len(diag)):                        # normalise to a divisibility chain (already one)
        for j in range(i + 1, len(diag)):
            g = gcd(diag[i], diag[j])
            diag[i], diag[j] = g, diag[i] * diag[j] // g
    return diag


def rank_and_torsion_Z(cols, nrows):
    """(rank over Z, invariant factors > 1, unit_pivots) of the integer matrix with sparse columns cols.
    If every pivot of the Z-echelon form is +-1 the cokernel is torsion-free (see the module docstring);
    otherwise the invariant factors are taken from the dense Smith normal form of the reduced columns."""
    piv = z_echelon(cols)
    r = len(piv)
    if all(abs(p[row]) == 1 for row, p in piv.items()):
        return r, [], True
    d = smith_invariants(dense([piv[row] for row in sorted(piv)], nrows))
    assert len(d) == r, ("Smith rank != echelon rank", len(d), r)
    return r, [x for x in d if x > 1], False


def kernel_basis_Q(cols, ncols):
    """Column reduction over Q of the matrix A with sparse integer columns cols, recording each reduced
    column as a combination u of the original columns.  Returns (free, K, pivots): free = sorted indices of
    the columns that reduced to zero; K[f] = the recorded combination, a sparse Fraction vector with
    A K[f] = 0, K[f][f] = 1 and K[f][f'] = 0 for every other free f' (its remaining support lies in pivot
    columns); {K[f]} is a basis of ker A.  pivots = dict lowest_row -> reduced pivot column."""
    pivots = {}
    K = {}
    for j in range(ncols):
        v = {r: Fraction(c) for r, c in cols[j].items()}
        u = {j: Fraction(1)}
        while v:
            r = min(v)
            p = pivots.get(r)
            if p is None:
                pivots[r] = (v, u)
                break
            pv, pu = p
            a = v[r] / pv[r]
            for i, c in pv.items():
                x = v.get(i, 0) - a * c
                if x:
                    v[i] = x
                else:
                    del v[i]
            for i, c in pu.items():
                x = u.get(i, 0) - a * c
                if x:
                    u[i] = x
                else:
                    del u[i]
        if not v:
            K[j] = u
    free = sorted(K)
    kc = set(free)
    for f in free:                                       # the identity pattern on the free columns
        assert K[f][f] == 1 and all(i == f or i not in kc for i in K[f])
    return free, K, {r: pv for r, (pv, pu) in pivots.items()}


def apply_sparse(cols, v):
    """A v for sparse columns cols and a sparse vector v (dict column -> coefficient)."""
    acc = {}
    for i, x in v.items():
        for r, c in cols[i].items():
            acc[r] = acc.get(r, 0) + x * c
    return {r: c for r, c in acc.items() if c}


# ----------------------------------------------------------------------------------------------
# the chain complex of P(n)
# ----------------------------------------------------------------------------------------------

class PChainComplex:
    def __init__(self, n):
        self.n = n
        self.T = frozenset(range(1, n + 1))
        mn, mx = min_partition(self.T), max_partition(self.T)
        self.basis = {}
        for ch in nondegenerate_chains(self.T):
            assert ch[0] == mn and ch[-1] == mx and not is_degenerate(ch)
            self.basis.setdefault(len(ch) - 1, []).append(ch)
        for k in self.basis:
            self.basis[k].sort(key=ckey)
        self.degrees = sorted(self.basis)
        self.top = max(self.degrees)
        assert self.top == n - 1
        self.index = {k: {ch: i for i, ch in enumerate(b)} for k, b in self.basis.items()}
        self.dim = {k: len(b) for k, b in self.basis.items()}
        self.dropped = {k: {"basepoint": 0, "degenerate": 0} for k in self.degrees}
        self.bd = {}                                          # k -> sparse columns of boundary_k
        for k in self.degrees:
            cols = []
            for ch in self.basis[k]:
                col = {}
                if k > 0:
                    for i in range(k + 1):
                        f = face(ch, i, self.T)
                        if f is None:
                            self.dropped[k]["basepoint"] += 1
                            continue
                        if is_degenerate(f):
                            self.dropped[k]["degenerate"] += 1
                            continue
                        r = self.index[k - 1][f]
                        col[r] = col.get(r, 0) + (-1) ** i
                cols.append({r: c for r, c in col.items() if c})
            self.bd[k] = cols
        self._top = None

    # -- basic checks ------------------------------------------------------------------------
    def check_dd(self):
        """boundary o boundary = 0."""
        for k in self.degrees:
            if k - 1 not in self.bd:
                continue
            for col in self.bd[k]:
                assert not apply_sparse(self.bd[k - 1], col), "dd != 0 in degree %d" % k

    def permutation(self, s, k):
        """img[i] = index of s.(basis[k][i]); the action of s on C_k is e_i -> e_{img[i]}."""
        return [self.index[k][act_chain(s, ch)] for ch in self.basis[k]]

    def check_equivariance(self, perms):
        """boundary_k P(s) = P(s) boundary_k on C_k for all s in perms and all k.  (Since s -> P(s) is a
        homomorphism by construction, equivariance for a generating set implies it for all of Sigma_n.)"""
        for s in perms:
            img = {k: self.permutation(s, k) for k in self.degrees}
            for k in self.degrees:
                if k - 1 not in self.bd:
                    continue
                for j, col in enumerate(self.bd[k]):
                    moved = {img[k - 1][r]: c for r, c in col.items()}
                    assert moved == self.bd[k][img[k][j]], ("boundary not equivariant", s, k, j)

    def fixed_counts(self, s):
        """Number of basis chains of C_k fixed by s, per degree k (= trace of s on C_k)."""
        return {k: sum(1 for ch in self.basis[k] if act_chain(s, ch) == ch) for k in self.degrees}

    # -- homology ----------------------------------------------------------------------------
    def homology(self, dense_snf=False):
        """dict k -> dict(dim, rank_Q, rank_Z, betti, torsion, unit_pivots): H~_k = Z^betti (+) sum Z/d for
        d in torsion.  rank_Q by Fraction elimination; rank_Z and torsion from the Z-echelon form of the
        boundaries (torsion of H~_k = torsion of C_k / im boundary_{k+1}).  With dense_snf=True the dense
        Smith normal form of every boundary matrix is computed as well and must agree."""
        rq = {k: rank_Q(self.bd[k]) for k in self.degrees}
        rz, tors, units = {}, {}, {}
        for k in self.degrees:
            if k - 1 in self.bd:
                rz[k], tors[k], units[k] = rank_and_torsion_Z(self.bd[k], self.dim[k - 1])
                assert rz[k] == rq[k], ("rank over Z != rank over Q", k, rz[k], rq[k])
                if dense_snf:
                    d = smith_invariants(dense(self.bd[k], self.dim[k - 1]))
                    assert len(d) == rq[k] and [x for x in d if x > 1] == tors[k], ("dense SNF disagrees", k, d)
            else:
                rz[k], tors[k], units[k] = 0, [], True
        out = {}
        for k in self.degrees:
            out[k] = dict(dim=self.dim[k], rank_Q=rq[k], rank_Z=rz[k],
                          betti=self.dim[k] - rq[k] - rq.get(k + 1, 0),
                          torsion=tors.get(k + 1, []), unit_pivots=units.get(k + 1, True))
        return out

    # -- top homology and the Sigma_n action ---------------------------------------------------
    def top_cycle_basis(self):
        """Basis of H~_{n-1}(P(n); Q) = ker boundary_{n-1} as explicit cycles.  Returns (free, K, integral):
        free = the free column indices f_1 < ... < f_m, K = list of the cycles K_1..K_m (sparse vectors on
        the basis chains of C_{n-1}, K_i[f_i] = 1, K_i[f_j] = 0 for j != i), integral = True iff all K_i
        are integer vectors (then they are a Z-basis of the integer cycles and are stored with int entries).
        Each K_i is verified to be a cycle."""
        k = self.top
        A = self.bd[k]
        free, Kd, _ = kernel_basis_Q(A, self.dim[k])
        K = [Kd[f] for f in free]
        integral = all(c.denominator == 1 for v in K for c in v.values())
        if integral:
            K = [{i: int(c) for i, c in v.items()} for v in K]
        for v in K:
            assert not apply_sparse(A, v), "kernel vector is not a cycle"
        assert len(K) == self.dim[k] - rank_Q(A), "kernel dimension != nullity"
        self._top = (free, K, integral)
        return self._top

    def coordinates(self, w, verify=True):
        """Coordinates of the sparse cycle w in the basis K: x_i = w[f_i].  With verify=True the vector is
        reconstructed from its coordinates and compared with w (this is what certifies w in span K)."""
        free, K, _ = self._top
        x = [w.get(f, 0) for f in free]
        if verify:
            recon = {}
            for xi, Ki in zip(x, K):
                if xi:
                    for i, c in Ki.items():
                        recon[i] = recon.get(i, 0) + xi * c
            assert {i: c for i, c in recon.items() if c} == {i: c for i, c in w.items() if c}, \
                "coordinates do not reconstruct the vector"
        return x

    def top_action_matrix(self, s, verify=True):
        """Matrix M(s) of s on H~_{n-1} in the basis K (row-major; column j = coordinates of s.K_j, i.e.
        M[i][j] = (s.K_j)[f_i]); integer entries when the basis is integral, Fractions otherwise."""
        free, K, _ = self._top
        img = self.permutation(s, self.top)
        m = len(K)
        M = [[0] * m for _ in range(m)]
        for j, v in enumerate(K):
            w = {img[i]: c for i, c in v.items()}
            for i, xi in enumerate(self.coordinates(w, verify=verify)):
                M[i][j] = xi
        return M

    def top_trace(self, s):
        """Trace of s on H~_{n-1}: sum_i (s.K_i)[f_i] = sum_i K_i[img^{-1}(f_i)], without forming M(s)."""
        free, K, _ = self._top
        img = self.permutation(s, self.top)
        inv = [0] * len(img)
        for i, x in enumerate(img):
            inv[x] = i
        return sum(Ki.get(inv[f], 0) for f, Ki in zip(free, K))


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(A[i], Bt[j])) for j in range(p)] for i in range(n)]


def trace(M):
    return sum(M[i][i] for i in range(len(M)))


# ----------------------------------------------------------------------------------------------
# characters of Sigma_n
# ----------------------------------------------------------------------------------------------

def moebius(m):
    mu, p = 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            mu = -mu
        p += 1
    if m > 1:
        mu = -mu
    return mu


def induced_lie_character(n):
    """Character of Lie(n) = Ind_{C_n}^{Sigma_n}(omega) (Klyachko), C_n = <c>, c = (1 2 ... n),
    omega(c^j) = exp(2 pi i j/n), by brute force from the induced-character formula
        chi(g) = (1/|C_n|) sum_{h in Sigma_n : h g h^{-1} in C_n} omega(h g h^{-1}).
    The complex sum is evaluated in floating point and must be an integer to 1e-9; it is cross-checked
    exactly: the multiplicity of c^j depends only on gcd(j, n), and the primitive m-th roots of unity sum to
    moebius(m).  Returns dict g -> integer value, for every g in Sigma_n."""
    import cmath
    perms = all_perms(n)
    c = tuple(list(range(2, n + 1)) + [1])
    Cn, p = {}, identity(n)
    for j in range(n):
        Cn[p] = j
        p = compose(c, p)
    assert len(Cn) == n and p == identity(n)
    roots = [cmath.exp(2j * cmath.pi * jj / n) for jj in range(n)]
    conj = [(h, inverse(h)) for h in perms]
    chi = {}
    for g in perms:
        counts = [0] * n
        for h, hinv in conj:
            j = Cn.get(compose(compose(h, g), hinv))
            if j is not None:
                counts[j] += 1
        z = sum(m * roots[jj] for jj, m in enumerate(counts)) / n
        assert abs(z.imag) < 1e-9 and abs(z.real - round(z.real)) < 1e-9, (g, z)
        val = int(round(z.real))
        by_d = {}
        for jj, m in enumerate(counts):
            d = gcd(jj, n)
            assert by_d.setdefault(d, m) == m
        exact = sum(m * moebius(n // d) for d, m in by_d.items())
        assert exact % n == 0 and exact // n == val, (g, exact, val)
        chi[g] = val
    return chi


def klyachko_closed_form(ct):
    """Reference only (not used in any assertion): mu(m) m^(d-1) (d-1)! for cycle type m^d, else 0."""
    from math import factorial
    n, m, d = sum(ct), ct[0], len(ct)
    if any(x != m for x in ct):
        return 0
    return moebius(m) * m ** (d - 1) * factorial(d - 1)


# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    import time
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(1, nmax + 1):
        t0 = time.time()
        cc = PChainComplex(n)
        cc.check_dd()
        hom = cc.homology(dense_snf=(n <= 5))
        free, K, integral = cc.top_cycle_basis()
        classes = conjugacy_classes(n)
        chi = {ct: cc.top_trace(cl[0]) for ct, cl in classes.items()}
        print("n = %d  (%.2fs)" % (n, time.time() - t0))
        for k in cc.degrees:
            h = hom[k]
            print("  degree %d: dim C_k = %4d  rank d_k (Q) = %4d  rank d_k (Z) = %4d  betti = %3d  torsion = %s  "
                  "faces dropped: basepoint %d, degenerate %d"
                  % (k, h["dim"], h["rank_Q"], h["rank_Z"], h["betti"], h["torsion"] or "none",
                     cc.dropped[k]["basepoint"], cc.dropped[k]["degenerate"]))
        print("  top cycle basis: %d cycles, integral: %s, max support %d; character on H~_%d: %s"
              % (len(K), integral, max(len(v) for v in K), n - 1,
                 ", ".join("%s: %d" % (format_cycle_type(ct), chi[ct]) for ct in classes)))
