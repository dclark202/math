"""Degree counting for the two cocomposition maps compared in Theorem 1.3 (check C11).

For a partition alpha = {T_1, ..., T_n} of {1..k} (blocks ordered by least element, |T_i| = k_i) both operad
structures on Map(|P(-)|, S) are dual to a based map

    c_alpha : |P(k)|  -->  |P(n)| ^ |P(k_1)| ^ ... ^ |P(k_n)|          (Theta used to identify B(t) = |P(t)|),

namely
  * `lev`: the levelled cut c^lev_alpha = kappa_{F^{1/2}} of Example A.50(2) = the half-cut cell psi_0 of
    Construction A.68 acting through the decomposition Psi of Lemma 4.4 (Proposition A.71, checked pointwise in
    ../checks_transport.py): on the cell (gamma', u) put P = max{i : u_i <= 1/2}; the value is the basepoint unless
    lambda'_P = alpha, and otherwise it is (front chain collapsed at alpha, 2u_1, ..., 2u_P) smashed with the
    restrictions of the back segment to the blocks, with heights 2u_{P+1}-1, ..., 2u_{k-1}-1, degenerate simplices
    normalized ((s_j beta, w) ~ (beta, sigma^j w): the repeated entry and the height between the repeats are dropped);
  * `tr`: Ching's cut c^tr_alpha (Definition A.31) transported through Theta: Theta(gamma', u), degrafting
    (Lemma A.28), bottom heights unchanged, top piece on T_i rescaled by [pi_i, 1] -> [0, 1] with
    pi_i = h(p(T_i)) (pi_i = 0 for T_i = {1..k}), then Theta^{-1} on each piece and the order-preserving
    identifications T_i = {1..k_i}.

Both maps are (piecewise) smooth bijections from parts of the open top cells of the source onto the open top
cells of the target, so the coefficient of a target top cell e in c_*(z), z = sum_{gamma'} z_{gamma'} gamma' a top
cycle, is the local degree sum over the preimages of a generic point y of e:
    c_*(z)[e] = sum_{(gamma', u) : c(gamma', u) = y}  z_{gamma'} * sign det D(c|_{gamma'})(u),
the sign taken with respect to the orientations fixed below.  The preimages are found exactly: on each source top
cell the combinatorial output (the target chains and which source level feeds which target coordinate) is
constant, and the coordinates are inverted by hand (`CellMap.invert`); every preimage found is re-verified by
evaluating the forward map.  All arithmetic is in fractions.Fraction.

ORIENTATION CONVENTIONS
  (O1) A nondegenerate p-simplex of P(t), i.e. a strict chain lambda_0 < ... < lambda_p, is oriented by the order of
       the chain (as in chains.py, where the boundary is sum (-1)^i d_i).  In height coordinates u_j = t_0 + ... +
       t_{j-1} the vertex lambda_i is the point (0^i, 1^{p-i}), so the chain orientation is (-1)^p times the
       orientation of R^p given by (u_1, ..., u_p).  Source and target cells both have dimension k-1, so this
       factor cancels and the local degree is the sign of the Jacobian of the map in height coordinates.
  (O2) A product cell gamma x beta_1 x ... x beta_n of the smash product carries the product orientation in the
       order (n-factor, then the k_i-factors in the order of the blocks), i.e. the coordinates are listed as
       (v_1..v_{n-1}, w^1_1..w^1_{k_1-1}, ..., w^n_1..w^n_{k_n-1}).  A factor |P(1)| = S^0 has no coordinates.
  (O3) H~_{k-1}(target) = H~_{n-1}(P(n)) (x) H~_{k_1-1}(P(k_1)) (x) ... (Kunneth; every factor is a wedge of
       spheres of a single dimension): the basis element K_a (x) K^1_{b_1} (x) ... (x) K^n_{b_n} is the cellular
       cycle sum K_a[gamma] prod_i K^i_{b_i}[beta_i] * (gamma x beta_1 x ... x beta_n) with the orientation (O2).
       With the cycle bases of chains.py (identity pattern on the free columns f) the coordinates of a top cycle
       of the target are read off at the product cells (f_a, f^1_{b_1}, ..., f^n_{b_n}), and the reconstruction
       is verified, which certifies that the computed image is a cycle.
  (O4) The levelling homeomorphism Theta and its inverse are used with the height coordinates of Lemma A.23:
       a vertex B has height u_{b(B)+1}, and the chain of a weighted tree lists the alive partitions at the
       distinct heights in increasing order.  This is a coordinate identification, not a reorientation:
       transported through Theta, c^tr sends (gamma', u) to a target point whose height coordinates are
       coordinate selections of u (bottom) and increasing affine images of coordinates of u (tops), so no sign
       is introduced by Theta itself.
  (O5) Sigma_k acts on cells by relabelling, with coefficient +1 (relabelling preserves the order of a chain and
       the heights).  On the target the induced map sends the factor at T_i to the factor at sigma(T_i); reordering
       the factors reorders the coordinate blocks, which contributes the Koszul sign
       prod_{i<i', sbar(i)>sbar(i')} (-1)^{(k_i-1)(k_{i'}-1)}.
  With these conventions the Jacobian of either map on a source cell, in the coordinates (u_1..u_{k-1}) ->
  (v; w^1; ...; w^n), is block triangular with positive diagonal blocks (coordinate selection for the bottom,
  positive multiples of the identity for the tops), so its sign is the sign of the permutation of {1..k-1}
  listing the source levels in target-coordinate order ("shuffle sign").  This is what `CellMap.sign` records;
  `jacobian_sign_fd` recomputes the sign from the forward map by exact central differences, independently of
  that analysis.
"""
import os
import sys
from fractions import Fraction
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, os.pardir))
from partitions import decompose, max_partition  # noqa: E402
from trees import tree_of_chain, levels, parent, root, alive  # noqa: E402
from chains import PChainComplex, act_partition, act_chain, is_degenerate, mat_mul  # noqa: E402

ZERO, ONE, HALF = Fraction(0), Fraction(1), Fraction(1, 2)
POINT = (frozenset([frozenset([1])]),)          # the non-basepoint 0-simplex of P(1)


# ----------------------------------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------------------------------

def perm_list_sign(seq):
    """Sign of the permutation given as a list of distinct integers (sign of the sorting permutation)."""
    inv = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                inv += 1
    return -1 if inv % 2 else 1


def det_sign_exact(M):
    """Sign of the determinant of a square matrix of Fractions (0 if singular); 1 for the empty matrix."""
    A = [[Fraction(x) for x in row] for row in M]
    n = len(A)
    sgn = 1
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            sgn = -sgn
        if A[c][c] < 0:
            sgn = -sgn
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return sgn


def blocks_of(alpha):
    """The blocks of a partition in the order used everywhere (least element first; = partitions.decompose)."""
    return sorted(alpha, key=min)


def profile_of(alpha):
    return tuple(len(b) for b in blocks_of(alpha))


def order_iso(Tb):
    """The order-preserving bijection T_b -> {1..|T_b|} as a dict."""
    return {x: i for i, x in enumerate(sorted(Tb), 1)}


def relabel_chain(chain, m):
    return tuple(frozenset(frozenset(m[x] for x in B) for B in lam) for lam in chain)


def normalize(chain, heights):
    """(s_j beta, w) ~ (beta, sigma^j w): drop every repeated entry of a degenerate chain together with the height
    of the step into it (height u_r sits between lambda_{r-1} and lambda_r)."""
    assert len(heights) == len(chain) - 1
    ch, hs = [chain[0]], []
    for r in range(1, len(chain)):
        if chain[r] != chain[r - 1]:
            ch.append(chain[r])
            hs.append(heights[r - 1])
    return tuple(ch), tuple(hs)


def split_blocks(chain):
    """B_j (1 <= j <= dim): the unique block of lambda_{j-1} absent from lambda_j; the height u_j is the height of
    B_j under Theta (b(B_j) = j-1).  Only for maximal chains, where exactly one block splits per step."""
    out = []
    for j in range(1, len(chain)):
        gone = [B for B in chain[j - 1] if B not in chain[j]]
        assert len(gone) == 1, "not a maximal chain"
        out.append(gone[0])
    return out


def all_partitions(k):
    from partitions import partitions_of
    return sorted(partitions_of(frozenset(range(1, k + 1))), key=lambda a: (len(a), sorted(sorted(b) for b in a)))


# ----------------------------------------------------------------------------------------------
# Theta and its inverse, exactly
# ----------------------------------------------------------------------------------------------

def theta_exact(chain, u):
    """Construction A.22: h(B) = u_{b(B)+1}; u is a tuple of Fractions of length dim(chain)."""
    F = tree_of_chain(chain)
    a, b = levels(chain)
    return F, {B: u[b[B]] for B in F if len(B) > 1}


def theta_inverse_exact(F, h):
    """Proof of Theorem A.24 (III): u = the distinct heights in increasing order, lambda_j = alive partition at
    u_{j+1}, lambda_last = max.  For a tree without internal vertices returns the 0-simplex."""
    T = root(F)
    if not h:
        return (max_partition(T),), ()
    ws = sorted(set(h.values()))
    return tuple(alive(F, h, w) for w in ws) + (max_partition(T),), tuple(ws)


# ----------------------------------------------------------------------------------------------
# the two forward maps on a cell (gamma', u) of |P(k)|; value None = basepoint, otherwise
# ((gamma, v), ((beta_1, w^1), ..., (beta_n, w^n))) with gamma a chain on {1..n}, beta_i on {1..k_i}
# ----------------------------------------------------------------------------------------------

def lev_forward(chain, u, alpha):
    """The levelled cut c^lev_alpha in chain coordinates (Construction A.68 through Psi, Proposition A.71)."""
    m = len(chain) - 1
    uu = [ZERO] + list(u) + [ONE]
    P = max(i for i in range(m + 1) if uu[i] <= HALF)
    if chain[P] != alpha:
        return None
    gamma, betas, blocks = decompose(chain, P, alpha)
    front = relabel_chain(gamma, {i: i + 1 for i in range(len(blocks))})
    assert not is_degenerate(front)
    v = tuple(2 * x for x in u[:P])
    back = tuple(2 * x - 1 for x in u[P:])
    pieces = []
    for beta, Tb in zip(betas, blocks):
        bn, wn = normalize(beta, back)
        pieces.append((relabel_chain(bn, order_iso(Tb)), wn))
    return (front, v), tuple(pieces)


def tr_forward(chain, u, alpha):
    """Ching's cut c^tr_alpha (Definition A.31), transported through Theta."""
    F, h = theta_exact(chain, u)
    T = root(F)
    blocks = blocks_of(alpha)
    if any(len(Tb) > 1 and Tb not in F for Tb in blocks):
        return None
    inside = {B: next((i for i, Tb in enumerate(blocks) if B <= Tb), None) for B in F}

    def q(B):
        return frozenset(i + 1 for i, Tb in enumerate(blocks) if Tb <= B)
    G = frozenset(q(B) for B in F if not any(B < Tb for Tb in blocks))
    h_bot = {q(B): h[B] for B in F if len(B) > 1 and inside[B] is None}
    bottom = theta_inverse_exact(G, h_bot)
    pieces = []
    for Tb in blocks:
        if len(Tb) == 1:
            pieces.append((POINT, ()))
            continue
        pi = ZERO if Tb == T else h[parent(F, Tb)]
        U = frozenset(B for B in F if B <= Tb)
        hi = {B: (h[B] - pi) / (1 - pi) for B in U if len(B) > 1}
        ch, w = theta_inverse_exact(U, hi)
        pieces.append((relabel_chain(ch, order_iso(Tb)), w))
    return bottom, tuple(pieces)


FORWARD = {"lev": lev_forward, "tr": tr_forward}


def target_key(y):
    """The target cell of a non-basepoint value: (gamma, beta_1, ..., beta_n)."""
    return (y[0][0],) + tuple(b for b, w in y[1])


def target_heights(y):
    """The target coordinates (v, (w^1, ..., w^n))."""
    return y[0][1], tuple(w for b, w in y[1])


def flat_heights(y):
    v, ws = target_heights(y)
    return list(v) + [x for w in ws for x in w]


# ----------------------------------------------------------------------------------------------
# the restriction of a map to one source top cell: combinatorial output + exact inversion
# ----------------------------------------------------------------------------------------------

class CellMap:
    """c restricted to (the top-cell-contributing region of) the source top cell `chain`.
    target = the target cell hit; levels = the source levels (1-based) in target-coordinate order; sign = the
    shuffle sign (= sign of the Jacobian, see the module docstring); invert(v, ws) = the unique candidate
    preimage u (a list of Fractions), to be tested for 0 < u_1 < ... < u_{k-1} < 1."""

    def __init__(self, structure, chain, target, levels_, invert, ref):
        self.structure, self.chain, self.target, self.levels, self.invert, self.ref = structure, chain, target, levels_, invert, ref
        assert sorted(levels_) == list(range(1, len(chain)))
        self.sign = perm_list_sign(levels_)


def lev_cellmaps(chain, alpha):
    m = len(chain) - 1
    blocks = blocks_of(alpha)
    n = len(blocks)
    P = n - 1                                          # lambda'_P = alpha forces P + 1 = n blocks
    if chain[P] != alpha:
        return []
    Bj = split_blocks(chain)
    J = [[] for _ in blocks]
    for j in range(P + 1, m + 1):
        i = next(i for i, Tb in enumerate(blocks) if Bj[j - 1] <= Tb)
        J[i].append(j)
    assert all(len(J[i]) == len(Tb) - 1 for i, Tb in enumerate(blocks))
    lv = list(range(1, P + 1)) + [j for Ji in J for j in Ji]
    ref = tuple([Fraction(j, 2 * (P + 1)) for j in range(1, P + 1)] +
                [HALF + Fraction(j - P, 2 * (m - P + 1)) for j in range(P + 1, m + 1)])
    y = lev_forward(chain, ref, alpha)
    assert y is not None

    def invert(v, ws):
        u = [None] * m
        for r in range(P):
            u[r] = v[r] / 2
        for i, Ji in enumerate(J):
            for s, j in enumerate(Ji):
                u[j - 1] = (ws[i][s] + 1) / 2
        return u
    return [CellMap("lev", chain, target_key(y), lv, invert, ref)]


def tr_cellmaps(chain, alpha):
    m = len(chain) - 1
    F = tree_of_chain(chain)
    T = root(F)
    blocks = blocks_of(alpha)
    if any(len(Tb) > 1 and Tb not in F for Tb in blocks):
        return []
    Bj = split_blocks(chain)
    level_of = {B: j for j, B in enumerate(Bj, 1)}
    JG = [j for j, B in enumerate(Bj, 1) if not any(B <= Tb for Tb in blocks)]
    J = [[j for j, B in enumerate(Bj, 1) if B <= Tb] for Tb in blocks]
    jp = [None if len(Tb) == 1 else (0 if Tb == T else level_of[parent(F, Tb)]) for Tb in blocks]
    assert all(j is None or j == 0 or j in JG for j in jp)
    lv = JG + [j for Ji in J for j in Ji]
    ref = tuple(Fraction(j, m + 1) for j in range(1, m + 1))
    y = tr_forward(chain, ref, alpha)
    assert y is not None

    def invert(v, ws):
        u = [None] * m
        for r, j in enumerate(JG):
            u[j - 1] = v[r]
        for i, Ji in enumerate(J):
            if not Ji:
                continue
            pi = ZERO if jp[i] == 0 else u[jp[i] - 1]
            for s, j in enumerate(Ji):
                u[j - 1] = pi + (1 - pi) * ws[i][s]
        return u
    return [CellMap("tr", chain, target_key(y), lv, invert, ref)]


CELLMAPS = {"lev": lev_cellmaps, "tr": tr_cellmaps}


def jacobian_sign_fd(cm, alpha, h=Fraction(1, 10 ** 7)):
    """Sign of det of the Jacobian of the forward map at the reference point of the cell, in the coordinates of
    (O1)-(O2), by exact central differences (exact for the affine levelled map; error O(h^2) for the rational tree
    map, negligible against |det| >= 1)."""
    fwd = FORWARD[cm.structure]
    m = len(cm.chain) - 1
    u0 = list(cm.ref)
    cols = []
    for j in range(m):
        up, um = list(u0), list(u0)
        up[j] += h
        um[j] -= h
        yp, ym = fwd(cm.chain, tuple(up), alpha), fwd(cm.chain, tuple(um), alpha)
        assert yp is not None and ym is not None and target_key(yp) == cm.target and target_key(ym) == cm.target
        fp, fm = flat_heights(yp), flat_heights(ym)
        cols.append([(a - b) / (2 * h) for a, b in zip(fp, fm)])
    M = [[cols[j][r] for j in range(m)] for r in range(m)]
    return det_sign_exact(M)


# ----------------------------------------------------------------------------------------------
# generic points, preimages, degrees
# ----------------------------------------------------------------------------------------------

class NonGeneric(Exception):
    pass


def random_heights(p, rng, D=10 ** 6):
    while True:
        vals = sorted(Fraction(rng.randint(1, D - 1), D) for _ in range(p))
        if len(set(vals)) == p:
            return tuple(vals)


def random_target_point(e, rng):
    return random_heights(len(e[0]) - 1, rng), tuple(random_heights(len(b) - 1, rng) for b in e[1:])


def preimages(cellmaps_of_e, e, v, ws, alpha):
    """All preimages of the target point (e; v, ws): the CellMaps whose inversion lands in the open source cell.
    Each candidate is re-verified by the forward map.  Raises NonGeneric at a tie."""
    out = []
    for cm in cellmaps_of_e:
        u = cm.invert(v, ws)
        assert all(x is not None for x in u)
        ok, prev = True, ZERO
        for x in u:
            if x == prev or x >= ONE:
                raise NonGeneric()
            if x < prev:
                ok = False
                break
            prev = x
        if not ok:
            continue
        y = FORWARD[cm.structure](cm.chain, tuple(u), alpha)
        assert y is not None and target_key(y) == e and target_heights(y) == (v, ws), "inversion does not reproduce y"
        out.append(cm)
    return out


class Structure:
    """All CellMaps of one structure for one partition alpha, grouped by target cell."""

    def __init__(self, structure, k, alpha, top_chains):
        self.structure, self.k, self.alpha = structure, k, alpha
        self.blocks = blocks_of(alpha)
        self.profile = tuple(len(b) for b in self.blocks)
        self.by_target = {}
        self.cellmaps = []
        for ch in top_chains:
            for cm in CELLMAPS[structure](ch, alpha):
                self.by_target.setdefault(cm.target, []).append(cm)
                self.cellmaps.append(cm)

    def degree_vector(self, e, cycles, index, rng, npoints=3, max_resample=50):
        """[sum_{preimages of y} z[chain] * sign for z in cycles] for npoints random generic y in the cell e;
        asserts the vectors agree and returns one of them together with the preimage counts."""
        vecs, counts = [], []
        cms = self.by_target.get(e, [])
        t = 0
        tries = 0
        while t < npoints:
            v, ws = random_target_point(e, rng)
            try:
                pre = preimages(cms, e, v, ws, self.alpha)
            except NonGeneric:
                tries += 1
                assert tries < max_resample
                continue
            vecs.append([sum(z.get(index[cm.chain], 0) * cm.sign for cm in pre) for z in cycles])
            counts.append(len(pre))
            t += 1
        assert all(vec == vecs[0] for vec in vecs), ("degree sum depends on the generic point", self.structure, e, vecs)
        return vecs[0], counts


# ----------------------------------------------------------------------------------------------
# homology: Kunneth bases and induced matrices
# ----------------------------------------------------------------------------------------------

class Complexes:
    """Cache of PChainComplex(t) with top cycle bases."""

    def __init__(self):
        self.cx = {}
        self.product_cycles = {}

    def __call__(self, t):
        if t not in self.cx:
            cc = PChainComplex(t)
            cc.top_cycle_basis()
            self.cx[t] = cc
        return self.cx[t]


def target_cells(cxs, profile):
    n = len(profile)
    factors = [cxs(n).basis[n - 1]] + [cxs(ki).basis[ki - 1] for ki in profile]
    return [tuple(cell) for cell in product(*factors)]


def target_basis(cxs, profile):
    """Index tuples (a, b_1, ..., b_n) of the Kunneth basis, lexicographic; a indexes the basis of H~(P(n))."""
    n = len(profile)
    ranges = [range(len(cxs(n)._top[1]))] + [range(len(cxs(ki)._top[1])) for ki in profile]
    return [tuple(ix) for ix in product(*ranges)]


def product_cycle_coefficient(cxs, profile, ix, e):
    """Coefficient of the product cell e in the Kunneth basis element ix (O3)."""
    n = len(profile)
    c = cxs(n)._top[1][ix[0]].get(cxs(n).index[n - 1][e[0]], 0)
    for i, ki in enumerate(profile):
        if not c:
            return 0
        c *= cxs(ki)._top[1][ix[i + 1]].get(cxs(ki).index[ki - 1][e[i + 1]], 0)
    return c


def product_cycle(cxs, profile, ix):
    """The Kunneth basis cycle ix (O3) as a sparse vector on product cells; its support is the product of the
    supports of the factor cycles, so this is cheap even when the target has thousands of cells."""
    key = (profile, ix)
    cache = cxs.product_cycles
    if key not in cache:
        n = len(profile)
        factors = [(cxs(n), ix[0])] + [(cxs(ki), ix[i + 1]) for i, ki in enumerate(profile)]
        vec = {(): 1}
        for cc, b in factors:
            K, basis = cc._top[1][b], cc.basis[cc.top]
            vec = {cell + (basis[j],): c * x for cell, c in vec.items() for j, x in K.items()}
        cache[key] = vec
    return cache[key]


def induced_matrix(structure, k, alpha, cxs, rng, npoints=3, stats=None):
    """The matrix of (c_alpha)_* : H~_{k-1}(P(k)) -> H~_{n-1}(P(n)) (x) H~_{k_1-1}(P(k_1)) (x) ... in the cycle
    bases of chains.py (columns) and the Kunneth basis (rows, order of `target_basis`).  Returns (M, info)."""
    cc = cxs(k)
    free, K, integral = cc._top
    assert integral
    index = cc.index[k - 1]
    S = Structure(structure, k, alpha, cc.basis[k - 1])
    profile = S.profile
    n = len(profile)
    cells = target_cells(cxs, profile)
    image = [{} for _ in K]                      # image[a] : target cell -> coefficient
    npre = {}
    for e in cells:
        vec, counts = S.degree_vector(e, K, index, rng, npoints)
        npre[e] = counts
        for a, c in enumerate(vec):
            if c:
                image[a][e] = c
    basis = target_basis(cxs, profile)
    frees = [cxs(n)._top[0]] + [cxs(ki)._top[0] for ki in profile]
    free_cells = {ix: tuple(cxs(n).basis[n - 1][frees[0][ix[0]]] for _ in [0]) +
                  tuple(cxs(ki).basis[ki - 1][frees[i + 1][ix[i + 1]]] for i, ki in enumerate(profile)) for ix in basis}
    M = [[0] * len(K) for _ in basis]
    for a in range(len(K)):
        coords = {ix: image[a].get(free_cells[ix], 0) for ix in basis}
        for r, ix in enumerate(basis):
            M[r][a] = coords[ix]
        # certify: the image is the corresponding combination of product cycles (hence a cycle).  This is a
        # strong consistency check on the degree sums: they must satisfy every cycle relation of the target.
        recon = {}
        for ix, x in coords.items():
            if x:
                for e, c in product_cycle(cxs, profile, ix).items():
                    recon[e] = recon.get(e, 0) + x * c
        recon = {e: c for e, c in recon.items() if c}
        assert recon == image[a], ("image of cycle %d is not in the span of the product cycles" % a, structure, alpha)
    info = dict(profile=profile, ncells=len(cells), ncellmaps=len(S.cellmaps), npre=npre, structure=S)
    if stats is not None:
        stats["cells"] = stats.get("cells", 0) + len(cells)
        stats["cellmaps"] = stats.get("cellmaps", 0) + len(S.cellmaps)
    return M, info


# ----------------------------------------------------------------------------------------------
# equivariance data (O5)
# ----------------------------------------------------------------------------------------------

def induced_block_data(alpha, sigma):
    """For sigma in Sigma_k (tuple, sigma[x-1] = sigma(x)): the block permutation sbar (0-based list, block i of
    alpha -> block sbar[i] of sigma.alpha), the permutations sigma_i of {1..k_i} induced through the
    order-preserving identifications, and the Koszul sign of the reordering of the smash factors."""
    blocks = blocks_of(alpha)
    salpha = act_partition(sigma, alpha)
    sblocks = blocks_of(salpha)
    sbar = [sblocks.index(frozenset(sigma[x - 1] for x in Tb)) for Tb in blocks]
    sigmas = []
    for i, Tb in enumerate(blocks):
        m_in, m_out = order_iso(Tb), order_iso(sblocks[sbar[i]])
        inv_in = {p: x for x, p in m_in.items()}
        sigmas.append(tuple(m_out[sigma[inv_in[p] - 1]] for p in range(1, len(Tb) + 1)))
    dims = [len(Tb) - 1 for Tb in blocks]
    koszul = 1
    for i in range(len(blocks)):
        for i2 in range(i + 1, len(blocks)):
            if sbar[i] > sbar[i2] and dims[i] % 2 and dims[i2] % 2:
                koszul = -koszul
    return salpha, sbar, sigmas, koszul


def act_target_point(y, alpha, sigma):
    """The induced map target(alpha) -> target(sigma.alpha) on points (no signs, (O5)): relabel the bottom chain
    by sbar, send the factor at T_i to the slot sbar(i) relabelled by sigma_i; all heights unchanged."""
    if y is None:
        return None
    salpha, sbar, sigmas, koszul = induced_block_data(alpha, sigma)
    (gamma, v), pieces = y
    sb = tuple(s + 1 for s in sbar)
    new_pieces = [None] * len(sbar)
    for i, (beta, w) in enumerate(pieces):
        new_pieces[sbar[i]] = (act_chain(sigmas[i], beta), w)
    return (act_chain(sb, gamma), v), tuple(new_pieces)


def target_action_matrix(alpha, sigma, cxs):
    """Matrix of the induced map H~(target of alpha) -> H~(target of sigma.alpha) in the Kunneth bases."""
    salpha, sbar, sigmas, koszul = induced_block_data(alpha, sigma)
    profile, sprofile = profile_of(alpha), profile_of(salpha)
    n = len(profile)
    Mn = cxs(n).top_action_matrix(tuple(s + 1 for s in sbar))
    Mi = [cxs(ki).top_action_matrix(sigmas[i]) for i, ki in enumerate(profile)]
    rows, cols = target_basis(cxs, sprofile), target_basis(cxs, profile)
    T = [[0] * len(cols) for _ in rows]
    for r, ix2 in enumerate(rows):
        for c, ix in enumerate(cols):
            val = koszul * Mn[ix2[0]][ix[0]]
            for i in range(n):
                if not val:
                    break
                val *= Mi[i][ix2[sbar[i] + 1]][ix[i + 1]]
            T[r][c] = val
    return T, salpha


def matrices_equal(A, B):
    return len(A) == len(B) and all(ra == rb for ra, rb in zip(A, B))


def fmt_matrix(M, indent="    "):
    if not M:
        return indent + "(empty)"
    w = max(len(str(x)) for row in M for x in row)
    return "\n".join(indent + "[" + " ".join(str(x).rjust(w) for x in row) + "]" for row in M)


def fmt_partition(alpha):
    return "|".join("".join(map(str, sorted(b))) for b in blocks_of(alpha))


__all__ = [n for n in dir() if not n.startswith("_")]
