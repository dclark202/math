"""Flags, positions, pieces (Sec. A.46), framed cut systems and their cut maps (Def A.48),
tree frames / labelled frames (Example A.50), Ching's cut c^tr (Def A.31).
A weighted tree is (F, h): F a laminar family (frozenset of frozensets), h: internal vertex -> float.
The basepoint is returned as None."""
from trees import internal, root, parent, children, leaf_parent, alive
from partitions import refines


def positions(flag):
    """All positions (j, B), 0 <= j <= l-1, B a block of flag[j]."""
    return [(j, B) for j in range(len(flag) - 1) for B in flag[j]]


def child_set(flag, pos):
    j, B = pos
    return frozenset(C for C in flag[j + 1] if C <= B)


def nonunary(flag):
    return [p for p in positions(flag) if len(child_set(flag, p)) >= 2]


def parent_position(flag, pos):
    j, B = pos
    if j == 0:
        return None
    return (j - 1, next(C for C in flag[j - 1] if B <= C))


def flag_present(F, flag):
    return all(B in F for lam in flag for B in lam if len(B) > 1)


def piece(F, h, pos, flag):
    """Raw piece (V_rho, h|raw): tree on the children set, ambient heights."""
    j, B = pos
    Cs = child_set(flag, pos)

    def q(D):
        return frozenset(C for C in Cs if C <= D)
    members = [D for D in F if D <= B and any(C <= D for C in Cs)]
    V = frozenset(q(D) for D in members)
    hp = {q(D): h[D] for D in members if len(q(D)) > 1}
    return V, hp


def rescale_piece(V, hp, S, E):
    """(A.24) on a raw piece; None at a threshold (S >= E, root height <= 0, a leaf-parent height >= 1)."""
    if not (E > S):
        return None
    hn = {D: (hp[D] - S) / (E - S) for D in hp}
    if not hn:
        return (V, hn)
    r = root(V)
    if not (hn[r] > 0):
        return None
    for x in r:
        if not (hn[leaf_parent(V, x)] < 1):
            return None
    return V, hn


def kappa(F, h, flag, frames):
    """Cut map of a framed cut system. frames: (F, h, pos) -> (S, E). Returns {pos: (V, h)} or None."""
    if not flag_present(F, flag):
        return None
    out = {}
    for pos in nonunary(flag):
        S, E = frames(F, h, pos)
        V, hp = piece(F, h, pos, flag)
        res = rescale_piece(V, hp, S, E)
        if res is None:
            return None
        out[pos] = res
    return out


def h_parent(F, h, B):
    p = parent(F, B)
    return 0.0 if p is None else h[p]


def tree_frames(F, h, pos):
    return (h_parent(F, h, pos[1]), 1.0)


def labelled_frames(flag, e):
    """Example A.50(2): S = max(h(p(B)), e_{pi(nu)}), E = e_nu; e: {pos: height} on layers 0..l-2; e = 1 on layer l-1."""
    l = len(flag) - 1

    def E_of(pos):
        return 1.0 if pos[0] == l - 1 else e[pos]

    def fr(F, h, pos):
        pp = parent_position(flag, pos)
        c = 0.0 if pp is None else E_of(pp)
        return (max(h_parent(F, h, pos[1]), c), E_of(pos))
    return fr


def uniform_labelling(flag, a):
    """e_{(j,B)} = a_{j+1}, a = (a_1, ..., a_{l-1})."""
    return {(j, B): a[j] for (j, B) in positions(flag) if j <= len(flag) - 3}


def two_flag(F, alpha):
    T = root(F)
    return (frozenset([T]), alpha, frozenset(frozenset([x]) for x in T))


def c_tr(F, h, alpha):
    """Ching's cut at the partition alpha (Def A.31) = kappa with tree frames on (min, alpha, max)."""
    return kappa(F, h, two_flag(F, alpha), tree_frames)


def levelled_cut(F, h, alpha, a=0.5):
    flag = two_flag(F, alpha)
    return kappa(F, h, flag, labelled_frames(flag, uniform_labelling(flag, [a])))


def flatten(D):
    """A vertex of a piece is a set of blocks; flatten it to the underlying set of elements."""
    return frozenset().union(*D)
