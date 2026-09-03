"""Windowed cut systems (Def A.54), tree cells and constant cells (Example A.55), substitution (Def A.56, (A.25)),
anchors (Def A.52), cut maps of cells, and the L-operad composition of labellings (Def A.62)."""
from cuts import positions, child_set, nonunary, parent_position, flag_present, piece, rescale_piece, flatten, kappa
from trees import parent, root


def hp(F, h, B, u):
    """Parent height with the convention h(p(root)) := u (Def A.52)."""
    p = parent(F, B)
    return u if p is None else h[p]


def clamp(x, u, v):
    return min(max(x, u), v)


class Cell:
    """A windowed cut system of shape `flag`: fn(F, h, pos, u, v) -> (S, E) at non-unary positions."""
    def __init__(self, flag, fn):
        self.flag, self.fn = flag, fn

    def frames(self, F, h, pos, u, v):
        return self.fn(F, h, pos, u, v)


def tree_cell(flag):
    return Cell(flag, lambda F, h, pos, u, v: (clamp(hp(F, h, pos[1], u), u, v), v))


def constant_cell(flag, e):
    """e: {pos: label} on layers 0..l-2; label 1 on layer l-1; label 0 at the parent of the root position."""
    l = len(flag) - 1

    def E_of(pos):
        return 1.0 if pos[0] == l - 1 else e[pos]

    def fn(F, h, pos, u, v):
        pp = parent_position(flag, pos)
        s = 0.0 if pp is None else E_of(pp)
        return (u + (v - u) * s, u + (v - u) * E_of(pos))
    return Cell(flag, fn)


def mix_cell(c1, c2, s):
    assert c1.flag == c2.flag

    def fn(F, h, pos, u, v):
        (S1, E1), (S2, E2) = c1.frames(F, h, pos, u, v), c2.frames(F, h, pos, u, v)
        return ((1 - s) * S1 + s * S2, (1 - s) * E1 + s * E2)
    return Cell(c1.flag, fn)


def merged_flag(of, j, inner):
    lp = next(len(c.flag) - 1 for c in inner.values())
    merged = list(of[:j])
    for r in range(1, lp):
        mu = set()
        for rho, c in inner.items():
            for D in c.flag[r]:
                mu.add(flatten(D))
        merged.append(frozenset(mu))
    merged += list(of[j:])
    return tuple(merged), lp


def substitute(outer, j, inner):
    """Composite of `outer` with the family `inner` (one cell per position rho = (j-1, B) of the outer flag, all of the
    same level) at level slot j (1-indexed), by (A.25): frames of the outer cell at unrefined positions; at an insert
    position, the inner frame evaluated on the raw piece at rho with the handed window (S_rho, E_rho)."""
    of = outer.flag
    merged, lp = merged_flag(of, j, inner)

    def fn(F, h, pos, u, v):
        layer, B = pos
        if layer <= j - 2:
            return outer.frames(F, h, pos, u, v)
        if layer >= j - 1 + lp:
            return outer.frames(F, h, (layer - (lp - 1), B), u, v)
        r = layer - (j - 1)
        B0 = next(C for C in of[j - 1] if B <= C)
        rho = (j - 1, B0)
        c = inner[rho]
        S_rho, E_rho = outer.frames(F, h, rho, u, v)
        V, hraw = piece(F, h, rho, of)
        Dp = next(D for D in c.flag[r] if flatten(D) == B)
        return c.frames(V, hraw, (r, Dp), S_rho, E_rho)
    return Cell(merged, fn)


def anchor(flag, pos):
    """Def A.52: None for a root-run position, else (j0 - 1, block of lambda_{j0-1} containing B)."""
    j, B = pos
    j0 = j
    while j0 > 0 and B in flag[j0 - 1]:
        j0 -= 1
    if j0 == 0:
        return None
    return (j0 - 1, next(C for C in flag[j0 - 1] if B <= C))


def kappa_cell(cell, F, h):
    """Cut map of the top-window slice (Def A.54)."""
    return kappa(F, h, cell.flag, lambda F_, h_, pos: cell.frames(F_, h_, pos, 0.0, 1.0))


def inner_flag_at(merged, j, lp, rho):
    """The inner flag at rho = (j-1, B): partitions of the children set, read off the merged flag."""
    _, B = rho
    Cs = child_set(merged[:j] + merged[j + lp - 1:], rho)  # children of B in the outer flag
    out = []
    for r in range(lp + 1):
        lam = merged[j - 1 + r]
        out.append(frozenset(frozenset(C for C in Cs if C <= D) for D in lam if D <= B))
    return tuple(out)


def merge_labellings(of, j, e_out, inner_e):
    """L-composition (Def A.62): keep outer labels, insert inner labels affinely into the window of their position."""
    merged, lp = merged_flag(of, j, {rho: Cell(fl, None) for rho, fl in inner_e['flags'].items()})
    l = len(of) - 1
    e = {}
    for pos in positions(merged):
        layer, B = pos
        if layer >= len(merged) - 2 + 1:
            continue
        if layer <= j - 2:
            e[pos] = e_out[pos]
        elif layer >= j - 1 + lp:
            e[pos] = e_out[(layer - (lp - 1), B)] if (layer - (lp - 1)) <= l - 2 else None
        else:
            r = layer - (j - 1)
            B0 = next(C for C in of[j - 1] if B <= C)
            rho = (j - 1, B0)
            pp = parent_position(of, rho)
            s = 0.0 if pp is None else (e_out[pp] if pp[0] <= l - 2 else 1.0)
            t = e_out[rho] if rho[0] <= l - 2 else 1.0
            fl = inner_e['flags'][rho]
            Dp = next(D for D in fl[r] if flatten(D) == B)
            lab = inner_e['labels'][rho]
            b = lab[(r, Dp)] if r <= len(fl) - 3 else 1.0
            e[pos] = s + (t - s) * b
    return merged, {k: v for k, v in e.items() if v is not None}


def anchored_labelled_cell(flag, e):
    """Windowed version of the labelled frames F^e of Example A.50(2): S = clamp(max(h(p(B)), u+(v-u)e_pi)), E = u+(v-u)e_nu.
    Anchored, hence (W3); (W1),(W2),(W4) clear."""
    l = len(flag) - 1

    def E_of(pos):
        return 1.0 if pos[0] == l - 1 else e[pos]

    def fn(F, h, pos, u, v):
        pp = parent_position(flag, pos)
        s = 0.0 if pp is None else E_of(pp)
        return (clamp(max(hp(F, h, pos[1], u), u + (v - u) * s), u, v), u + (v - u) * E_of(pos))
    return Cell(flag, fn)


def leg1_cell(cell):
    """Leg 1 of Prop A.57: S_flat = min(max(S, h(p(B))), E), ceilings unchanged."""
    def fn(F, h, pos, u, v):
        S, E = cell.frames(F, h, pos, u, v)
        return (min(max(S, hp(F, h, pos[1], u)), E), E)
    return Cell(cell.flag, fn)


def leg2_cell(cell, s):
    """Leg 2 of Prop A.57: E_s = (1-s)E + s v, S_s = min(max(S, h(p(B))), E_s)."""
    def fn(F, h, pos, u, v):
        S, E = cell.frames(F, h, pos, u, v)
        Es = (1 - s) * E + s * v
        return (min(max(S, hp(F, h, pos[1], u)), Es), Es)
    return Cell(cell.flag, fn)


def leg3_cell(cell, s):
    """Leg 3 of Prop A.57: straight line from the Leg-2 endpoint to the tree cell."""
    end2 = leg2_cell(cell, 1.0)
    tr = tree_cell(cell.flag)
    return mix_cell(end2, tr, s)
