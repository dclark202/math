"""Batch 4 checks on Def A.56 / Prop A.57 / Thm A.60 / Lemma A.63:
(W2),(W3),(W4) of composites; associativity of substitution; (A.26) cut(composite) = cut then cut;
o(theta) o o(theta') = o(merged); sigma respects composition (merged labelling)."""
import random
from partitions import *
from trees import *
from cuts import *
from wcs import *
random.seed(4)
EPS = 1e-9

def random_labelling(flag):
    """Strictly increasing along the tree, values in (0,1), on layers 0..l-2."""
    l = len(flag) - 1
    e = {}
    for j in range(l - 1):
        for B in flag[j]:
            pp = parent_position(flag, (j, B))
            lo = 0.0 if pp is None else e[pp]
            e[(j, B)] = lo + (1 - lo) * (0.2 + 0.6 * random.random())
    return e

def random_cell(flag, kind=None):
    kind = kind or random.choice(["tree", "const", "leg1", "leg2", "leg3"])
    if kind == "tree": return tree_cell(flag), kind
    c = constant_cell(flag, random_labelling(flag))
    if kind == "const": return c, kind
    if kind == "leg1": return leg1_cell(c), kind
    if kind == "leg2": return leg2_cell(c, random.random()), kind
    return leg3_cell(c, random.random()), kind

def check_axioms(cell, F, h, u, v, st):
    for pos in nonunary(cell.flag):
        S, E = cell.frames(F, h, pos, u, v)
        assert u - EPS <= S <= E + EPS and E <= v + EPS, ("W2", pos, S, E, u, v)
        a = anchor(cell.flag, pos)
        ok = S >= hp(F, h, pos[1], u) - EPS
        if not ok and a is not None:
            Sa, Ea = cell.frames(F, h, a, u, v)
            ok = Ea <= S + EPS
        if not ok:
            print("W3 FAILURE at", pos, "S,E =", (S, E), "anchor", a, "hp(B) =", hp(F, h, pos[1], u), "window", (u, v))
            print("  flag:", [sorted(map(sorted, lam)) for lam in cell.flag])
            print("  heights:", {tuple(sorted(B)): round(x, 4) for B, x in h.items()})
            for q in nonunary(cell.flag):
                print("   frames", q[0], sorted(q[1]), "->", tuple(round(t, 4) for t in cell.frames(F, h, q, u, v)), "hp", round(hp(F, h, q[1], u), 4), "anc", anchor(cell.flag, q))
            print("  DEBUG:", DEBUG)
            raise SystemExit(1)
        # (W4)
        aa, bb = 0.5 + random.random(), random.random() - 0.5
        h2 = {B: aa * x + bb for B, x in h.items()}
        S2, E2 = cell.frames(F, h2, pos, aa * u + bb, aa * v + bb)
        assert abs(S2 - (aa * S + bb)) < 1e-9 and abs(E2 - (aa * E + bb)) < 1e-9, ("W4", pos)
        st['axioms'] += 1

def deep(X):
    """Deep-flatten nested frozensets to the underlying set of elements."""
    out = set()
    for y in X:
        if isinstance(y, frozenset): out |= deep(y)
        else: out.add(y)
    return frozenset(out)

def check_A26(outer, j, inner, comp, F, h, st):
    of = outer.flag
    lp = next(len(c.flag) - 1 for c in inner.values())
    direct = kappa_cell(comp, F, h)
    first = kappa_cell(outer, F, h)
    if first is None:
        assert direct is None, "A.26: outer degenerate but composite not"
        st['A26_base'] += 1
        return
    # apply the inner cuts to the rescaled pieces
    assembled = {}
    degenerate = False
    for pos in nonunary(of):
        layer, B = pos
        if layer != j - 1:
            newpos = pos if layer <= j - 2 else (layer + lp - 1, B)
            assembled[newpos] = first[pos]
            continue
        V, hV = first[pos]
        sub = kappa_cell(inner[pos], V, hV)
        if sub is None:
            degenerate = True
            break
        for (r, Dp), (Vs, hs) in sub.items():
            assembled[(j - 1 + r, flatten(Dp))] = (Vs, hs)
    if degenerate:
        assert direct is None, "A.26: inner degenerate but composite not"
        st['A26_base'] += 1
        return
    assert direct is not None, "A.26: composite degenerate but factors not"
    assert set(direct) == set(assembled), (set(direct) ^ set(assembled))
    for pos in direct:
        hd = direct[pos][1]; ha = assembled[pos][1]
        keyd2 = {deep(X): val for X, val in hd.items()}
        keya2 = {deep(X): val for X, val in ha.items()}
        assert keyd2.keys() == keya2.keys(), (pos, keyd2.keys(), keya2.keys())
        for K in keyd2:
            assert abs(keyd2[K] - keya2[K]) < 1e-9, (pos, K, keyd2[K], keya2[K])
    st['A26'] += 1

def wrap_cell(cell):
    """Transport a cell along the atom bijection a -> {a} (the atoms of the flag get wrapped once)."""
    def wrap_set(X):
        return frozenset(frozenset([a]) for a in X)
    def unwrap_set(X):
        return frozenset(next(iter(a)) for a in X)
    flag = tuple(frozenset(wrap_set(D) for D in lam) for lam in cell.flag)
    def fn(F, h, pos, u, v):
        F0 = frozenset(unwrap_set(X) for X in F)
        h0 = {unwrap_set(X): val for X, val in h.items()}
        return cell.frames(F0, h0, (pos[0], unwrap_set(pos[1])), u, v)
    return Cell(flag, fn)

st = dict(axioms=0, A26=0, A26_base=0, assoc=0, o_comp=0, sigma_comp=0)
for n in (3, 4, 5):
    T = frozenset(range(n))
    for F in all_trees(T):
        h0 = random_weighting(F)
        chain, _ = theta_inverse(F, h0)
        k = len(chain) - 1
        if k < 3:
            continue
        for _ in range(3):
            # merged flag: random subchain containing min and max, length L in [3, k]
            L = random.randint(3, k)
            mids = sorted(random.sample(range(1, k), L - 1))
            merged = (chain[0],) + tuple(chain[i] for i in mids) + (chain[-1],)
            j = random.randint(1, L - 1)
            lp = random.randint(2, L - j + 1) if L - j + 1 >= 2 else 2
            lp = min(lp, L - j + 1)
            outer_flag = merged[:j] + merged[j + lp - 1:]
            outer, okind = random_cell(outer_flag)
            inner = {}; KINDS = {}
            for B in outer_flag[j - 1]:
                fl = inner_flag_at(merged, j, lp, (j - 1, B))
                inner[(j - 1, B)], kk = random_cell(fl); KINDS[tuple(sorted(B))] = kk
            comp = substitute(outer, j, inner)
            DEBUG = dict(okind=okind, j=j, lp=lp, outer_flag=[sorted(map(sorted, lam)) for lam in outer_flag], inner_kinds=KINDS)
            assert comp.flag == merged
            for _ in range(3):
                h = random_weighting(F)
                u, v = 0.0, 1.0
                check_axioms(comp, F, h, u, v, st)
                u2 = random.random() * 0.3; v2 = 0.7 + random.random() * 0.3
                hh = {B: u2 + (v2 - u2) * x for B, x in h.items()}
                check_axioms(comp, F, hh, u2, v2, st)
                check_A26(outer, j, inner, comp, F, h, st)
            # o and sigma respect composition
            if okind == "tree" and all(True for _ in inner):
                tinner = {rho: tree_cell(c.flag) for rho, c in inner.items()}
                tcomp = substitute(tree_cell(outer_flag), j, tinner)
                tmerged = tree_cell(merged)
                for _ in range(3):
                    h = random_weighting(F); u, v = random.random() * 0.3, 0.7 + random.random() * 0.3
                    for pos in nonunary(merged):
                        a1, a2 = tcomp.frames(F, h, pos, u, v), tmerged.frames(F, h, pos, u, v)
                        assert abs(a1[0] - a2[0]) < 1e-9 and abs(a1[1] - a2[1]) < 1e-9, ("o comp", pos, a1, a2)
                    st['o_comp'] += 1
            e_out = random_labelling(outer_flag)
            flags = {rho: c.flag for rho, c in inner.items()}
            labels = {rho: random_labelling(c.flag) for rho, c in inner.items()}
            cinner = {rho: constant_cell(flags[rho], labels[rho]) for rho in inner}
            ccomp = substitute(constant_cell(outer_flag, e_out), j, cinner)
            m2, e_m = merge_labellings(outer_flag, j, e_out, {'flags': flags, 'labels': labels})
            assert m2 == merged
            cm = constant_cell(merged, e_m)
            for _ in range(3):
                h = random_weighting(F); u, v = random.random() * 0.3, 0.7 + random.random() * 0.3
                for pos in nonunary(merged):
                    a1, a2 = ccomp.frames(F, h, pos, u, v), cm.frames(F, h, pos, u, v)
                    assert abs(a1[0] - a2[0]) < 1e-9 and abs(a1[1] - a2[1]) < 1e-9, ("sigma comp", pos, a1, a2)
                st['sigma_comp'] += 1
        # associativity: 4-level merged flag, two bracketings (A o_2 B) o_3 C  vs  A o_2 (B o_2 C)
        for _ in range(3):
            if k < 4:
                break
            mids = sorted(random.sample(range(1, k), 3))
            m4 = (chain[0],) + tuple(chain[i] for i in mids) + (chain[-1],)
            A_flag = (m4[0], m4[1], m4[4])
            A, _ = random_cell(A_flag)
            B_cells = {}
            for Bk in m4[1]:
                B_cells[(1, Bk)], _ = random_cell(inner_flag_at((m4[0], m4[1], m4[2], m4[4]), 2, 2, (1, Bk)))
            AB = substitute(A, 2, B_cells)
            assert AB.flag == (m4[0], m4[1], m4[2], m4[4])
            C_cells = {}
            for Ck in m4[2]:
                C_cells[(2, Ck)], _ = random_cell(inner_flag_at(m4, 3, 2, (2, Ck)))
            ABC1 = substitute(AB, 3, C_cells)
            BC = {}
            for Bk in m4[1]:
                Bc = B_cells[(1, Bk)]
                sub = {(1, Dp): wrap_cell(C_cells[(2, flatten(Dp))]) for Dp in Bc.flag[1]}
                BC[(1, Bk)] = substitute(Bc, 2, sub)
            ABC2 = substitute(A, 2, BC)
            assert ABC1.flag == m4 == ABC2.flag
            for _ in range(3):
                h = random_weighting(F); u, v = random.random() * 0.3, 0.7 + random.random() * 0.3
                for pos in nonunary(m4):
                    a1, a2 = ABC1.frames(F, h, pos, u, v), ABC2.frames(F, h, pos, u, v)
                    assert abs(a1[0] - a2[0]) < 1e-9 and abs(a1[1] - a2[1]) < 1e-9, ("assoc", pos, a1, a2)
                st['assoc'] += 1
print(st)
print("ALL WCS CHECKS PASSED")
