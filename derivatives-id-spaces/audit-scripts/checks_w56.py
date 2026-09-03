"""Check axioms (W5), (W6) with the EXPLICIT constants derived in the repaired proofs:
tree cells (λ=1, K=0); constant cells (λ_ν = e_ν − e_π, K=0); composites (seam/inner root-run: λ'' = λ_ρ λ'_root,
K'' = (K_ρ+1)/λ'_root; insert positions: inner constants; unrefined: outer constants); Legs 1–2: λ_s = min(λ, 1/(1/λ+K), 1),
K_s = K; Leg 3 at t>0: λ = t, K = 0."""
import random, math
from partitions import *
from trees import *
from cuts import *
from wcs import *
random.seed(21)
EPS = 1e-9

def random_labelling(flag):
    l = len(flag) - 1; e = {}
    for j in range(l - 1):
        for B in flag[j]:
            pp = parent_position(flag, (j, B)); lo = 0.0 if pp is None else e[pp]
            e[(j, B)] = lo + (1 - lo) * (0.2 + 0.6 * random.random())
    return e

def tree_with_consts(flag):
    return tree_cell(flag), {p: 1.0 for p in nonunary(flag)}, {p: 0.0 for p in nonunary(flag)}

def const_with_consts(flag):
    e = random_labelling(flag); l = len(flag) - 1
    def E_of(pos): return 1.0 if pos[0] == l - 1 else e[pos]
    lam = {}
    for p in nonunary(flag):
        pp = parent_position(flag, p); s = 0.0 if pp is None else E_of(pp)
        lam[p] = E_of(p) - s
    return constant_cell(flag, e), lam, {p: 0.0 for p in nonunary(flag)}

def composite_with_consts(outer, inner, j):
    """outer = (cell, lam, K); inner = {rho: (cell, lam, K)}."""
    oc, olam, oK = outer
    of = oc.flag
    icells = {rho: t[0] for rho, t in inner.items()}
    comp = substitute(oc, j, icells)
    merged, lp = merged_flag(of, j, icells)
    lam, K = {}, {}
    for pos in nonunary(merged):
        layer, B = pos
        if layer <= j - 2:
            lam[pos], K[pos] = olam[pos], oK[pos]
        elif layer >= j - 1 + lp:
            q = (layer - (lp - 1), B); lam[pos], K[pos] = olam[q], oK[q]
        else:
            r = layer - (j - 1)
            B0 = next(C for C in of[j - 1] if B <= C); rho = (j - 1, B0)
            ic, ilam, iK = inner[rho]
            Dp = next(D for D in ic.flag[r] if flatten(D) == B)
            ipos = (r, Dp)
            if B == B0:      # seam / inner root-run position: inner constant at that inner root-run position
                lam[pos] = olam[rho] * ilam[ipos]
                K[pos] = (oK[rho] + 1) / ilam[ipos]
            else:
                lam[pos], K[pos] = ilam[ipos], iK[ipos]
    return comp, lam, K

def legs_with_consts(cellc, s, kind):
    c, lam, K = cellc
    if kind == "leg1":
        return leg1_cell(c), {p: min(lam[p], 1 / (1 / lam[p] + K[p]), 1.0) for p in lam}, dict(K)
    if kind == "leg2":
        return leg2_cell(c, s), {p: min(lam[p], 1 / (1 / lam[p] + K[p]), 1.0) for p in lam}, dict(K)
    return leg3_cell(c, s), {p: s for p in lam}, {p: 0.0 for p in lam}   # s > 0

def check(cellc, F, h, u, v, st, name):
    c, lam, K = cellc
    for pos in nonunary(c.flag):
        S, E = c.frames(F, h, pos, u, v); w = E - S
        a = anchor(c.flag, pos)
        if a is None:
            assert w >= lam[pos] * (v - u) - EPS, (name, "W5 root", pos, w, lam[pos], v - u)
        else:
            Sa, Ea = c.frames(F, h, a, u, v)
            rhs = lam[pos] * min(v - u, E - hp(F, h, pos[1], u), Ea - Sa)
            assert w >= rhs - EPS, (name, "W5", pos, w, rhs)
            assert Ea - E <= K[pos] * w + EPS, (name, "W6", pos, Ea - E, K[pos], w)
        st[name] = st.get(name, 0) + 1

def squeeze(F, h):
    h = dict(h)
    for B in internal(F):
        r = random.random(); p = parent(F, B)
        if r < 0.25 and p is not None: h[B] = h[p] + 10 ** (-random.randint(1, 8)) * (h[B] - h[p])
        elif r < 0.4: h[B] = 1 - 10 ** (-random.randint(1, 8)) * (1 - h[B])
    def fix(B, lo):
        if len(B) == 1: return
        h[B] = max(h[B], lo)
        for C in children(F, B): fix(C, h[B])
    fix(root(F), 0.0); return h

st = {}
for n in (3, 4, 5):
    T = frozenset(range(n))
    for F in all_trees(T):
        h0 = random_weighting(F); chain, _ = theta_inverse(F, h0); k = len(chain) - 1
        if k < 3: continue
        for _ in range(2):
            L = random.randint(3, k); mids = sorted(random.sample(range(1, k), L - 1))
            merged = (chain[0],) + tuple(chain[i] for i in mids) + (chain[-1],)
            j = random.randint(1, L - 1); lp = min(random.randint(2, L), L - j + 1)
            of = merged[:j] + merged[j + lp - 1:]
            mk = lambda kind, fl: tree_with_consts(fl) if kind == "tree" else const_with_consts(fl)
            cells = {}
            for ko in ("tree", "const"):
                for ki in ("tree", "const"):
                    inner = {(j - 1, B): mk(ki, inner_flag_at(merged, j, lp, (j - 1, B))) for B in of[j - 1]}
                    cells[f"{ko}*{ki}"] = composite_with_consts(mk(ko, of), inner, j)
            cells["tree"] = mk("tree", merged); cells["const"] = mk("const", merged)
            base = dict(cells)
            for name, cc in base.items():
                for kind in ("leg1", "leg2", "leg3"):
                    s = 0.05 + 0.9 * random.random()
                    cells[f"{kind}({name})"] = legs_with_consts(cc, s, kind)
            # second-order composites: leg cells substituted into leg cells
            inner = {(j - 1, B): legs_with_consts(mk("const", inner_flag_at(merged, j, lp, (j - 1, B))), 0.3, "leg2") for B in of[j - 1]}
            cells["leg2(const)*leg2(const)"] = composite_with_consts(legs_with_consts(mk("const", of), 0.6, "leg2"), inner, j)
            for name, cc in cells.items():
                for _ in range(5):
                    h = random_weighting(F)
                    if random.random() < 0.7: h = squeeze(F, h)
                    u, v = (0.0, 1.0) if random.random() < 0.5 else (random.random() * 0.3, 0.7 + random.random() * 0.3)
                    hh = {B: u + (v - u) * x for B, x in h.items()}
                    check(cc, F, hh, u, v, st, name)
print({k: v for k, v in sorted(st.items())})
print("ALL (W5)/(W6) CONSTANT CHECKS PASSED")
