"""Batch 5 checks: Prop A.71 (transport). Chain side: the levelled action of the multi-cut cells (Lemma 4.4 decomposition
iterated + the rescalings of Construction A.69), and for a general labelling the peeled composite of Def A.70;
tree side: kappa_{F^e} composed with Theta.  Also Lemma A.64 (kappa_{W^e} = kappa_{F^e}) and the coface compatibility of
Construction A.69."""
import random
from partitions import *
from trees import *
from cuts import *
from wcs import constant_cell, kappa_cell
random.seed(5)
EPS = 1e-9


def collapse_segment(segment, B, children):
    """Partitions of a segment, restricted to B and collapsed at the children (each block -> set of children in it)."""
    out = []
    for lam in segment:
        out.append(frozenset(frozenset(C for C in children if C <= D) for D in lam if D <= B))
    return tuple(out)


def piece_of_chain(chain, heights, flag_present_test=None):
    """Theta of a (possibly degenerate) chain with heights, reduced; None if basepoint."""
    F, h = theta(chain, heights)
    F, h = reduce(F, h)
    if not h:
        return (F, h)
    if is_basepoint(F, h):
        return None
    return (F, h)


def chain_cut_general(chain, u, flag, e, layer0=0):
    """Levelled cut of the cell (chain, u) along the flag with labelling e (labels on positions (j,B), j <= l-2;
    label 1 on the top layer), by peeling from the root (Def A.70). Returns {pos: (V,h)} with pos relative to layer0."""
    m = len(chain) - 1
    l = len(flag) - 1
    T = next(iter(flag[0]))
    if l == 1:
        # single position (0, T): the whole chain collapsed at the children (= singletons) with heights u
        children = flag[1]
        seg = collapse_segment(chain, T, children)
        pc = piece_of_chain(seg, u)
        return None if pc is None else {(layer0, T): pc}
    a = e[(0, T)]
    uu = [0.0] + list(u) + [1.0]
    P = max(i for i in range(m + 1) if uu[i] <= a)
    if chain[P] != flag[1]:
        return None
    children = flag[1]
    front = collapse_segment(chain[:P + 1], T, children)
    pc = piece_of_chain(front, [x / a for x in u[:P]])
    if pc is None:
        return None
    out = {(layer0, T): pc}
    for B in flag[1]:
        subflag = tuple(frozenset(D for D in lam if D <= B) for lam in flag[1:])
        if len(B) == 1:
            continue
        sub_e = {}
        for (j, D), val in e.items():
            if j >= 1 and D <= B:
                sub_e[(j - 1, D)] = (val - a) / (1 - a)
        seg = tuple(frozenset(D for D in lam if D <= B) for lam in chain[P:])
        hs = [(x - a) / (1 - a) for x in u[P:]]
        sub = chain_cut_general(seg, hs, subflag, sub_e, layer0 + 1)
        if sub is None:
            return None
        out.update(sub)
    return out


def canon(res):
    """Normalise a dict pos -> (V, h) to {pos: {element-set: height}} after reduction."""
    out = {}
    for pos, (V, h) in res.items():
        V2, h2 = reduce(V, h)
        if not h2:
            continue   # unary position (trivial piece): not a factor of B(gamma)
        out[(pos[0], deep(pos[1]))] = {deep(D): round(v, 9) for D, v in h2.items()}
    return out


def deep(X):
    if isinstance(X, frozenset):
        s = set()
        for y in X:
            s |= deep(y) if isinstance(y, frozenset) else {y}
        return frozenset(s)
    return frozenset([X])


def random_labelling(flag):
    l = len(flag) - 1
    e = {}
    for j in range(l - 1):
        for B in flag[j]:
            pp = parent_position(flag, (j, B))
            lo = 0.0 if pp is None else e[pp]
            e[(j, B)] = lo + (1 - lo) * (0.2 + 0.6 * random.random())
    return e


st = dict(uniform=0, general=0, nonbase=0, unrelated=0, degenerate=0, A64=0, tie=0)
for t in (2, 3, 4, 5):
    T = frozenset(range(t))
    chains = nondegenerate_chains(T)
    for trial in range(400 if t == 5 else 150):
        gamma = random.choice(chains)
        m = len(gamma) - 1
        u = sorted(random.random() for _ in range(m))
        # flag: (i) random subchain of gamma, (ii) unrelated random chain, (iii) degenerate gamma
        mode = random.choice(["sub", "sub", "unrelated", "degenerate"])
        if mode == "degenerate":
            jj = random.randrange(m + 1)
            gamma_d = degeneracy(gamma, jj)
            u_d = sorted(u + [random.random()])
            gamma_use, u_use = gamma_d, u_d
            mode = "sub"
            st['degenerate'] += 1
        else:
            gamma_use, u_use = gamma, u
        if mode == "sub":
            L = random.randint(1, m) if m >= 1 else 1
            mids = sorted(random.sample(range(1, m), min(L - 1, max(0, m - 1))))
            flag = (gamma[0],) + tuple(gamma[i] for i in mids) + (gamma[-1],)
        else:
            flag = random.choice(chains)
            st['unrelated'] += 1
        l = len(flag) - 1
        F, h = theta(gamma_use, u_use)
        F, h = reduce(F, h)
        # uniform labelling
        a = sorted(random.random() for _ in range(l - 1))
        e_u = {(j, B): a[j] for (j, B) in positions(flag) if j <= l - 2}
        tree_side = kappa(F, h, flag, labelled_frames(flag, e_u))
        chain_side = chain_cut_general(gamma_use, u_use, flag, e_u)
        assert (tree_side is None) == (chain_side is None), ("uniform", gamma_use, u_use, flag, a)
        if tree_side is not None:
            assert canon(tree_side) == canon(chain_side), ("uniform values", canon(tree_side), canon(chain_side))
            st['nonbase'] += 1
        st['uniform'] += 1
        # general labelling
        e_g = random_labelling(flag)
        tree_side = kappa(F, h, flag, labelled_frames(flag, e_g))
        chain_side = chain_cut_general(gamma_use, u_use, flag, e_g)
        assert (tree_side is None) == (chain_side is None), ("general", gamma_use, u_use, flag, e_g)
        if tree_side is not None:
            assert canon(tree_side) == canon(chain_side), ("general values",)
        st['general'] += 1
        # Lemma A.64: constant cell (stacked) vs labelled frames (anchored) have the same cut map
        w = kappa_cell(constant_cell(flag, e_g), F, h)
        assert (w is None) == (tree_side is None)
        if w is not None:
            assert canon(w) == canon(tree_side)
        st['A64'] += 1
        # ties: force a height equal to a cut height
        if l >= 2 and m >= 1:
            i = random.randrange(m)
            u_t = list(u_use); u_t[i] = a[0] if l == 2 else e_g[(0, T)]
            u_t = sorted(u_t)
            Ft, ht = theta(gamma_use, u_t); Ft, ht = reduce(Ft, ht)
            for lab in (e_u, e_g):
                ts = kappa(Ft, ht, flag, labelled_frames(flag, lab))
                cs = chain_cut_general(gamma_use, u_t, flag, lab)
                assert (ts is None) == (cs is None), ("tie", gamma_use, u_t, flag, lab)
                if ts is not None and canon(ts) != canon(cs):
                    print("TIE MISMATCH")
                    print(" chain:", [sorted(map(sorted, lam)) for lam in gamma_use])
                    print(" u_t:", [round(x, 4) for x in u_t])
                    print(" flag:", [sorted(map(sorted, lam)) for lam in flag])
                    print(" lab:", {(j, tuple(sorted(B))): round(v, 4) for (j, B), v in lab.items()})
                    print(" tree side :", {(p, tuple(sorted(B))): {tuple(sorted(D)): v for D, v in d.items()} for (p, B), d in canon(ts).items()})
                    print(" chain side:", {(p, tuple(sorted(B))): {tuple(sorted(D)): v for D, v in d.items()} for (p, B), d in canon(cs).items()})
                    raise SystemExit(1)
            st['tie'] += 1
print(st)

# Construction A.69: coface compatibility of the multi-cut map chi_a on nabla^m (heights coordinates)
def chi(u, a):
    """Blockwise rescaled tuples for the multidegree determined by u and a (ties: u_i == a_j goes to the earlier block)."""
    uu = [0.0] + list(u) + [1.0]
    m = len(u)
    aa = [0.0] + list(a) + [1.0]
    Ps = [0] + [max(i for i in range(m + 1) if uu[i] <= aj) for aj in a] + [m]
    return tuple(tuple((uu[i] - aa[j - 1]) / (aa[j] - aa[j - 1]) for i in range(Ps[j - 1] + 1, Ps[j] + 1)) for j in range(1, len(aa)))

def coface(u, i):
    m = len(u)
    if i == 0: return [0.0] + list(u)
    if i == m: return list(u) + [1.0]
    return list(u[:i]) + [u[i - 1]] + list(u[i:])

cnt = 0
for _ in range(3000):
    k = random.randint(0, 6); l = random.randint(2, 4)
    u = sorted(random.random() for _ in range(k))     # u in nabla^k
    a = sorted(random.random() for _ in range(l - 1))
    base = chi(u, a)
    aa = [0.0] + list(a) + [1.0]
    for i in range(k + 1):                            # cofaces delta^i : nabla^k -> nabla^{k+1}
        img = chi(coface(u, i), a)
        val = 0.0 if i == 0 else (1.0 if i == k else u[i - 1])
        b = max(j for j in range(l) if aa[j] <= val) if val < 1.0 else l - 1
        exp = list(map(list, base))
        r = (val - aa[b]) / (aa[b + 1] - aa[b])
        blk = exp[b]
        if i == 0: blk.insert(0, 0.0)
        elif i == k: blk.append(1.0)
        else: blk.insert(sum(1 for x in blk if x < r - 1e-12), r)
        ok = len(img) == len(exp) and all(len(x) == len(y) and all(abs(p - q) < 1e-9 for p, q in zip(x, y)) for x, y in zip(img, exp))
        assert ok, (u, a, i, img, exp)
        cnt += 1
print("coface compatibility cases:", cnt)
print("ALL TRANSPORT CHECKS PASSED")
