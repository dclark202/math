"""Do the cells along the three legs of Prop A.57 (started from constant cells) have continuous cut maps?
We approach (i) flag-absence walls h(B) -> h(p(B)), (ii) the basepoint (root height -> 0, a maximal vertex -> 1),
and measure the distance-to-basepoint d of the value; continuity needs d -> 0 (or exact basepoint)."""
import random
from partitions import *
from trees import *
from cuts import *
from wcs import *
random.seed(7)

def dist(res):
    if res is None:
        return 0.0
    d = 1.0
    for pos, (V, hn) in res.items():
        r = root(V)
        d = min(d, hn[r], min(1 - hn[leaf_parent(V, x)] for x in r))
    return d

def random_labelling(flag):
    l = len(flag) - 1
    e = {}
    for j in range(l - 1):
        for B in flag[j]:
            pp = parent_position(flag, (j, B))
            lo = 0.0 if pp is None else e[pp]
            e[(j, B)] = lo + (1 - lo) * (0.2 + 0.6 * random.random())
    return e

EPSS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
bad = {}
tested = {}
for n in (3, 4, 5):
    T = frozenset(range(n))
    for F in all_trees(T):
        for _ in range(3):
            h = random_weighting(F)
            chain, _ = theta_inverse(F, h)
            k = len(chain) - 1
            if k < 2:
                continue
            L = random.randint(2, k)
            mids = sorted(random.sample(range(1, k), L - 1))
            flag = (chain[0],) + tuple(chain[i] for i in mids) + (chain[-1],)
            e = random_labelling(flag)
            c = constant_cell(flag, e)
            s1, s2 = random.random(), random.random()
            cells = {"tree": tree_cell(flag), "const": c, "leg1": leg1_cell(c), f"leg2": leg2_cell(c, s1), f"leg3": leg3_cell(c, s2),
                     "tree*leg1": None}
            # a composite: tree outer cell refined by leg1 cells at slot 1 (root), when the flag has >= 3 levels
            cells.pop("tree*leg1")
            if L >= 3:
                j = 1
                outer_flag = flag[:j] + flag[j + 1:]
                inner = {}
                for B in outer_flag[j - 1]:
                    fl = inner_flag_at(flag, j, 2, (j - 1, B))
                    inner[(j - 1, B)] = leg1_cell(constant_cell(fl, random_labelling(fl)))
                cells["tree*leg1"] = substitute(tree_cell(outer_flag), j, inner)
                inner2 = {rho: leg2_cell(constant_cell(cc.flag, random_labelling(cc.flag)), random.random()) for rho, cc in inner.items()}
                cells["leg1*leg2"] = substitute(leg1_cell(constant_cell(outer_flag, random_labelling(outer_flag))), j, inner2)
            # sequences
            seqs = []
            flag_blocks = {B for lam in flag for B in lam if len(B) > 1}
            for B in internal(F):
                if B == T or B not in flag_blocks:
                    continue
                p = parent(F, B)
                seqs.append(("wall " + str(sorted(B)), [{**h, B: h[p] + eps * (h[B] - h[p])} for eps in EPSS + [0.0]]))
            seqs.append(("root->0", [{**h, T: eps * h[T]} for eps in EPSS + [0.0]]))
            top = max(internal(F), key=lambda B: h[B])
            seqs.append(("top->1", [{**h, top: 1 - eps * (1 - h[top])} for eps in EPSS + [0.0]]))
            for name, cell in cells.items():
                for sname, hs in seqs:
                    ds = [dist(kappa_cell(cell, F, hh)) for hh in hs]
                    tested[name] = tested.get(name, 0) + 1
                    if ds[-1] > 1e-9:
                        print("NOT BASEPOINT AT LIMIT:", name, sname, "n =", n, "ds =", [round(d, 4) for d in ds])
                        print("  flag:", [sorted(map(sorted, lam)) for lam in flag])
                        hh = hs[-1]
                        print("  heights at limit:", {tuple(sorted(B)): round(x, 4) for B, x in hh.items()})
                        for q in nonunary(cell.flag):
                            print("   pos", q[0], sorted(q[1]), "frames", tuple(round(t, 4) for t in cell.frames(F, hh, q, 0.0, 1.0)), "hp", round(hp(F, hh, q[1], 0.0), 4))
                        res = kappa_cell(cell, F, hh)
                        print("  value:", {(q[0], tuple(sorted(q[1]))): {tuple(sorted(map(tuple, map(sorted, D)))): round(v, 4) for D, v in hn.items()} for q, (V, hn) in res.items()})
                        raise SystemExit(1)
                    if ds[-2] > 1e-3 and ds[-3] > 1e-3:
                        bad.setdefault(name, []).append((n, sname, [round(d, 4) for d in ds]))
print("sequences tested per cell type:", tested)
for name in tested:
    b = bad.get(name, [])
    print(f"{name:10s}: {len(b)} discontinuous sequences" + (f"; e.g. {b[0]}" if b else ""))
