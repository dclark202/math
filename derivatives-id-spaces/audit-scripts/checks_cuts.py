"""Batch 3 checks: Def A.31 formulas; Example A.50(1) kappa_{F^tr} = c^tr; Example A.50(2) alive-partition
characterisation of the levelled cut; coassociativity of Ching's cuts (nested)."""
import random
from partitions import *
from trees import *
from cuts import *
random.seed(2)

EPS = 1e-12
stats = dict(formula=0, ex2=0, ex2_nonbase=0, coassoc=0)

for n in range(2, 6):
    T = frozenset(range(n))
    for F in all_trees(T):
        for _ in range(4):
            h = random_weighting(F)
            chain, u = theta_inverse(F, h)
            k = len(chain) - 1

            # --- Def A.31 explicit formulas (A.21)/(A.22) vs kappa with tree frames, alpha = any level of the chain
            for j in range(k + 1):
                alpha = chain[j]
                res = c_tr(F, h, alpha)
                blocks = sorted(alpha, key=sorted)
                assert res is not None, (F, h, alpha)          # alpha is present and no pi_i = 1 (h reduced, h(v(x)) < 1)
                if len(blocks) >= 2:
                    Vb, hb = res[(0, T)]
                    for D, v in hb.items():                       # bottom tree keeps heights, (A.21)
                        assert abs(v - h[flatten(D)]) < EPS
                for B in blocks:
                    if len(B) < 2:
                        continue
                    pi = h_parent(F, h, B)
                    V, hB = res[(1, B)]
                    for D, v in hB.items():                       # (A.22)
                        assert abs(v - (h[flatten(D)] - pi) / (1 - pi)) < EPS
                stats['formula'] += 1

            # --- Example A.50(2): levelled cut at a is non-basepoint iff alive(a) == alpha; pieces rescaled by [0,a],[a,1]
            for _ in range(4):
                j = random.randrange(k + 1)
                alpha = chain[j]
                a = random.random()
                res = levelled_cut(F, h, alpha, a)
                predicted = (alive(F, h, a) == alpha)
                assert (res is not None) == predicted, (F, h, alpha, a, res, predicted)
                stats['ex2'] += 1
                if res is not None:
                    stats['ex2_nonbase'] += 1
                    if len(alpha) >= 2:
                        for D, v in res[(0, T)][1].items():
                            assert abs(v - h[flatten(D)] / a) < EPS
                    for B in alpha:
                        if len(B) < 2:
                            continue
                        for D, v in res[(1, B)][1].items():
                            assert abs(v - (h[flatten(D)] - a) / (1 - a)) < EPS

            # --- coassociativity: cut at chain[i], then cut each top piece at chain[j] restricted, vs the 3-flag tree cut
            if k >= 2:
                for _ in range(3):
                    i = random.randrange(1, k)
                    j = random.randrange(i, k + 1)
                    flag3 = (chain[0], chain[i], chain[j], chain[-1])
                    direct = kappa(F, h, flag3, tree_frames)
                    first = c_tr(F, h, chain[i])
                    assert direct is not None and first is not None
                    for B in chain[i]:
                        if len(B) < 2:
                            continue
                        V, hB = first[(1, B)]                     # weighted tree on the singletons of B
                        beta = frozenset(frozenset(frozenset([x]) for x in C) for C in chain[j] if C <= B)
                        sub = c_tr(V, hB, beta)
                        assert sub is not None
                        if len(beta) >= 2:                        # bottom of the sub-cut = direct piece at (1, B)
                            hd = direct[(1, B)][1]
                            hs = sub[(0, root(V))][1]
                            assert all(abs(hs[D] - hd[frozenset(flatten(C) for C in D)]) < EPS for D in hs)
                        for C in beta:                            # top pieces of the sub-cut = direct pieces at (2, C)
                            if len(C) < 2:
                                continue
                            hd = direct[(2, flatten(C))][1]
                            hs = sub[(1, C)][1]
                            assert all(abs(hs[D] - hd[frozenset(flatten(x) for x in D)]) < EPS for D in hs), (hs, hd)
                    stats['coassoc'] += 1

print(stats)
print("ALL CUT CHECKS PASSED")
