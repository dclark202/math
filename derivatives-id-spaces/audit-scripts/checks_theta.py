"""Batch 1 checks: Lemma A.20, interval property/(A.20), Lemma A.23, Theorem A.24 (relations, bijectivity)."""
import random, sys
from partitions import *
from trees import *
random.seed(1)

def check_n(n, samples=3):
    T = frozenset(range(n))
    chains = nondegenerate_chains(T)
    st = dict(n=n, chains=len(chains), cells=0, degen=0, inner=0, outer=0, surj=0)
    for ch in chains:
        k = len(ch) - 1
        F = tree_of_chain(ch)
        assert is_tree(F, T)
        a, b = levels(ch)
        # interval property and (A.20)
        for B in F:
            assert all(B in ch[j] for j in range(a[B], b[B] + 1))
            if B != T:
                assert a[B] == b[parent(F, B)] + 1, ("A.20 fails", ch, B)
        for _ in range(samples):
            u = sorted(random.random() for _ in range(k))
            F, h = theta(ch, u)
            assert is_weighting(F, h) and is_reduced(F, h)                       # Lemma A.23(1)
            assert h[T] > 0 and all(h[leaf_parent(F, x)] < 1 for x in T)
            uu = [0.0] + u + [1.0]
            for j in range(k):                                                   # Lemma A.23(2)
                for t in (uu[j + 1], (uu[j] + uu[j + 1]) / 2):
                    assert alive(F, h, t) == ch[j], (ch, u, j, t)
            assert sorted(set(h.values())) == u                                  # Lemma A.23(3)
            ch2, u2 = theta_inverse(F, h)                                        # Thm A.24 (III)
            assert ch2 == ch and list(u2) == u
            st['cells'] += 1
            # degeneracies: Theta(s_j ch, u') == Theta(ch, sigma^j u'), u' in nabla^{k+1}
            up = sorted(random.random() for _ in range(k + 1))
            for j in range(k + 1):
                assert theta(degeneracy(ch, j), up) == theta(ch, up[:j] + up[j + 1:])
            st['degen'] += 1
            # inner faces: reduce(Theta(ch, delta^i u'')) == Theta(d_i ch, u''), u'' in nabla^{k-1}
            um = sorted(random.random() for _ in range(k - 1))
            for i in range(1, k):
                di = face(ch, i, T)
                assert di is not None
                delta = um[:i] + [um[i - 1]] + um[i:]
                assert reduce(*theta(ch, delta)) == theta(di, um), (ch, i)
            st['inner'] += 1
            # outer faces
            d0 = face(ch, 0, T); F0, h0 = theta(ch, [0.0] + um)
            if d0 is None: assert is_basepoint(F0, h0)
            else: assert (F0, h0) == theta(d0, um)
            dk = face(ch, k, T); Fk, hk = theta(ch, um + [1.0])
            if dk is None: assert is_basepoint(Fk, hk)
            else: assert (Fk, hk) == theta(dk, um)
            st['outer'] += 1
    # surjectivity / Lemma A.20 on random reduced weighted trees
    trees = all_trees(T)
    for F in trees:
        for _ in range(samples):
            h = random_weighting(F)
            assert is_reduced(F, h)
            ch, u = theta_inverse(F, h)
            assert is_chain(ch) and ch[0] == min_partition(T) and ch[-1] == max_partition(T)
            assert all(ch[j] != ch[j + 1] for j in range(len(ch) - 1))
            assert theta(ch, u) == (F, h)
            # Lemma A.20: alive(t) is a partition, refinement in t, jumps only at heights
            prev = None
            for t in sorted([random.random() for _ in range(6)] + list(h.values())):
                if t <= 0: continue
                A = alive(F, h, t)
                assert frozenset().union(*A) == T and sum(len(B) for B in A) == n
                if prev is not None: assert refines(A, prev)
                prev = A
            st['surj'] += 1
    st['trees'] = len(trees)
    return st

for n in range(2, 6):
    print(check_n(n))
print("ALL THETA CHECKS PASSED")
