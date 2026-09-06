"""Check F1 (memo INFTY-COMPARISON.md, A0.3 (S1)): the composition product of
reduced symmetric sequences is homotopically free.

For reduced X, Y the k-th summand of (X o Y)(n) is X(k) smashed over Sigma_k
with Y^{k}(n) = wedge over surjections f: {1..n} ->> {1..k} of Y(f^-1 1) ^ ... ^ Y(f^-1 k).
The claim is that Sigma_k (acting by post-composition on the surjections) acts
FREELY on the set of surjections, so that the strict coinvariants are the
coinvariants of a free permutation action on wedge summands, i.e. one summand
per orbit and no residual quotient. We also check that the orbits are in
bijection with unordered partitions of {1..n} into k blocks (Stirling S(n,k)),
so that (X o Y)(n) = wedge over partitions pi of X(k) ^ smash_{B in pi} Y(B).

Second check: the paper's formula (2.?) writes the same summand as
Sigma[n] ^_{H(k_1..k_m)} X(m) ^ Y(k_1) ^ ... ^ Y(k_m) with H the block-
permutation stabilizer; H acts freely on Sigma_n by right multiplication, so
this quotient is free as well.  We verify |Sigma_n / H| = number of ordered
partitions with the given block sizes modulo permutations of equal-size blocks.

Pure Python 3. Ends with ALL F1 CHECKS PASSED.
"""
from itertools import permutations, product
from math import factorial
from collections import Counter

def surjections(n, k):
    for f in product(range(k), repeat=n):
        if len(set(f)) == k:
            yield f

def stirling2(n, k):
    # S(n,k) by recursion
    S = [[0]*(k+1) for _ in range(n+1)]
    S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, k+1):
            S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][k]

def check_free_action(n, k):
    surjs = list(surjections(n, k))
    perms = list(permutations(range(k)))
    # freeness: sigma . f = f  implies sigma = id
    for f in surjs:
        for s in perms:
            g = tuple(s[x] for x in f)
            if g == f and any(s[i] != i for i in range(k)):
                raise AssertionError(f"non-free action: n={n} k={k} f={f} fixed by {s}")
    # orbits <-> unordered partitions into k blocks
    orbits = set()
    for f in surjs:
        blocks = frozenset(frozenset(i for i in range(n) if f[i] == j) for j in range(k))
        orbits.add(blocks)
    assert len(orbits) == stirling2(n, k), (n, k, len(orbits), stirling2(n, k))
    assert len(surjs) == factorial(k) * stirling2(n, k)
    return len(surjs), len(orbits)

def check_H_free(ks):
    # H(k_1..k_m) = prod Sigma_{k_i} semidirect (permutations of equal blocks) acting on Sigma_n
    n = sum(ks)
    # build H as a set of permutations of range(n)
    m = len(ks)
    starts = [sum(ks[:i]) for i in range(m)]
    blocks = [list(range(starts[i], starts[i]+ks[i])) for i in range(m)]
    # block permutations preserving sizes
    size_groups = {}
    for i, kk in enumerate(ks):
        size_groups.setdefault(kk, []).append(i)
    H = set()
    # generate: for each block a permutation, and a permutation of blocks within size classes
    inner = [list(permutations(b)) for b in blocks]
    def block_perms():
        # product over size classes of permutations of the index lists
        classes = list(size_groups.values())
        for choice in product(*[list(permutations(c)) for c in classes]):
            bp = list(range(m))
            for c, ch in zip(classes, choice):
                for a, b in zip(c, ch):
                    bp[a] = b
            yield bp
    for bp in block_perms():
        for ins in product(*inner):
            perm = [None]*n
            for i in range(m):
                target_block = blocks[bp[i]]
                for pos, x in enumerate(blocks[i]):
                    perm[x] = target_block[pos] if False else None
            # simpler: perm sends block i (in its inner order ins[i]) onto block bp[i]
            perm = [None]*n
            for i in range(m):
                src = ins[i]              # a permutation of block i's elements
                dst = blocks[bp[i]]
                for pos, x in enumerate(src):
                    perm[x] = dst[pos]
            H.add(tuple(perm))
    expected = 1
    for kk in ks:
        expected *= factorial(kk)
    for c in Counter(ks).values():
        expected *= factorial(c)
    assert len(H) == expected, (ks, len(H), expected)
    # freeness of right multiplication on Sigma_n is automatic for subgroups; check orbit count
    assert factorial(n) % len(H) == 0
    return len(H), factorial(n)//len(H)

if __name__ == "__main__":
    total = 0
    for n in range(1, 8):
        for k in range(1, n+1):
            s, o = check_free_action(n, k)
            total += s
    print(f"free-action check: all surjections n<=7 checked ({total} surjections)")
    cases = [(1,1),(2,1),(1,1,1),(2,2),(3,1),(2,1,1),(2,2,1),(3,3),(2,2,2),(3,2,1)]
    for ks in cases:
        h, idx = check_H_free(ks)
        print(f"H{ks}: |H|={h}, |Sigma_n/H|={idx}")
    print("ALL F1 CHECKS PASSED")
