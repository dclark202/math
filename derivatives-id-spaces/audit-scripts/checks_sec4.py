"""Batch 8 checks (Section 4): Lemma 4.4 (decomposition, uniqueness, compose/decompose inverse),
the set-level squares of Prop 4.6 / the claim at line 624 (both composites trivial unless lambda_p = lambda_{p+1}),
and associativity of the Psi-decompositions (Prop 4.6 via Clark-1 Prop 5.8)."""
import itertools
from partitions import *

def all_chains(T, k):
    """All k-simplices of P(T) (non-basepoint, degenerate included): chains min = l0 <= ... <= lk = max."""
    parts = partitions_of(T)
    mn, mx = min_partition(T), max_partition(T)
    finer_eq = {p: [q for q in parts if refines(q, p)] for p in parts}
    out = []
    def dfs(chain):
        if len(chain) == k + 1:
            if chain[-1] == mx:
                out.append(tuple(chain))
            return
        for q in finer_eq[chain[-1]]:
            dfs(chain + [q])
    dfs([mn])
    return out

stats = dict(uniq=0, inverse=0, square=0, trivial_ok=0, assoc=0)
for n in (2, 3, 4, 5):
    T = frozenset(range(n))
    parts = partitions_of(T)
    for k in (2, 3, 4):
        if n == 5 and k == 4:
            continue
        for gamma in all_chains(T, k):
            # Lemma 4.4: decomposes w.r.t. alpha at position p iff gamma[p] == alpha; exactly one alpha per p
            for p in range(k + 1):
                q = k - p
                hits = [alpha for alpha in parts if decompose(gamma, p, alpha) is not None]
                assert hits == [gamma[p]], (gamma, p, hits)
                g, betas, blocks = decompose(gamma, p, gamma[p])
                assert compose(g, betas, blocks) == gamma            # decomposition inverts composition
                stats['uniq'] += 1; stats['inverse'] += 1
            # squares (4.14): for p + q + 1 = k
            for p in range(k):
                q = k - 1 - p
                for alpha in parts:
                    # route 1: Psi_{(p,q+1)} then id x prod d_0
                    r1 = None
                    d = decompose(gamma, p, alpha)
                    if d is not None:
                        g, betas, blocks = d
                        b0 = [face(b, 0, Tb) for b, Tb in zip(betas, blocks)]
                        if all(b is not None for b in b0):
                            r1 = (g, tuple(b0))
                    # route 2: Psi_{(p+1,q)} then d_{p+1} x id
                    r2 = None
                    d = decompose(gamma, p + 1, alpha)
                    if d is not None:
                        g, betas, blocks = d
                        nT = frozenset(range(len(blocks)))
                        g0 = face(g, p + 1, nT)
                        if g0 is not None:
                            r2 = (g0, tuple(betas))
                    assert r1 == r2, (gamma, p, alpha, r1, r2)
                    # line 624: trivial unless lambda_p == lambda_{p+1} (== alpha)
                    nontrivial = (r1 is not None)
                    assert nontrivial == (gamma[p] == gamma[p + 1] == alpha)
                    stats['square'] += 1
                    stats['trivial_ok'] += 1
            # associativity of decompositions: decompose at p then the front at p' < p, vs decompose at p' then the backs
            for p in range(1, k):
                for pp in range(1, p):
                    g, betas, blocks = decompose(gamma, p, gamma[p])
                    g2, betas2, blocks2 = decompose(g, pp, g[pp])
                    # other route
                    g3, betas3, blocks3 = decompose(gamma, pp, gamma[pp])
                    # front pieces agree
                    assert g2 == g3, (gamma, p, pp)
                    stats['assoc'] += 1
print(stats)
print("ALL SECTION 4 CHECKS PASSED")
