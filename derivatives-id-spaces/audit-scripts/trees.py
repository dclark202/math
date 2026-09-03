"""Trees as laminar families (Def A.6), weightings by heights (Def A.10), reduction (Lemma A.16),
alive partitions (Def A.19), the levelling map Theta (Constr. A.22) and its inverse (proof of Thm A.24)."""
import random
from partitions import max_partition

def is_tree(F, T):
    T = frozenset(T)
    if T not in F or any(frozenset([x]) not in F for x in T): return False
    if any((not B) or not B <= T for B in F): return False
    return all(A <= B or B <= A or not (A & B) for A in F for B in F)

def internal(F): return [B for B in F if len(B) > 1]
def root(F): return max(F, key=len)
def parent(F, B):
    c = [C for C in F if B < C]
    return min(c, key=len) if c else None
def children(F, B):
    return [C for C in F if C < B and not any(C < D < B for D in F)]
def leaf_parent(F, x): return parent(F, frozenset([x]))

def is_weighting(F, h, eps=1e-12):
    for B in internal(F):
        if not (-eps <= h[B] <= 1 + eps): return False
        p = parent(F, B)
        if p is not None and h[p] > h[B] + eps: return False
    return True

def is_reduced(F, h):
    return all(h[B] > h[parent(F, B)] for B in internal(F) if parent(F, B) is not None)

def is_basepoint(F, h):
    T = root(F)
    return h[T] == 0 or any(h[leaf_parent(F, x)] == 1 for x in T)

def reduce(F, h):
    """Lemma A.16: keep root, leaves, and internal non-root B with h(B) > h(p_F(B))."""
    T = root(F)
    keep = frozenset(B for B in F if len(B) == 1 or B == T or h[B] > h[parent(F, B)])
    return keep, {B: h[B] for B in keep if len(B) > 1}

def alive(F, h, t):
    """Def A.19 with h(p(T)) := 0 and h({x}) := 1."""
    def H(B): return 1.0 if len(B) == 1 else h[B]
    def Hp(B):
        p = parent(F, B); return 0.0 if p is None else h[p]
    return frozenset(B for B in F if H(B) >= t and Hp(B) < t)

def tree_of_chain(chain):
    return frozenset(B for lam in chain for B in lam)

def levels(chain):
    """a(B), b(B): first/last index j with B a block of chain[j]."""
    a, b = {}, {}
    for j, lam in enumerate(chain):
        for B in lam:
            a.setdefault(B, j); b[B] = j
    return a, b

def theta(chain, u):
    """Construction A.22: h(B) = u_{b(B)+1}, u = (u_1..u_k), u_0 = 0, u_{k+1} = 1."""
    F = tree_of_chain(chain)
    a, b = levels(chain)
    uu = [0.0] + list(u) + [1.0]
    return F, {B: uu[b[B] + 1] for B in F if len(B) > 1}

def theta_inverse(F, h):
    """Proof of Thm A.24 (III): u = distinct heights, lambda_j = alive at w_{j+1}, lambda_k = max."""
    ws = sorted(set(h[B] for B in internal(F)))
    chain = tuple(alive(F, h, w) for w in ws) + (max_partition(root(F)),)
    return chain, tuple(ws)

def all_trees(T):
    """All trees on T (laminar families with root T and all singletons)."""
    from partitions import partitions_of
    T = frozenset(T)
    if len(T) == 1:
        return [frozenset([T])]
    out = []
    for lam in partitions_of(T):
        if len(lam) < 2: continue
        # choose a tree on each block; cartesian product
        choices = [all_trees(B) for B in sorted(lam, key=sorted)]
        def prod(i, acc):
            if i == len(choices):
                out.append(frozenset([T]).union(*acc)); return
            for t in choices[i]:
                prod(i + 1, acc + [t])
        prod(0, [])
    return out

def random_weighting(F, strict=True):
    """Random order-preserving weighting, built root-down: h(B) in (h(p(B)), 1)."""
    T = root(F)
    h = {}
    def go(B, lo):
        if len(B) == 1: return
        h[B] = lo + (1 - lo) * random.random()
        for C in children(F, B): go(C, h[B])
    go(T, 0.0)
    return h
