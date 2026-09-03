"""Chains of partitions of a finite set (the pointed simplicial set P(n), Def 2.4),
profiles, faces/degeneracies, and the decomposition Psi of Lemma 4.4.
A partition is a frozenset of frozensets (blocks). Chains run coarse -> fine."""

def set_partitions(elems):
    elems = list(elems)
    if not elems:
        yield []
        return
    first = elems[0]
    for smaller in set_partitions(elems[1:]):
        for i in range(len(smaller)):
            yield smaller[:i] + [smaller[i] | {first}] + smaller[i+1:]
        yield [frozenset([first])] + smaller

def partitions_of(T):
    return [frozenset(p) for p in set_partitions(T)]

def refines(fine, coarse):
    return all(any(b <= c for c in coarse) for b in fine)

def leq(lam, mu):
    """Paper's order (line 283): lam <= mu iff mu refines lam."""
    return refines(mu, lam)

def min_partition(T): return frozenset([frozenset(T)])
def max_partition(T): return frozenset(frozenset([x]) for x in T)

def nondegenerate_chains(T):
    """Non-basepoint nondegenerate simplices of P(T): strict chains min < ... < max."""
    T = frozenset(T)
    parts = partitions_of(T)
    mn, mx = min_partition(T), max_partition(T)
    finer = {p: [q for q in parts if q != p and refines(q, p)] for p in parts}
    out = []
    def dfs(chain):
        if chain[-1] == mx:
            out.append(tuple(chain)); return
        for q in finer[chain[-1]]:
            dfs(chain + [q])
    dfs([mn])
    return out

def is_chain(chain):
    return all(leq(chain[i], chain[i+1]) for i in range(len(chain)-1))

def face(chain, i, T):
    """d_i removes entry i; None = basepoint (Def 2.4)."""
    c = chain[:i] + chain[i+1:]
    if not c or c[0] != min_partition(T) or c[-1] != max_partition(T):
        return None
    return c

def degeneracy(chain, j):
    """s_j repeats entry j."""
    return chain[:j+1] + chain[j:]

def profile(chain):
    """|alpha|: for each level j>=1 and each block of lambda_{j-1}, the number of blocks of lambda_j inside it."""
    out = [len(chain[0])]
    for j in range(1, len(chain)):
        out.append(tuple(sorted(len([b for b in chain[j] if b <= c]) for c in chain[j-1])))
    return tuple(out)

def collapse(chain_prefix, blocks):
    """Collapse block blocks[i] to the point i (chain_prefix consists of partitions coarser than {blocks})."""
    out = []
    for lam in chain_prefix:
        out.append(frozenset(frozenset(i for i, Tb in enumerate(blocks) if Tb <= C) for C in lam))
    return tuple(out)

def restrict(chain_suffix, Tb):
    """Restrict partitions finer than {..,Tb,..} to the block Tb."""
    return tuple(frozenset(C for C in lam if C <= Tb) for lam in chain_suffix)

def decompose(chain, p, alpha):
    """Psi_{k,(p,q)} at the factor alpha (a partition): (gamma, betas) if chain[p] == alpha else None (Lemma 4.4)."""
    if chain[p] != alpha:
        return None
    blocks = sorted(alpha, key=lambda b: sorted(b))
    return collapse(chain[:p+1], blocks), [restrict(chain[p:], Tb) for Tb in blocks], blocks

def compose(gamma, betas, blocks):
    """gamma o (beta_1..beta_n) of the proof of Lemma 4.4: lambda'_0<=..<=lambda'_p (= union of blocks) then unions of the betas."""
    q = len(betas[0]) - 1
    front = [frozenset(frozenset().union(*[blocks[i] for i in D]) for D in lam) for lam in gamma[:-1]]
    back = [frozenset().union(*[beta[r] for beta in betas]) for r in range(q + 1)]
    return tuple(front + back)
