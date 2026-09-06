"""Check for F2 (memo INFTY-COMPARISON.md, A1): is the cooperad structure on the
constant symmetric sequence S_ (S in every arity) determined at the level of
pi_0?  A coalgebra structure on S_ for the composition product has, in arity n,
a comultiplication S -> (+)_{pi in Part(n)} S, i.e. an integer a(pi) per
partition, Sigma_n-invariant (so a depends only on the multiset of block sizes,
the *type* of pi).  Counit: a(one block) = a(discrete) = 1.  Coassociativity:
for every chain rho <= pi (rho finer),
    a(pi) * prod_{B in pi} a(rho|_B)  =  a(rho) * a(pi/rho).
We solve these multiplicative equations in log-space over Q (a = exp(x), x a
Q-vector), compute the dimension of the solution space, and compare with the
space of 'coboundaries' x(pi) = h(k(pi)) + sum_B h(|B|) - h(n) (h(1)=0), which
are the structures obtained from the standard one by the endomorphism of S_
that is 'multiplication by h(n)' in arity n.  Ends with a summary line.
"""
from fractions import Fraction
from itertools import combinations
from functools import lru_cache

def set_partitions(elems):
    elems = list(elems)
    if not elems:
        yield []
        return
    first = elems[0]
    for smaller in set_partitions(elems[1:]):
        for i, block in enumerate(smaller):
            yield smaller[:i] + [[first] + block] + smaller[i+1:]
        yield [[first]] + smaller

def ptype(blocks):
    return tuple(sorted(len(b) for b in blocks))

def finer(rho, pi):
    # rho finer than pi: every block of rho inside a block of pi
    return all(any(set(b) <= set(B) for B in pi) for b in rho)

def induced_type_on_block(rho, B):
    return tuple(sorted(len(b) for b in rho if set(b) <= set(B)))

def quotient_type(rho, pi):
    return tuple(sorted(sum(1 for b in rho if set(b) <= set(B)) for B in pi))

def rank_and_nullity(rows, nvars):
    # Gaussian elimination over Q; returns rank
    m = [list(r) for r in rows]
    rank = 0
    col = 0
    nrows = len(m)
    for col in range(nvars):
        piv = None
        for r in range(rank, nrows):
            if m[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        pv = m[rank][col]
        m[rank] = [v / pv for v in m[rank]]
        for r in range(nrows):
            if r != rank and m[r][col] != 0:
                f = m[r][col]
                m[r] = [a - f*b for a, b in zip(m[r], m[rank])]
        rank += 1
    return rank

def analyse(N):
    # variables: one per (n, type) with 2 <= n <= N, excluding the two trivial types
    types = {}
    for n in range(1, N+1):
        for p in set_partitions(range(n)):
            t = (n, ptype(p))
            if t not in types:
                types[t] = None
    var_index = {}
    for t in sorted(types):
        n, tp = t
        if len(tp) == 1 or len(tp) == n:   # one block, or discrete: forced = 1
            continue
        var_index[t] = len(var_index)
    nv = len(var_index)
    def vec_for(n, tp):
        v = [Fraction(0)]*nv
        t = (n, tp)
        if t in var_index:
            v[var_index[t]] = Fraction(1)
        return v
    rows = []
    for n in range(2, N+1):
        parts = list(set_partitions(range(n)))
        for rho in parts:
            for pi in parts:
                if rho is pi or not finer(rho, pi):
                    continue
                lhs = vec_for(n, ptype(pi))
                for B in pi:
                    v = vec_for(len(B), induced_type_on_block(rho, B))
                    lhs = [a+b for a, b in zip(lhs, v)]
                rhs = vec_for(n, ptype(rho))
                v = vec_for(len(rho), quotient_type(rho, pi))
                rhs = [a+b for a, b in zip(rhs, v)]
                row = [a-b for a, b in zip(lhs, rhs)]
                if any(row):
                    rows.append(row)
    rank = rank_and_nullity(rows, nv) if rows else 0
    sol_dim = nv - rank
    # coboundaries: h(2..N) free, x(n,tp) = h(k) + sum h(|B|) - h(n)
    hvars = list(range(2, N+1))
    cob = []
    for hn in hvars:
        v = [Fraction(0)]*nv
        for t, idx in var_index.items():
            n, tp = t
            k = len(tp)
            val = (1 if k == hn else 0) + sum(1 for s in tp if s == hn) - (1 if n == hn else 0)
            v[idx] = Fraction(val)
        cob.append(v)
    cob_rank = rank_and_nullity(cob, nv) if cob else 0
    return nv, len(rows), sol_dim, cob_rank

if __name__ == "__main__":
    for N in range(3, 7):
        nv, neq, sol, cob = analyse(N)
        print(f"arities <= {N}: {nv} unknown coefficients, {neq} equations, "
              f"solution space dim = {sol}, coboundary space dim = {cob}")
    print("F2 pi_0 ANALYSIS DONE")
