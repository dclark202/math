"""Counterexamples to Prop A.49 / Prop A.58 (continuity of kappa for every framed / windowed system),
and sanity checks that the paper's own cells are continuous at the same points.
We work in arity 2 (B(2) = [0,1]/{0,1}) and arity 3 (triangle x = h(3) <= y = h({1,2}))."""
from cuts import *
from trees import *

T2 = frozenset([0, 1]); T3 = frozenset([0, 1, 2])
B12 = frozenset([0, 1])
F2 = frozenset([T2, frozenset([0]), frozenset([1])])
F3 = frozenset([T3, B12, frozenset([0]), frozenset([1]), frozenset([2])])
flag2 = (frozenset([T2]), frozenset([frozenset([0]), frozenset([1])]))
alpha3 = frozenset([B12, frozenset([2])])
flag3 = two_flag(F3, alpha3)
nu = (1, B12)

def describe(res):
    """Coarse description of a cut value: None or dict pos -> rescaled heights."""
    if res is None:
        return "BASEPOINT"
    return {pos: {tuple(sorted(tuple(sorted(C)) for C in D)): round(v, 4) for D, v in hn.items()} for pos, (V, hn) in res.items()}

print("=== Counterexample 1 (arity 2, flag (min,max), single position): S = h - g, E = h + g, g = min(|h-1/2|, h, 1-h)")
def frames1(F, h, pos):
    hv = h[T2]; g = min(abs(hv - 0.5), hv, 1 - hv)
    return (hv - g, hv + g)
# axioms: values in [0,1], S <= E, S >= h(p(root)) = 0  -> (F1)-(F3) hold; (W2)-(W4) hold for the windowed version
for hv in [0.3, 0.45, 0.499, 0.4999, 0.5, 0.5001, 0.51, 0.7]:
    print(f"  h = {hv:<7} ->", describe(kappa(F2, {T2: hv}, flag2, frames1)))

print("=== Counterexample 2 (arity 3, anchored floor S = x = h(p(B)), ceiling E = min(1, 2y - x)); wall y -> x")
def frames2(F, h, pos):
    x, y = h[T3], h[B12]
    if pos == (0, T3): return (0.0, 1.0)
    return (x, min(1.0, 2 * y - x))          # S = h(p(B)) (anchored), E - S = y - x -> 0 at the wall
for eps in [0.2, 0.05, 0.01, 0.001, 0.0]:
    x, y = 0.3, 0.3 + eps
    print(f"  (x,y) = ({x},{y}) ->", describe(kappa(F3, {T3: x, B12: y}, flag3, frames2)))
print("  at the wall y = x the point of B(3) is the corolla, where the flag {1,2} is absent -> value must be BASEPOINT;")
print("  the values along y -> x stay at (bottom: x, top: 1/2).  So kappa is discontinuous although (F1)-(F3) and (W1)-(W4) hold.")

print("=== Counterexample 3 (arity 3, anchored, ceiling 1: S = max(x, 2y - 1)); approach the basepoint y -> 1")
def frames3(F, h, pos):
    x, y = h[T3], h[B12]
    if pos == (0, T3): return (0.0, 1.0)
    return (max(x, 2 * y - 1), 1.0)
for y in [0.8, 0.95, 0.99, 0.999, 1.0]:
    print(f"  (x,y) = (0.3,{y}) ->", describe(kappa(F3, {T3: 0.3, B12: y}, flag3, frames3)))
print("  y = 1 is the basepoint of B(3) (zero leaf edges); the values do not tend to the basepoint.")

print("=== The paper's cells at the same limits (tree frames, uniform labelled frames a=1/2, straight line s=0.3)")
lab = labelled_frames(flag3, uniform_labelling(flag3, [0.5]))
def mix(F, h, pos, s=0.3):
    (S1, E1), (S2, E2) = tree_frames(F, h, pos), lab(F, h, pos)
    return ((1 - s) * S1 + s * S2, (1 - s) * E1 + s * E2)
for name, fr in [("tree", tree_frames), ("labelled a=1/2", lab), ("mix s=0.3", mix)]:
    print(" ", name)
    for (x, y) in [(0.3, 0.31), (0.3, 0.301), (0.3, 0.3), (0.3, 0.99), (0.3, 0.999), (0.3, 1.0), (0.999, 0.9995), (0.4, 0.999)]:
        print(f"    (x,y) = ({x},{y}) ->", describe(kappa(F3, {T3: x, B12: y}, flag3, fr)))
