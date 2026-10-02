#!/usr/bin/env python3
"""Verifier for the refutation of Linden-Winter Conjecture 6.

Conjecture 6 (Linden-Winter, quant-ph/0406162, Commun. Math. Phys. 259:129-138, 2005):
there exist positive constants k1, k2, k3 such that for every quadripartite quantum
state,

    f_k = k1*I(A;C|B) + k2*I(C;B|A) + k3*I(A;B|D) + I(C;D) - I(C;AB)  >=  0.

This script certifies that NO positive constants exist, using the one-parameter family
of classical distributions W_L (classical distributions are diagonal quantum states, so
a classical counterexample refutes the quantum conjecture):

    atoms (A,B,C,D):   (0,0,0,0) (0,0,1,0) (1,1,0,1) (1,1,1,1) (0,1,0,1) (0,1,1,0)
    weights W_L:        2*10^(4L+3),  27*10^(4L+3) - 27*10^(3L+3),
                        2*10^(4L),    27*10^(4L),
                        2*10^(L+3),   27*10^(L+3)

Everything is exact or one-sidedly rounded:
  * probabilities are exact integer ratios (fractions.Fraction);
  * the two structural identities I(A;C|B) = 0 and I(A;B|D) = 0 are verified as exact
    integer statements about the support and weights -- no floating point at all;
  * every entropy is bracketed [lo, hi] with decimal directed rounding (ROUND_FLOOR /
    ROUND_CEILING contexts around a correctly-rounded Decimal.ln), so every certified
    sign and every certified lower bound is a THEOREM about the printed digits, not a
    floating-point estimate.

What is certified on each run:
  1. x = I(A;C|B) = 0 and z = I(A;B|D) = 0 exactly, for every L checked      [integer]
  2. y = I(C;B|A) > 0 and e = I(C;D) - I(C;AB) < 0, for every L checked      [interval]
  3. R_L = -e/y lower bounds for L = 2..16, and R_{L+1}/R_L >= 1.9           [interval]
  4. explicit kills f_k(W_L) < 0 for a menu of k reaching k = (10^24)        [interval]
  5. any k you pass on the command line:  verify.py 3.5e9 1e12 7  finds an L
     with a certified kill for k = (3.5e9, 1e12, 7)                          [interval]

Because x = z = 0 exactly, f_k = k2*y + e on the family, so a kill certificate for
(0, k2, 0) covers every (k1, k2, k3) with the same k2: the k1 and k3 terms multiply
exact zeros. For k2 = 0 and any other k_i > 0, f_k = e < 0 -- one row settles it.

Requires Python 3.9+. Stdlib only. Exit code 0 iff every certificate holds.
"""

import sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext, localcontext
from fractions import Fraction

# Atoms in (A, B, C, D) order.
SUPPORT = [(0, 0, 0, 0), (0, 0, 1, 0), (1, 1, 0, 1), (1, 1, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0)]
A, B, C, D = 0, 1, 2, 3

FAILURES = []


def check(label, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {label}")
    if not ok:
        FAILURES.append(label)
    return ok


def weights(L):
    return [
        2 * 10 ** (4 * L + 3),
        27 * 10 ** (4 * L + 3) - 27 * 10 ** (3 * L + 3),
        2 * 10 ** (4 * L),
        27 * 10 ** (4 * L),
        2 * 10 ** (L + 3),
        27 * 10 ** (L + 3),
    ]


# ----------------------------------------------------------------------------- exact part

def structural_zeros(w):
    """x = z = 0 as exact integer statements. Returns (x_is_zero, z_is_zero).

    z = I(A;B|D): on every D=0 atom A is constant, and on every D=1 atom B is
    constant, so conditioned on D one side of the mutual information is
    deterministic and z vanishes identically -- a support property alone.

    x = I(A;C|B): on every B=0 atom A is constant (x's B=0 block vanishes), and on
    the B=1 block {1101, 1111, 0101, 0110} independence of A and C is exactly the
    2x2 determinant identity w(1101)*w(0110) == w(1111)*w(0101).
    """
    live = [(s, wi) for s, wi in zip(SUPPORT, w) if wi != 0]
    z_zero = len({s[A] for s, _ in live if s[D] == 0}) <= 1 and \
             len({s[B] for s, _ in live if s[D] == 1}) <= 1
    x_b0 = len({s[A] for s, _ in live if s[B] == 0}) <= 1
    by = {s: wi for s, wi in live}
    x_b1 = by.get((1, 1, 0, 1), 0) * by.get((0, 1, 1, 0), 0) == \
           by.get((1, 1, 1, 1), 0) * by.get((0, 1, 0, 1), 0)
    return x_b0 and x_b1, z_zero


# -------------------------------------------------------------------------- interval part

def _ln_bracket(q, prec):
    """[lo, hi] on ln(q) for an exact Fraction q > 0, widened one ulp each way.

    Every arithmetic op stays inside the directed context. Decimal applies the AMBIENT
    context to any operation performed outside a localcontext, silently truncating to
    its precision -- performing the ulp widening outside destroyed these brackets from
    L = 11 on while leaving L <= 10 correct, which is exactly the wrong kind of bug: the
    small cases that get eyeballed pass, and the deep rows that carry the theorem fail.
    """
    with localcontext() as ctx:
        ctx.prec = prec
        ctx.rounding = ROUND_FLOOR
        lo = (Decimal(q.numerator) / Decimal(q.denominator)).ln()
        ulp = Decimal(10) ** (lo.adjusted() - prec + 2)
        lo = lo - ulp
        ctx.rounding = ROUND_CEILING
        hi = (Decimal(q.numerator) / Decimal(q.denominator)).ln()
        hi = hi + ulp
    return lo, hi


def entropy_bracket(masses, prec):
    """[lo, hi] on the Shannon entropy (nats) of exact Fraction masses."""
    lo = hi = Decimal(0)
    for q in masses:
        if q <= 0:
            continue
        l_lo, l_hi = _ln_bracket(q, prec)
        with localcontext() as ctx:
            ctx.prec = prec
            ctx.rounding = ROUND_FLOOR
            qd = Decimal(q.numerator) / Decimal(q.denominator)
            lo += qd * (-l_hi)
            ctx.rounding = ROUND_CEILING
            qd = Decimal(q.numerator) / Decimal(q.denominator)
            hi += qd * (-l_lo)
    return lo, hi


def coords(L, prec):
    """Interval brackets for y = I(C;B|A) and e = I(C;D) - I(C;AB) at W_L (nats)."""
    w = weights(L)
    total = sum(w)
    probs = [Fraction(wi, total) for wi in w]

    def marginal(idx):
        m = {}
        for atom, p in zip(SUPPORT, probs):
            key = tuple(atom[i] for i in idx)
            m[key] = m.get(key, Fraction(0)) + p
        return [v for v in m.values() if v > 0]

    S = {tuple(i): entropy_bracket(marginal(list(i)), prec) for i in
         [(A,), (D,), (A, B), (A, C), (C, D), (A, B, C)]}

    def add(u, v):
        with localcontext() as ctx:
            ctx.prec = prec
            ctx.rounding = ROUND_FLOOR
            lo = u[0] + v[0]
            ctx.rounding = ROUND_CEILING
            hi = u[1] + v[1]
        return lo, hi

    def sub(u, v):
        with localcontext() as ctx:
            ctx.prec = prec
            ctx.rounding = ROUND_FLOOR
            lo = u[0] - v[1]
            ctx.rounding = ROUND_CEILING
            hi = u[1] - v[0]
        return lo, hi

    y = sub(add(S[(A, B)], S[(A, C)]), add(S[(A,)], S[(A, B, C)]))
    e = sub(add(S[(A, B, C)], S[(D,)]), add(S[(A, B)], S[(C, D)]))
    return y, e


def ratio_bracket(L, prec):
    """[lo, hi] on R_L = -e/y (dimensionless: the ln-2 factors cancel)."""
    y, e = coords(L, prec)
    if not (y[0] > 0 and e[1] < 0):
        return None
    with localcontext() as ctx:
        ctx.prec = prec
        ctx.rounding = ROUND_FLOOR
        lo = (-e[1]) / y[1]
        ctx.rounding = ROUND_CEILING
        hi = (-e[0]) / y[0]
    return lo, hi


def prec_for(L):
    return 6 * L + 160


def kill(k1, k2, k3, L, prec):
    """Certified f_k(W_L) < 0? Returns the certified upper bound, or None.

    x = z = 0 exactly (verified separately), so f_k = k2*y + e.
    """
    k2 = Fraction(k2)
    if k2 < 0:
        raise ValueError("The certificate requires a nonnegative k2")
    y, e = coords(L, prec)
    with localcontext() as ctx:
        ctx.prec = prec
        ctx.rounding = ROUND_CEILING
        # Round the exact rational upward before multiplying the positive y bound.
        k2_upper = Decimal(k2.numerator) / Decimal(k2.denominator)
        hi = k2_upper * y[1] + e[1]
    return hi if hi < 0 else None


def find_kill(k1, k2, k3, l_max=40):
    for L in range(2, l_max + 1):
        w = weights(L)
        xz = structural_zeros(w)
        if not (xz[0] and xz[1]):
            return None
        hi = kill(k1, k2, k3, L, prec_for(L))
        if hi is not None:
            return L, hi
    return None


# --------------------------------------------------------------------------------- driver

def main():
    print(__doc__.splitlines()[0])
    print()

    # 1. structural zeros, exact, for every L used anywhere below
    ok = all(all(structural_zeros(weights(L))) for L in range(2, 25))
    check("x = I(A;C|B) = 0 and z = I(A;B|D) = 0 exactly for L = 2..24 [integer identity]", ok)

    # 2 + 3. signs and the ladder
    print()
    print(f"  {'L':>3} {'R_L certified lower bound':>30} {'ratio to previous (lower bound)':>34}")
    prev_hi = None
    ratios_ok = True
    signs_ok = True
    for L in range(2, 17):
        prec = prec_for(L)
        y, e = coords(L, prec)
        signs_ok &= y[0] > 0 and e[1] < 0
        r = ratio_bracket(L, prec)
        if r is None:
            print(f"  {L:>3} {'NOT CERTIFIED — interval spans zero':>30}")
            signs_ok = ratios_ok = False
            prev_hi = None
            continue
        with localcontext() as ctx:
            ctx.prec = prec
            ctx.rounding = ROUND_FLOOR
            step = None if prev_hi is None else r[0] / prev_hi
        print(f"  {L:>3} {str(r[0])[:28]:>30} {(str(step)[:14] if step is not None else '-'):>34}")
        if step is not None:
            ratios_ok &= step >= Decimal("1.9")
        prev_hi = r[1]
    check("y > 0 and e < 0 certified for L = 2..16", signs_ok)
    check("R_{L+1}/R_L >= 1.9 certified for L = 2..15 (measured ~10)", ratios_ok)

    # 4. explicit kills. Each certifies f_k < 0 for EVERY k sharing that k2,
    #    because the k1 and k3 terms multiply exact zeros.
    print()
    menu = [
        ((1, 0, 0), 2), ((0, 0, 1), 2), ((1, 1, 1), 2),
        ((10**9, 10**9, 10**9), 15), ((10**15, 10**15, 10**15), 18),
        ((10**30, 100, 10**30), 5), ((0, 10**20, 0), 22), ((0, 10**24, 0), 26),
    ]
    for k, L in menu:
        hi = kill(*k, L, prec_for(L))
        check(f"f_k(W_{L}) < 0 certified for k = {k}   [upper bound {hi}]", hi is not None)

    # 5. caller-supplied k
    if len(sys.argv) == 4:
        k = tuple(Fraction(a) for a in sys.argv[1:4])
        print()
        got = find_kill(*k)
        if got:
            L, hi = got
            check(f"f_k(W_{L}) < 0 certified for caller k = {tuple(map(str, k))} [bound {hi}]", True)
        else:
            check(f"no kill found for caller k = {tuple(map(str, k))} within L <= 40", False)

    print()
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} certificate(s) FAILED")
        return 1
    print("RESULT: all certificates hold — no positive (k1, k2, k3) satisfies Conjecture 6")
    return 0


if __name__ == "__main__":
    getcontext().prec = 50
    raise SystemExit(main())
