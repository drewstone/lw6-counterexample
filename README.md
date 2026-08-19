# Linden–Winter Conjecture 6 is false

Linden and Winter ([*Commun. Math. Phys.* 259:129–138 (2005)](https://arxiv.org/abs/quant-ph/0406162)) proved:
if a four-party quantum state satisfies `I(A;C|B) = I(C;B|A) = I(A;B|D) = 0`, then `I(C;D) ≥ I(C;AB)`.
Their **Conjecture 6** asked for the stable version: positive constants `k1, k2, k3` with

```
k1·I(A;C|B) + k2·I(C;B|A) + k3·I(A;B|D) + I(C;D) − I(C;AB) ≥ 0
```

for **all** quadripartite quantum states. Open since 2004.

**This repository refutes it.** A one-parameter family of classical distributions
(six atoms on four bits, integer weights in closed form) has `I(A;C|B) = I(A;B|D) = 0`
*exactly* — two combinatorial identities you can check by hand — while the ratio
`[I(C;AB) − I(C;D)] / I(C;B|A)` grows like `2·10^L`. No constants work.
Classical distributions are diagonal quantum states, so the quantum conjecture falls with them.

## The family

Atoms `(A,B,C,D)` and weights, for integer `L ≥ 2`:

| atom | weight |
|------|--------|
| 0000 | `2·10^(4L+3)` |
| 0010 | `27·10^(4L+3) − 27·10^(3L+3)` |
| 1101 | `2·10^(4L)` |
| 1111 | `27·10^(4L)` |
| 0101 | `2·10^(L+3)` |
| 0110 | `27·10^(L+3)` |

Why the two zeros are exact: every `D=0` atom has `A=0` and every `D=1` atom has `B=1`,
so `I(A;B|D) = 0`; and in the `B=1` block, `w(1101)·w(0110) = w(1111)·w(0101)`
(both are `54·10^(5L+3)`), so `I(A;C|B) = 0`.

## Verify it yourself

```
python3 verify.py                 # all certificates: identities, ladder, kills to k2 = 1e24
python3 verify.py 3.5e9 1e12 7    # find a certified refuting state for YOUR (k1,k2,k3)
```

Python 3.9+, standard library only. Probabilities are exact rationals; every logarithm is
bracketed with directed rounding (`ROUND_FLOOR`/`ROUND_CEILING`), so every certified sign
is a theorem about the printed digits, not a floating-point estimate. Runs in a few
minutes; exits non-zero if any certificate fails. CI runs it on every push.

## Paper

[`paper/main.pdf`](paper/main.pdf) — statement, the family, hand proofs of the two
identities, the asymptotic mechanism (`y = Θ(δ⁵)` vs `−e = Θ(δ⁴)`, so the ratio
diverges), the certificate table, and what survives (the constrained theorem is
untouched; the named route to an unconstrained inequality is closed).

## Provenance

The family was found by an automated research system searching supports and weight
patterns; the certificates were produced and cross-checked by four independently written
interval-arithmetic implementations with no shared code, of which `verify.py` is the
fourth. The identities and the asymptotics are human-checkable by hand.
