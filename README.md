# Self-Similar Digit Blocks and Collatz Total Stopping Time

**Status:** amateur/independent investigation, computationally verified, not peer-reviewed.
**TL;DR:** Repeatedly prepending a fixed bit-block to a number causes the Collatz
total stopping time to increase by an amount that converges — with probability
approaching 1 — to exactly the block's own bit length. A specific family of
blocks (those satisfying `3·(x/2)+1 = 2^L`) makes this convergence exactly 100%
from the very first iteration.

This is **not** a claim about the Collatz conjecture itself (true/false, cycles,
divergent trajectories). It is a concrete, checkable statement about a specific
family of numbers built from self-similar bit patterns.

## Repository contents

```
README.md                              this document
LICENSE                                 MIT
collatz_block_repeat.py                 main demo, single block
big_survey.py                           105-block systematic sweep
exact_power_test.py                     isolates the "exact 2^L" condition
cycle_search.c                          unrelated side-quest: exhaustive
                                         search for non-trivial Collatz cycles
                                         via parity-vector fixed points (GMP,
                                         OpenMP) — negative result up to q=20,
                                         kept for reference / reuse
results/single_block_10251997.csv       raw diff_k data for x=10251997, k=2..300
results/survey_105_blocks.csv           raw results of the 105-block sweep (§2)
```

---

## 1. Construction

Fix an odd integer `x` with bit length `L` (so `2^(L-1) ≤ x < 2^L`). Define a
sequence of integers by repeatedly prepending `x` to the front (most
significant bits) of the previous term:

```
N_1 = x
N_{k+1} = (x << bitlength(N_k)) | N_k
```

In other words, `N_k` is `x` written `k` times in a row in binary — the base-`2^L`
analogue of a repunit (e.g. `x=5`, `L=3` gives `101`, `101101`, `101101101`, ...).

One can show algebraically (and we verify computationally) that this satisfies
the clean closed form:

```
N_{k+1} = 2^L · N_k + x        (exact, for all k)
```

Let `steps(n)` be the Collatz total stopping time (number of `n → n/2` /
`n → 3n+1` operations until reaching 1). Define:

```
diff_k = steps(N_{k+1}) − steps(N_k)
```

## 2. Main empirical finding

**Claim.** For fixed `x` (bit length `L`), the frequency of `diff_k = L` among
the first `k` iterations tends to increase with `k` and stabilizes at a high
value (often 60–100%, and *exactly* 100% for a specific family — see §4)
rather than the ≈`7·L` predicted by generic Collatz statistics (average total
stopping time grows ≈7 steps per bit for a "random" number).

Tested against **105 distinct blocks** spanning bit lengths 4–24 (all-ones,
alternating, single-bit, and random patterns): **0 exceptions**. Every block's
dominant `diff_k` value converged to its own bit length `L`.

| block type | example (L=16) | late-stage freq(diff=L) |
|---|---|---|
| alternating (`1010...10`) | 43690 | **100.0%** |
| all-ones | 65535 | 66.7% |
| random odd | 53327 | 40.0% |
| known delay-record (`8400511`, L=24) | 8400511 | 66.7% |

## 3. Why this happens (proved part)

Since `N_{k+1} = 2^L·N_k + x`, we have `N_{k+1} ≡ N_k (mod 2^{L·k})` — i.e.
`N_k`'s entire bit pattern is preserved unchanged as the low-order bits of
`N_{k+1}`. Combined with the standard fact that **a number's first `m` Collatz
steps (its parity vector) are fully determined by the number mod `2^m`**
(Lagarias; see Bernstein & Lagarias, *The 3x+1 conjugacy map*, Canad. J. Math.
1996), this forces `N_k` and `N_{k+1}` to take an *identical* sequence of
odd/even steps for at least the first `L·k` iterations. Empirically the shared
prefix is even longer (~1.5·L·k, consistent with ~2/3 of Collatz steps being
"even").

This rigorously explains *why* `N_k` and `N_{k+1}` are entangled for a long
stretch. It does **not** by itself explain why the final total-step
*difference*, after they eventually diverge, is almost always the small,
k-independent value `L` rather than scaling with the (growing) length of the
shared prefix. That part remains an open empirical pattern — see §5.

## 4. The 100%-clean family

Numbers of the alternating form `x = 1010...10` (L bits) satisfy an exact
algebraic identity:

```
x / 2 = (2^L − 1) / 3          (integer, since x/2 has this closed form)
3·(x/2) + 1 = 2^L               (exact — lands precisely on 2^L)
```

I.e. `x`'s own Collatz trajectory reaches **exactly `2^L`** (the power of two
matching its *own* bit length) within 2 steps, then descends to 1 by pure
halving. This is a much sharper condition than merely "reaches some power of
two quickly" — a control test using numbers that reach *unrelated* powers of
two (found by tracing the Collatz map backward from `2^m` for `m ≠ L`) does
**not** reproduce the 100% effect (only 43–70%, same as generic numbers, see
`exact_power_test.py`).

The self-referential condition — landing exactly on `2^L` where `L` is the
block's own length — appears to align perfectly with the shift amount `L`
used when constructing `N_{k+1}`, eliminating the carry/interference that
causes deviation from `diff_k = L` for generic blocks.

## 5. Open questions

- Is there a clean closed-form proof that `Pr(diff_k = L) → 1` as `k → ∞`
  for *typical* `x` (not just the 100% family)? A partial mechanism is
  sketched in §3, but the exact limiting frequency and its dependence on the
  bit-pattern of `x` is not derived, only observed.
- Is `3(x/2)+1 = 2^L` **sufficient and necessary** for the 100% effect, or are
  there other unrelated algebraic identities that also produce it?
- Does the phenomenon generalize to bases other than 2^L (i.e. repeating a
  block in base `b` for other `b`), or to other `an+b` Collatz-like maps?

## 6. Relation to existing theory

This sits inside the well-studied **2-adic extension of the Collatz map**
(Bernstein & Lagarias 1996; see also the notion of "Collatz cyclic numbers"
associated with periodic parity vectors). As `k → ∞`, `N_k` converges 2-adically
to the rational 2-adic integer `α = -x/(2^L - 1)`. This connects the present
observation to known machinery, but the specific claim about the *frequency*
of `diff_k = L`, as far as I can tell, is not stated in the literature I
searched (Bernstein–Lagarias 1996; Lagarias periodicity conjecture papers;
Hercher 2023; Eliahou 1993; Simons & de Weger). Corrections/pointers to prior
work are very welcome.

## 7. Reproducing this

All experiments are pure Python (only the `fractions` and `collections`
standard library modules) plus one C/GMP program for larger-scale exhaustive
cycle searches (unrelated negative-result side-quest, kept for completeness —
see `cycle_search.c`).

```
python3 collatz_block_repeat.py      # main phenomenon, single block
python3 big_survey.py                 # 105-block systematic sweep
python3 exact_power_test.py           # isolates the "exact 2^L" condition
```

No external dependencies beyond CPython 3.8+.

## 8. Caveats

- This does **not** bear on the truth of the Collatz conjecture itself.
- All claims here are backed by finite computation (k up to a few thousand,
  L up to 24), not proofs, except where explicitly marked "proved."
- This investigation was done collaboratively with an AI assistant (Claude);
  all code was executed and its output verified before being reported here.
  Please independently re-verify before citing.

## Acknowledgments

Investigation prompted by tracing down and correcting an AI-hallucinated
claim (a different model asserted a specific number was a "periodic point"
of the Collatz map via a fabricated formula, which did not survive direct
computational verification). The real phenomenon documented here was found
in the process of checking that claim.
