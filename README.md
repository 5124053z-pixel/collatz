# Self-Similar Digit Blocks and Collatz Total Stopping Time

Status: amateur/independent investigation, computationally verified, not peer-reviewed. TL;DR: Repeatedly prepending a fixed bit-block to a number causes the Collatz total stopping time to increase by an amount that converges — with probability approaching 1 — to exactly the block's own bit length. A specific family of blocks (those satisfying 3·(x/2)+1 = 2^L) makes this convergence exactly 100% from the very first iteration.

This is not a claim about the Collatz conjecture itself (true/false, cycles, divergent trajectories). It is a concrete, checkable statement about a specific family of numbers built from self-similar bit patterns.

Update: this investigation led to a bigger and more surprising finding than the original construction — see §5d. The well-known claim that "about half" of consecutive integers n, n+1 share the same Collatz total stopping time (OEIS A006577) appears, based on sampling out to n ~ 10^27000, to be a small-n snapshot of a quantity climbing steadily toward 100%, not a stable ~50% limit.

Further update: an Addendum (below) reframes §5's open questions using coupling theory (coalescing Markov chains) and reports new computational findings, including corrections of two errors made mid-investigation. See the Addendum section at the end of this document.

## Repository contents

```
README.md                              this document
LICENSE                                 MIT
collatz_block_repeat.py                 main demo, single block
big_survey.py                           105-block systematic sweep
exact_power_test.py                     isolates the "exact 2^L" condition
random_m_test.py                        §5a: tests whether self-similarity of
                                         N_k is even necessary (it isn't)
merging_residue_classes.py              §5b: classifies which residues mod 2^k
                                         provably force n, n+1 to merge
merging_classes_fast.c                  §5b: fast C/OpenMP version for larger k
general_merging_test.py                 §5c: extends §5b's merging-class
                                         analysis to general (L, x)
large_scale_sampling.c                  §5d: samples agreement rate at huge
                                         bit lengths (GMP + OpenMP)
cycle_search.c                          unrelated side-quest: exhaustive
                                         search for non-trivial Collatz cycles
                                         via parity-vector fixed points (GMP,
                                         OpenMP) — negative result up to q=20,
                                         kept for reference / reuse
results/single_block_10251997.csv       raw diff_k data for x=10251997, k=2..300
results/survey_105_blocks.csv           raw results of the 105-block sweep (§2)

Addendum tooling (coupling-theoretic follow-up, see Addendum section):
coupling_experiment3.py                 exhaustive early-merge check, N up to 1M
coupling_experiment4.py                 exhaustive early-merge check, N up to 5M
coupling_experiment5.py                 memoized exhaustive check, N up to 6x10^7
coupling_scaling.c                      C + GMP + OpenMP, windowed random
                                         sampling across bit-length scales
                                         (32 to 65536+ bits)
fit_gamma.py                            fits coupling-time decay rate gamma
                                         per bit-length window from
                                         coupling_scaling.c output
fast_slow_test.py                       tests whether merges are reproduced
                                         by fixing n's low K bits (sharp
                                         transition control test)
minimal_K_test.py                       binary search for smallest modulus
                                         2^K explaining each fast merger
hypothesis2_test.py                     correlates merge time against n's own
                                         trajectory length and 2-adic/run-
                                         length features

Follow-up (§10-12, after the Addendum tooling above was found to be lost):
k_cluster_analysis.py                   independent re-derivation of §6's
                                         K-distribution; finds a false-
                                         positive risk in minimal-K search
                                         and resolves the §4-connection
                                         question (negatively)
quantile_scaling_analysis.py            independent re-derivation of §5's
                                         quantile power-law scaling
                                         (narrow bit range)
quantile_scaling_large.py               same, at the original's full
                                         64-65536 bit range; matches the
                                         original exponents closely
coupling_exhaustive_verify.py           independent re-derivation of §2's
                                         exhaustive early-merge check;
                                         reproduces the original table
                                         exactly and extends past it
base_generalization_test.py             §13: control test, generic (non
                                         2/3-power) multipliers show no
                                         effect
power_of_3_test.py                      §13: main finding, c=3^p gives
                                         diff=-p, generalizing the
                                         phenomenon beyond base 2
power_of_3_merging_classes.py           §13: §5b-style provable merging
                                         classes for the (m, 3^p*m+x)
                                         pairing; every class found
                                         gives exactly diff=-p
```

## 1. Construction

Fix an odd integer x with bit length L (so 2^(L-1) ≤ x < 2^L). Define a sequence of integers by repeatedly prepending x to the front (most significant bits) of the previous term:

```
N_1 = x
N_{k+1} = (x << bitlength(N_k)) | N_k
```

In other words, N_k is x written k times in a row in binary — the base-2^L analogue of a repunit (e.g. x=5, L=3 gives 101, 101101, 101101101, ...).

One can show algebraically (and we verify computationally) that this satisfies the clean closed form:

```
N_{k+1} = 2^L · N_k + x        (exact, for all k)
```

Let steps(n) be the Collatz total stopping time (number of n → n/2 / n → 3n+1 operations until reaching 1). Define:

```
diff_k = steps(N_{k+1}) − steps(N_k)
```

## 2. Main empirical finding

Claim. For fixed x (bit length L), the frequency of diff_k = L among the first k iterations tends to increase with k and stabilizes at a high value (often 60–100%, and exactly 100% for a specific family — see §4) rather than the ≈7·L predicted by generic Collatz statistics (average total stopping time grows ≈7 steps per bit for a "random" number).

Tested against 105 distinct blocks spanning bit lengths 4–24 (all-ones, alternating, single-bit, and random patterns): 0 exceptions. Every block's dominant diff_k value converged to its own bit length L.

| block type | example (L=16) | late-stage freq(diff=L) |
|---|---|---|
| alternating (1010...10) | 43690 | 100.0% |
| all-ones | 65535 | 66.7% |
| random odd | 53327 | 40.0% |
| known delay-record (8400511, L=24) | 8400511 | 66.7% |

## 3. Why this happens (proved part)

Since N_{k+1} = 2^L·N_k + x, we have N_{k+1} ≡ N_k (mod 2^{L·k}) — i.e. N_k's entire bit pattern is preserved unchanged as the low-order bits of N_{k+1}. Combined with the standard fact that a number's first m Collatz steps (its parity vector) are fully determined by the number mod 2^m (Lagarias; see Bernstein & Lagarias, The 3x+1 conjugacy map, Canad. J. Math. 1996), this forces N_k and N_{k+1} to take an identical sequence of odd/even steps for at least the first L·k iterations. Empirically the shared prefix is even longer (~1.5·L·k, consistent with ~2/3 of Collatz steps being "even").

This rigorously explains why N_k and N_{k+1} are entangled for a long stretch. It does not by itself explain why the final total-step difference, after they eventually diverge, is almost always the small, k-independent value L rather than scaling with the (growing) length of the shared prefix. That part remains an open empirical pattern — see §5.

## 4. The 100%-clean family

Numbers of the alternating form x = 1010...10 (L bits) satisfy an exact algebraic identity:

```
x / 2 = (2^L − 1) / 3          (integer, since x/2 has this closed form)
3·(x/2) + 1 = 2^L               (exact — lands precisely on 2^L)
```

I.e. x's own Collatz trajectory reaches exactly 2^L (the power of two matching its own bit length) within 2 steps, then descends to 1 by pure halving. This is a much sharper condition than merely "reaches some power of two quickly" — a control test using numbers that reach unrelated powers of two (found by tracing the Collatz map backward from 2^m for m ≠ L) does not reproduce the 100% effect (only 43–70%, same as generic numbers, see exact_power_test.py).

The self-referential condition — landing exactly on 2^L where L is the block's own length — appears to align perfectly with the shift amount L used when constructing N_{k+1}, eliminating the carry/interference that causes deviation from diff_k = L for generic blocks.

## 5. Open questions

- Is there a clean closed-form proof that Pr(diff_k = L) → 1 as k → ∞ for typical x (not just the 100% family)? A partial mechanism is sketched in §3, but the exact limiting frequency and its dependence on the bit-pattern of x is not derived, only observed.
- Is 3(x/2)+1 = 2^L sufficient and necessary for the 100% effect, or are there other unrelated algebraic identities that also produce it?
- Does the phenomenon generalize to bases other than 2^L (i.e. repeating a block in base b for other b), or to other an+b Collatz-like maps?
- (see §5a) Is the self-similar block-repeat construction even necessary for the effect, or is something much more general going on?

### 5a. Generalization: the self-similarity turns out not to matter

A follow-up experiment shows the phenomenon has nothing to do with the self-similar repunit-like construction of N_k. Fix any odd x (bit length L) and take a uniformly random large integer m (not built from x at all). Then:

```
diff = steps(2^L * m + x) - steps(m)
```

still has L as its single most frequent value, and this mode sharpens (the frequency of diff = L increases) as the bit length of m grows — see random_m_test.py. Sample results (200 random trials each):

| x (bit length L) | m bit length | freq(diff = L) |
|---|---|---|
| alternating, L=12 | 100 | 75.0% |
| alternating, L=12 | 500 | 89.0% |
| random, L=12 | 100 | 32.5% |
| random, L=12 | 500 | 59.0% |
| all-ones, L=12 | 100 | 18.0% |
| all-ones, L=12 | 500 | 53.5% |

Since steps(2^L·m) = L + steps(m) exactly and trivially (multiplying by 2^L just prepends L guaranteed halving steps before the trajectory of m continues unchanged), this reduces the whole phenomenon to a cleaner question:

For a fixed small offset x and a large "clean" (2^L-divisible) number 2^L·m, why does adding x leave the total stopping time unchanged more often than any other specific outcome?

Connection to a known (but apparently unproven) fact. The case L=1, x=1 is exactly the classical, previously-documented observation that steps(n) and steps(n+1) coincide roughly half the time (OEIS A006577, comment: "It seems that about half of the terms satisfy a(i) = a(i+1)"; up to 10,000,000, 4,964,705 of the terms satisfy this). Concretely: 2m and 2m+1 are literally consecutive integers, and steps(2m) = 1 + steps(m) trivially, so "diff = L = 1" here is precisely "steps(2m) = steps(2m+1)", i.e. a same-total-stopping-time pair of neighbors — verified directly in random_m_test.py.

So the family of questions here — parameterized by (L, x) instead of just "distance 1" — appears to be a genuine generalization of that classical, apparently-still-unproven ~50% phenomenon. I could not find an existing rigorous proof or heuristic derivation of even the base L=1 case in the literature searched. (Terras's theorem characterizes the ordinary stopping time — first drop below the starting value — by finite congruence classes, but total stopping time is not known to reduce to finite congruence data, which is presumably why this is hard.)

### 5b. A proved partial result for the L=1 base case

Pushing on §5a's classical n, n+1 case directly: a clean, fully provable mechanism explains a substantial chunk of the ~50% agreement rate.

Theorem (verified). For all n ≡ 4 (mod 8) with n ≥ 12, write n = 8j+4. Then:

```
n:    8j+4 → 4j+2 → 2j+1 → 6j+4        (3 steps)
n+1:  8j+5 → 24j+16 → 12j+8 → 6j+4     (3 steps)
```

Both trajectories reach the exact same value 6j+4 after exactly 3 steps. Since the Collatz map is deterministic, everything after that point is identical, so steps(n) = 3 + steps(6j+4) = steps(n+1) exactly, with no exceptions (the sole edge case, n=4, terminates at 1 before completing 3 steps).

Generalization. Call a residue class r (mod 2^k) a merging class if n and n+1 provably reach the same value after the same fixed number of steps for every n ≡ r (mod 2^k) (checked computationally against multiple representatives per class; see merging_residue_classes.py). The fraction of merging classes, as k grows:

| modulus 2^k | merging classes | fraction |
|---|---|---|
| 8 | 1 | 12.5% |
| 16 | 3 | 18.8% |
| 32 | 8 | 25.0% |
| 64 | 18 | 28.1% |
| 128 | 39 | 30.5% |
| 256 | 82 | 32.0% |
| 512 | 170 | 33.2% |
| 1024 | 351 | 34.3% |
| 2048 | 721 | 35.2% |
| 4096 | 1476 | 36.0% |
| 8192 | 3012 | 36.8% |
| 16384 | 6130 | 37.4% |
| 32768 | 12450 | 38.0% |
| 65536 | 25241 | 38.5% |
| 131072 | 51105 | 39.0% |
| 262144 | 103358 | 39.4% |
| 524288 | 208840 | 39.8% |
| 1048576 | 421643 | 40.2% |
| 2097152 | 850737 | 40.6% |
| 4194304 | 1715546 | 40.9% |
| 8388608 | 3457791 | 41.2% |
| 16777216 | 6966495 | 41.5% |
| 33554432 | 14030369 | 41.8% |
| 67108864 | 28247507 | 42.1% |
| 134217728 | 56854178 | 42.4% |

Important correction, found after more computation (thanks to a reader running this on their own machine, up to k=27). The increments do not shrink at a constant geometric ratio. increment(k) × k is close to constant (≈7.2, slowly drifting down) across k=18..27, which looks more like increment(k) ≈ C/k — i.e. logarithmic-type growth, frac(k) ≈ frac(k0) + C·ln(k/k0), not geometric convergence to a fixed limit. A naive geometric extrapolation (which an earlier version of this README used to guess "44-46%") is therefore not reliable, and the true asymptotic behavior of this fraction as k → ∞ is currently unknown.

A further complication: the classical "~50%" itself is not a settled constant. Checking steps(n) == steps(n+1) directly for n up to 10,000,000 (not just the merging-class abstraction) shows the raw agreement rate drifting steadily upward with n — 45.0% below 100k, 47.7% below 1M, 49.6% below 10M, and already past 50% (50.4%) in the window n ∈ [9,000,000, 10,000,000]. So the OEIS comment's "~50%" appears to be a snapshot of a quantity that was itself still climbing at the point it was measured, not a stable limiting value.

Comparing the two quantities directly at matching scale (same k, i.e. N = 2^k vs modulus 2^k) shows they are related but not equal, and the gap between them is not constant either:

| k | N=2^k | raw agreement rate | merging-class fraction | gap |
|---|---|---|---|---|
| 16 | 65,536 | 44.52% | 38.51% | 6.01pt |
| 18 | 262,144 | 46.31% | 39.43% | 6.88pt |
| 20 | 1,048,576 | 47.77% | 40.21% | 7.56pt |
| 22 | 4,194,304 | 48.95% | 40.90% | 8.05pt |
| 24 | 16,777,216 | 50.03% | 41.52% | 8.50pt |

The fraction of raw agreement explained by provable merging actually decreases slightly with scale (86.5% → 83.0% across this range), which is the opposite of what I originally hoped to find.

Honest bottom line. Both the classical adjacent-agreement rate and the merging-class fraction appear to drift upward indefinitely (at least within the computationally reachable range, k up to 27), without settling at any constant found so far, and a naive extrapolation is not trustworthy given the apparent logarithmic-type growth. Whether either quantity has a well-defined limit at all — and if so, what it is, and whether the two are asymptotically related by a clean formula — is, as far as I can tell, a genuinely open question. This investigation stops here for now with that question explicitly unresolved, rather than with a (likely wrong) guessed number.

### 5c. Synthesis: it's all about how "simple" x's own trajectory is

Extending §5b's merging-class analysis from the classical L=1, x=1 case to general (L, x) ties §4 and §5b together into one picture (see general_merging_test.py):

| x | own trajectory | merging-class fraction (up to mod 512) |
|---|---|---|
| alternating, x=2730 (L=12) | reaches 2^L in 2 steps (§4) | converges cleanly to exactly 1/2 (count = 2^(k-1)-1 for modulus 2^k) |
| generic, x=2905 (L=12) | no special structure | zero merging classes found even at modulus 512 |

So the same property that made the alternating family special in §4 — its own Collatz trajectory being unusually short and predictable — is also exactly what makes algebraically-provable merging classes abundant and easy to find for that x. For "generic" x with a messy, unpredictable own trajectory, provable merges are much rarer (or require far larger moduli to detect), and the observed agreement between steps(m) and steps(2^L m+x) for such x is presumably dominated by the harder, still-unexplained "coincidental" mechanism from §5b rather than by literal trajectory merging.

This suggests a single underlying informal principle across §2–§5b:

The more predictable/short x's own Collatz trajectory is, the more of the diff = L phenomenon can be explained by provable algebraic merging, and the higher its observed frequency ceiling.

This is stated as an informal empirical pattern, not a theorem — turning it into one (e.g. a precise statement relating some complexity measure of x's own trajectory to the growth rate of the merging-class fraction) is the natural next step, left open here.

**Note (see Addendum §7):** a follow-up test found that this "simple own trajectory" principle, as tested directly on individual n (rather than on a fixed repeated block x), does *not* transfer — n's own total stopping time shows essentially zero correlation with merge speed against n+1. What does transfer is a related but distinct idea: *local* 2-adic bit structure near the LSB, rather than the global shape of the trajectory. See Addendum §7 for details.

### 5d. The classical "~50%" is not the limit — it approaches 100%

Switching from exhaustive enumeration to random sampling (much cheaper, lets us reach vastly larger scales) resolves the open question from §5b/§5c in an unexpected direction. Sampling steps(n) == steps(n+1) for random n of increasing bit length (see large_scale_sampling.c, GMP + OpenMP, run by a reader on their own machine):

| bits | ~magnitude | agreement rate |
|---|---|---|
| 10 | 10^3 | 21.3% |
| 53 | 10^16 | 49.6% |
| 264 | 10^79 | 68.6% |
| 1299 | 10^391 | 83.4% |
| 3756 | 10^1131 | 89.4% |
| 10857 | 10^3268 | 94.3% |
| 31377 | 10^9445 | 96.2% |
| 90680 | 10^27297 | 98.0% |

The agreement rate climbs well past 50% and keeps going, reaching 98% by n ~ 10^27297. Tracking 1 − rate (distance from 100%) across these points, the ratio between consecutive levels stays roughly constant (~0.7–0.85, itself not drifting toward 1 the way §5b's merging-class increments did), consistent with a power-law approach to 100% (roughly 1 − rate ~ C / sqrt(bits)) rather than the slow logarithmic growth seen in the merging-class fraction.

Revised conclusion. The OEIS A006577 comment's "about half of the terms satisfy a(i) = a(i+1)" appears to be describing a transient, small-n snapshot of a quantity that is not settling near 50% at all — it looks like Pr(steps(n) = steps(n+1)) → 1 as n → ∞. This is consistent with (and arguably the L=1, x=1 special case of) the general pattern from §5a: for a fixed shift, agreement becomes more and more likely as the numbers involved get larger. The ~50% figure that motivated this whole investigation was, in hindsight, simply the value at a scale where the convergence to 100% is still in its early stages.

This is sampling-based, not exhaustive, so it isn't a proof — but the trend across nearly 30,000 orders of magnitude, with a consistent power-law-looking decay rate, is hard to explain as a plateau at some value below 100%.

**Note (see Addendum §2):** independent exhaustive verification (up to N=6×10⁷) and independent sampling (up to 65536 bits) both corroborate this rising trend and additionally show that, at every scale tested, agreement is *fully* explained by literal early coalescence of the two trajectories (zero exceptions found) — see the Addendum for the coupling-theoretic framing and further findings, including a correction of a definitional error made along the way.

## 6. Relation to existing theory

This sits inside the well-studied 2-adic extension of the Collatz map (Bernstein & Lagarias 1996; see also the notion of "Collatz cyclic numbers" associated with periodic parity vectors). As k → ∞, N_k converges 2-adically to the rational 2-adic integer α = -x/(2^L - 1). This connects the present observation to known machinery, but the specific claim about the frequency of diff_k = L, as far as I can tell, is not stated in the literature I searched (Bernstein–Lagarias 1996; Lagarias periodicity conjecture papers; Hercher 2023; Eliahou 1993; Simons & de Weger). Corrections/pointers to prior work are very welcome.

## 7. Reproducing this

All experiments are pure Python (only the fractions and collections standard library modules) plus one C/GMP program for larger-scale exhaustive cycle searches (unrelated negative-result side-quest, kept for completeness — see cycle_search.c).

```
python3 collatz_block_repeat.py      # main phenomenon, single block
python3 big_survey.py                 # 105-block systematic sweep
python3 exact_power_test.py           # isolates the "exact 2^L" condition
python3 random_m_test.py              # §5a: self-similarity isn't necessary
python3 merging_residue_classes.py    # §5b: classifies merging residue classes
python3 general_merging_test.py       # §5c: extends §5b to general (L, x)
```

No external dependencies beyond CPython 3.8+.

## 8. Caveats

- This does not bear on the truth of the Collatz conjecture itself.
- All claims here are backed by finite computation (k up to a few thousand, L up to 24), not proofs, except where explicitly marked "proved."
- This investigation was done collaboratively with an AI assistant (Claude); all code was executed and its output verified before being reported here. Please independently re-verify before citing.

## Acknowledgments

Investigation prompted by tracing down and correcting an AI-hallucinated claim (a different model asserted a specific number was a "periodic point" of the Collatz map via a fabricated formula, which did not survive direct computational verification). The real phenomenon documented here was found in the process of checking that claim.

---

# Addendum: Coupling-Theoretic Follow-up to §5

Status: amateur/independent investigation, computationally verified, not peer-reviewed. Continuation of the original README's §5 (open questions on the `steps(n) == steps(n+1)` agreement rate). This addendum reframes the problem using probabilistic **coupling theory** (coalescing Markov chains / mixing-time arguments) and reports new computational findings, including corrections of errors made mid-investigation.

## 0. Motivation

The original README left an open question at §5c: is the "raw" agreement rate `P(steps(n) = steps(n+1))` fully explained by literal trajectory merging, or is some of it "coincidental" (both n and n+1 reaching 1 in the same number of steps via genuinely different intermediate values)? §5b/§5d's merging-class analysis suggested only a partial (and slowly, unclearly growing) fraction was explained this way.

This addendum tests that question directly by simulation rather than by exhaustive residue-class enumeration, framing it in coupling-theoretic language: define

    tau_couple(n) = min { t >= 1 : T^t(n) = T^t(n+1) }

using the raw (unaccelerated) Collatz map T, searched strictly before either trajectory reaches 1. This is the "coalescing coupling" of the two trajectories. §5b's merging-class theorem is a special case: it proves `tau_couple` is uniformly bounded (= 3) on the whole residue class n ≡ 4 (mod 8).

## 1. A definitional trap (documented for honesty)

An initial version of this experiment defined `tau_couple` by letting the raw map run indefinitely, *including after reaching 1* (where it cycles 1 → 4 → 2 → 1 → ...). This is wrong: it conflates literal early merging with accidental phase-alignment in the post-convergence cycle, and gave a spuriously high "100% coupling" reading that mixed two different phenomena.

A second version fixed the cycling issue but capped the search at `steps(n) + small_buffer`, which is a **tautology**: if `steps(n) = steps(n+1) = s`, then by definition `T^s(n) = T^s(n+1) = 1`, so a merge is *guaranteed* to be found by time `s` regardless of whether anything interesting happened earlier. Reporting "100% merge fraction" from this setup is not a finding — it restates the agreement condition.

**Corrected definition:** search only `t ∈ [1, s)`, strictly before either trajectory reaches 1. A merge found here is a genuine early coincidence of values, not a trivial simultaneous arrival at 1. All results below use this corrected definition.

## 2. Result: early merging appears to fully explain agreement, at every scale tested

**Exhaustive check** (Python), `n = 1 .. N`, computing `steps(n)`, `steps(n+1)`, and (when they agree) searching for a strict early merge:

| N | agreeing pairs | early-merge fraction | exceptions |
|---|---|---|---|
| 1,000,000 | 477,245 | 100.0000% | 0 |
| 5,000,000 | 2,454,559 | 100.0000% | 0 |
| 15,000,000 | 7,492,334 | 100.0000% | 0 |
| 60,000,000 | 30,547,761 | 100.0000% | 0 |

**Windowed random sampling** (C + GMP, `coupling_scaling.c`), single n per bit-length window:

| bits | samples | agree_frac | early_merge_frac | mean τ |
|---|---|---|---|---|
| 32 | 3000 | 40.87% | 100.0000% | 44.4 |
| 4096 | 2000 | 90.20% | 100.0000% | 2296.5 |
| 16384 | 800 | 95.63% | 100.0000% | 5011.9 |
| 65536 | 150 | 98.00% | 100.0000% | 13899.4 |

Zero exceptions across ~30.5 million exhaustively-checked pairs *and* every sampled bit-length window up to 65536 bits (n ~ 10^19728). No "coincidental, non-merging" agreement was found anywhere. The agreement-rate values also cross-check cleanly against the original README's §5d table (e.g. 98.0% at ~65536 bits, matching). If this holds in general:

    steps(n) = steps(n+1)   <=>   tau_couple(n) < steps(n)

i.e. agreement of total stopping times appears *equivalent* to (not just partially explained by) literal early coalescence — collapsing §5c's distinction between "provable merging" and "coincidental agreement" into a single mechanism, at every scale tested so far.

**Caveats:** this is exhaustive/sampled verification, not a proof; a counterexample could exist beyond the tested range. Mean early-merge time grows with scale (see §4 below), so the tail is getting heavier, not disappearing.

## 3. Coupling-time distribution: exponential-looking tail, but it's a two-component mixture

Conditional on early merging, a naive single-exponential fit to `P(tau_couple > t | merged)` looks clean at fixed small scale (N=300,000: γ≈0.0421, R²=0.9987), with minimum observed merge time 3 (matching the §5b theorem exactly). But this single-exponential picture turns out to be too simple — see §5.

## 4. γ(bits) scaling: roughly γ ∝ 1/bits, but with a caveat

Fitting the tail decay rate γ per bit-length window (log-linear regression on the survival function, windows 32–65536 bits):

| bits | γ | γ×bits |
|---|---|---|
| 32 | 0.02049 | 0.656 |
| 512 | 0.00098 | 0.502 |
| 4096 | 0.00013 | 0.532 |
| 16384 | 0.00003 | 0.492 |
| 65536 | 0.00001 | 0.655 |

Log-log regression of γ against bits gives exponent **β = −1.018** (the γ×bits column above is roughly flat across 3 orders of magnitude), consistent with γ ~ 1/bits.

However, fitting `mean(tau_couple)` against bits directly gives exponent **α = 0.728** (R² = 0.992) — not the α ≈ 1 a pure single-rate exponential model (mean = 1/γ) would predict. This mismatch was the motivation for §5.

## 5. Resolving the mismatch: the distribution is a two-component mixture

Quantile analysis (p10, p50, p90, p99, max of tau_couple) per bit-length window, each fit to a power law vs. bits:

| quantile | power-law exponent α | R² |
|---|---|---|
| p10 | **0.058** (essentially flat) | 0.738 |
| p50 (median) | 0.384 | 0.909 |
| p90 | 0.661 | 0.962 |
| p99 | **0.935** (nearly linear) | 0.998 |
| max | 0.951 | 0.9998 |

The exponent increases smoothly from ~0 to ~1 across the quantile range. Concretely, p10 is **numerically almost constant** (6.0, 6.0, 6.0, 6.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.8, 9.0) across the *entire* 32-to-65536-bit range — an 11-order-of-magnitude change in n's scale with essentially no change in the fastest 10%'s merge time. Meanwhile p99 and the max grow almost linearly with bits.

**Conclusion:** `tau_couple` is not a single exponential. It is a mixture of (a) a bits-independent fast component (roughly the bottom ~90% of the distribution, driven by small, fixed residue classes — see §6) and (b) a heavy tail whose scale grows roughly linearly with bit length. The previously-observed mean exponent (0.73) is what you get when averaging over a mixture whose bulk is O(1) and whose tail is O(bits); it sits strictly between the two component behaviors, as expected for a mixture, and is not itself a fundamental exponent.

## 6. Hypothesis 1 (confirmed): the fast component is explained by tiny, scale-independent residue classes

Direct test: for 150 "fast" mergers (tau_couple < 20) sampled at a fixed, large bit length (4096 bits), binary-search the smallest modulus 2^K such that fixing n's low K bits (with fresh random high bits) reliably reproduces the same merge time across multiple independent trials.

Result: **all 150/150 samples were fully explained by K ≤ 14** (mean K = 5.23):

| K | count |
|---|---|
| 3 | 48 |
| 4 | 24 |
| 5 | 29 |
| 6 | 15 |
| 7 | 11 |
| 8–14 | 23 |

K=3 (modulus 2^3=8) is the single most common value, at 48/150 (32%) — this is *exactly* the §5b theorem's residue class n ≡ 4 (mod 8), now shown empirically to be the dominant single mechanism behind the fast component of the distribution, not merely an isolated example. This directly confirms that §5's "fast, bits-independent" component is driven by small, fixed residue classes generalizing the §5b mechanism, matching the near-flat p10 behavior found in §5.

A related control test (fixing n's low 24 bits and checking whether merge time is reproduced) shows a sharp transition: merges with tau_couple below ~24–40 are ~100% reproduced by 24 bits of information; merges above that are essentially never reproduced. This is expected given the classical fact that a trajectory's first t steps are determined by n mod 2^t (Bernstein–Lagarias) — merges taking t steps generically need on the order of t bits to pin down — but confirms the mechanism is tight (no large hidden dependence on far-away bits for the fast component).

## 7. Hypothesis 2 (mostly refuted, with a clean positive finding instead): n's *own* trajectory length is irrelevant; local 2-adic structure is what matters

Tested whether §5c's principle ("the simpler x's own trajectory, the more merging is explained") transfers to individual n — i.e. does a short/simple own total stopping time `steps(n)` predict fast merging with n+1? Sampled 4000 pairs at 2048 bits, all with `steps(n)=steps(n+1)`:

- `corr(steps(n)/bits, log(tau_couple)) = +0.020` — **essentially zero**. Fastest/median/slowest 10% groups all have the same mean `steps(n)/bits ≈ 7.2`, matching the well-known generic ~7 steps/bit average with no distinguishing signal.
- Instead, **local 2-adic structure at the low-order bits** correlates meaningfully: `corr(v2(n+1), log τ) = +0.363`, `corr(run_length(n+1), log τ) = +0.289` (v2 = 2-adic valuation / trailing-zero count).
- The fastest 10% of all samples had **exactly** v2(n) = 2.000 (zero variance) and merge_time exactly 3 — i.e. the fastest decile is entirely the n ≡ 4 (mod 8) class, consistent with §6.

**An apparent asymmetry** — v2(n+1) correlating much more strongly than v2(n) (+0.363 vs. −0.021 pooled over all samples) — turned out to be a **pooling artifact**, not a real effect: since exactly one of n, n+1 is even, v2(n)=0 whenever n is odd and vice versa, so pooling both cases together dilutes whichever variable is "trivially zero" in each half. Comparing like-for-like (within the n-even subgroup only, where v2(n) is meaningful, vs. within the n-odd subgroup only, where v2(n+1) is meaningful) gives symmetric correlations: **+0.322 and +0.306 respectively** — the same effect, as expected by the n ↔ n+1 symmetry of the setup.

Partial correlation of `steps(n)/bits` against `log(tau_couple)`, controlling for v2(n+1) (and separately for run-length), stays near zero (+0.016, +0.015) — confirming the null result is not merely masked by the v2 confound. Stratifying by v2(n+1) value and computing within-stratum correlation of `steps(n)` with merge time also stays near zero in every stratum.

**Conclusion:** merging speed is governed by *local* low-order bit structure (2-adic valuation / run-length near the LSB, i.e. exactly the kind of information the Lagarias congruence fact says determines early steps), not by any *global* property of n's own trajectory shape. §5c's "simple trajectory" principle, formulated for the specific repeated-block construction of the original README (§1–§4), does not appear to generalize to arbitrary (n, n+1) pairs in the form tested here.

## 8. Tooling produced this session

- `coupling_experiment3.py` / `4.py` / `5.py` — Python, exhaustive verification up to N=6×10⁷ (memoized total-stopping-time cache), used for §2.
- `coupling_scaling.c` — C + GMP + OpenMP, windowed random sampling across bit-length scales (32–65536 bits, extensible), outputs per-window CSV of `(agree, merge_time)`. Build: `gcc -O3 -fopenmp -o coupling_scaling coupling_scaling.c -lgmp -lm`.
- `fit_gamma.py` — fits γ per bit-length window from the C program's output (log-linear regression on the survival-function tail).
- `fast_slow_test.py` — tests whether merges below/above a threshold are reproduced by fixing n's low K bits (§6, sharp-transition control test).
- `minimal_K_test.py` — binary search for the smallest modulus 2^K explaining each fast merger (§6 main result).
- `hypothesis2_test.py` — correlates merge time against n's own trajectory length and 2-adic/run-length features at fixed bit length; produces the partial-correlation and stratified-correlation analysis (§7).

## 9. Summary: established vs. still open

**Established (exhaustive up to N=6×10⁷; sampled up to 65536 bits with zero exceptions):**
- Every observed `steps(n)=steps(n+1)` pair coalesces to a literal common value strictly before reaching 1, at every scale tested.
- `tau_couple` is a two-component mixture: a bits-independent fast bulk (~90%, explained by tiny fixed residue classes, dominated by the §5b n≡4(mod 8) mechanism) plus a heavy tail growing roughly linearly in bit length (~top 1%, p99 exponent ≈ 0.94).
- Merging speed depends on *local* 2-adic structure near the LSB of n (and n+1), not on any global property of n's own Collatz trajectory length; the apparent v2(n) vs v2(n+1) asymmetry is a pooling artifact and disappears under like-for-like comparison.

**Open:**
- Whether the 100% early-merge / zero-exception pattern holds for all n (currently: exhaustive to 6×10⁷, sampled with zero exceptions to 65536 bits).
- A precise theoretical model for the two-component mixture — e.g. can the fast/slow split point and the p99 ≈ bits^0.94 exponent be derived (rather than fit) from the finite-state "carry" structure of the 3x+1 map read bit-serially (Bernstein–Lagarias)?
- Why the fast component's modulus needs (§6, K values 3–14) cluster where they do, and whether this connects to §4's "exact 2^L" family in the original README — **partially answered, see §10.**

## 10. Follow-up (2026-07-26): independent re-verification, a methodology caveat, and the §4-connection question resolved (negatively)

**Note on provenance.** The tooling behind §6–§8 (`coupling_experiment3/4/5.py`, `coupling_scaling.c`, `fit_gamma.py`, `fast_slow_test.py`, `minimal_K_test.py`, `hypothesis2_test.py`) was written and run in a session whose source files were never saved to disk — only the README text and result tables survived. Separately, four files referenced since §5b/§5c/§5d (`merging_residue_classes.py`, `merging_classes_fast.c`, `general_merging_test.py`, `large_scale_sampling.c`) existed locally but had never actually been uploaded to this repository; they have now been added. The lost §6–§8 tooling has not been reconstructed file-for-file, but its central claim is re-derived independently below using a new script, `k_cluster_analysis.py`.

**Re-derivation of §6's K-distribution.** Sampling 400 random 2000-bit n and isolating fast mergers (tau_couple < 20) reproduces §6's shape closely: K=3 is the single most common minimal modulus (47/145 = 32.4%, vs. the original 32%), with counts falling off through K=14, matching the original run's K ≤ 14 range. This is independent evidence the original §6 finding was not a fluke of that particular (lost) implementation.

**New finding: the minimal-K binary search can false-positive.** Cross-checking each minimal-K residue against the exhaustive, provably-correct classifier from `merging_residue_classes.py` (which checks *every* representative up to the modulus, not a handful of random trials), 129/133 (97%) of K≤9 cases matched exactly — but 4 did not, all of the form "K=3, r≡2 (mod 8)" merging at t=5. The exhaustive classifier shows r≡2 (mod 8) is *not* actually a uniform merging class (only r≡4 (mod 8) is, per §5b) — r≡2 (mod 8) splits into a merging sub-class (r≡2 mod 16) and a non-merging one (r≡10 mod 16). With only ~5 random trials per K in the binary search, all 4 exceptions happened to sample only the merging sub-class by chance, so the search reported "K=3 suffices" when the true minimal K was 4. **Practical implication:** §6's reported K-values (and by extension its mean K≈5.23) are very slightly biased low by this sampling artifact; the effect is small (3% of cases here) but not zero, and a rerun with more trials per K would be needed to fully correct it.

**The §4-connection question, resolved negatively (beyond the trivial cases).** For each fast-merging residue r found above, checking how many steps r's own trajectory takes to first hit a power of two: the two smallest classes, r=4 (K=3) and r=2 (K=4), hit a power of two in 0 steps — trivially, since 4 and 2 already *are* powers of two. r=5 (K=5) reaches 2^4 in 1 step, consistent with the §4 "exact 2^L" family (5 = `101`₂ is itself an alternating-block number). But every larger-K class checked (r=14, 22, 49, 65, 99, ...) takes 11–23 steps to reach any power of two — no faster or cleaner than a generic number. So the K-clustering is **not** generally explained by §4's mechanism: only the smallest 1–2 classes coincide with it, apparently because the smallest possible residues are close to powers of two by sheer smallness, not because of the §4 mechanism itself. This is consistent with, and sharpens, Addendum §7's separate finding that global trajectory simplicity doesn't predict merge speed — here we see directly that it doesn't explain the K-clustering either.

See `k_cluster_analysis.py` for the full script (self-contained, pure Python, run: `python k_cluster_analysis.py`).

## 11. Follow-up: independent re-verification of the quantile-scaling mixture model (§5)

Addendum §5's claim — that `tau_couple`'s quantiles scale as power laws in bit length, with the exponent rising smoothly from ~0 at p10 to ~1 at p99 — is another claim whose original tooling (`coupling_scaling.c`, `fit_gamma.py`) was lost (§10). Re-derived from scratch in pure Python (`quantile_scaling_analysis.py`), sampling 400 random n per bit-length window across 64–4096 bits (a narrower range than the original's 32–65536 bits, since this uses plain Python big-int simulation rather than GMP/OpenMP, but wide enough — a factor of 64 — to see the trend):

| bits | p10 | p50 | p90 | p99 | max | agreement rate |
|---|---|---|---|---|---|---|
| 64 | 3 | 13 | 181 | 369 | 453 | 63.0% |
| 128 | 3 | 14 | 502 | 774 | 894 | 69.0% |
| 256 | 3 | 18 | 678 | 1396 | 1808 | 76.2% |
| 512 | 3 | 30 | 1657 | 3319 | 3773 | 80.8% |
| 1024 | 3 | 60 | 2723 | 6479 | 6928 | 87.2% |
| 2048 | 3 | 58 | 4982 | 13065 | 14835 | 88.5% |
| 4096 | 3 | 79 | 6438 | 25712 | 29626 | 92.2% |

Power-law exponents (quantile ~ bits^alpha) from log-log regression, compared to the original (lost-code) values:

| quantile | this re-run (alpha, R²) | original Addendum §5 |
|---|---|---|
| p10 | **0.000**, exactly flat (literally 3 at every window tested — even flatter than originally reported) | 0.058 |
| p50 | 0.487, R²=0.944 | 0.384 |
| p90 | 0.860, R²=0.979 | 0.661 |
| p99 | 1.026, R²=0.999 | 0.935 |

The exponents are not numerically identical to the original run (expected: this samples a 64x bit-length range vs. the original's 2048x range, and only 400 samples/window vs. presumably more), but the **qualitative structure is fully reproduced independently**: p10 is exactly bits-independent (here, literally constant at 3 — the n≡4 (mod 8) mechanism from §5b, its minimum possible value), and the exponent rises monotonically through p50 and p90 to essentially linear (alpha≈1) at p99. This is strong independent support for the two-component mixture picture — an O(1) fast bulk plus a heavy tail scaling roughly linearly in bit length — even without the original code surviving.

See `quantile_scaling_analysis.py` (self-contained, pure Python, run: `python quantile_scaling_analysis.py`).

**Update: re-run at the original's full 64–65536 bit range (`quantile_scaling_large.py`).** The initial 64–4096 bit re-run above was intentionally narrow (a factor of 64 in scale) and its exponents came out systematically higher than the original's, most likely from that narrow range. Re-running with the same six-window design but spanning the original's full 32768x range (64, 256, 1024, 4096, 16384, 65536 bits; 800→120 samples, shrinking at larger scale to keep runtime reasonable — the whole run took 65s):

| bits | p10 | p50 | p90 | p99 | max | agreement rate |
|---|---|---|---|---|---|---|
| 64 | 3 | 11 | 193 | 419 | 462 | 62.8% |
| 256 | 3 | 17 | 664 | 1396 | 1841 | 75.2% |
| 1024 | 3 | 61 | 2613 | 6028 | 7495 | 84.8% |
| 4096 | 5 | 92 | 4244 | 20115 | 20898 | 91.3% |
| 16384 | 3 | 67 | 13875 | 77027 | 92525 | 97.5% |
| 65536 | 3 | 124 | 11926 | 274650 | 354391 | 98.3% |

| quantile | this re-run (full range) | narrow-range re-run (above) | original Addendum §5 |
|---|---|---|---|
| p10 | 0.011, R²=0.017 (flat, noise-level) | 0.000 | 0.058 |
| p50 | **0.343**, R²=0.837 | 0.487 | **0.384** |
| p90 | **0.623**, R²=0.933 | 0.860 | **0.661** |
| p99 | **0.941**, R²=1.000 | 1.026 | **0.935** |

At matching scale, the exponents land almost exactly on the original (lost-code) values — p50 within 0.04, p90 within 0.04, p99 within 0.006 — confirming the narrow-range run's higher exponents were indeed a range artifact, not a real discrepancy. This is about as clean an independent reproduction of a "lost" numerical result as this project has managed: different code, different session, matching to 2-3 significant figures. See `quantile_scaling_large.py`.

## 12. Follow-up: independent re-verification of the exhaustive §2 claim

Unlike §5b–§5d, the Addendum's §2 exhaustive table (N=1e6/5e6/15e6/60e6, 100.0000% early-merge, 0 exceptions) had **no surviving verification code at all** — `coupling_experiment3/4/5.py` were never saved to disk (§10). Reconstructed from scratch (`coupling_exhaustive_verify.py`, memoized total-stopping-time cache) and re-run at every N from the original table, plus an extension beyond it:

| N | agreeing pairs (this re-run) | agreeing pairs (original) | early-merge | exceptions | wall time |
|---|---|---|---|---|---|
| 1,000,000 | 477,245 | 477,245 | 100.0000% | 0 | 2.0s |
| 5,000,000 | 2,454,559 | 2,454,559 | 100.0000% | 0 | 11.6s |
| 15,000,000 | 7,492,334 | 7,492,334 | 100.0000% | 0 | 36.5s |
| 60,000,000 | 30,547,761 | 30,547,761 | 100.0000% | 0 | 2m27s |
| 100,000,000 | 51,242,281 | *(new — extends past the original's largest tested N)* | 100.0000% | 0 | 4m20s |

Every single agreeing-pair count from the original (lost-code) table matches **exactly**, digit for digit, and the zero-exceptions / 100.0000% early-merge result reproduces cleanly at every scale, including the new N=100,000,000 point which goes beyond anything in the original table. This closes the biggest provenance gap left by the lost tooling: §2's central claim is no longer resting on unreproducible code, and now has independent exhaustive verification to a new high-water mark of N=10^8.

See `coupling_exhaustive_verify.py` (self-contained, pure Python, memoized; run: `python coupling_exhaustive_verify.py <N>`).

## 13. New finding: the phenomenon generalizes to powers of 3 (answering the original README's §5 open question)

The original README's §5 left open: *"Does the phenomenon generalize to bases other than 2^L ... or to other an+b Collatz-like maps?"* This had never been tested. It now has an answer.

**Control: generic non-power-of-2/3 multipliers show no effect.** Building `N_{k+1} = c*N_k + x` for decimal blocks (c=10^3, 10^4), a large prime (c=1,000,003), and an arbitrary composite (c=12,345) — analogous to the original §1 construction but in a base unrelated to the Collatz map's own arithmetic — produces no concentration at all: `diff_k` values scatter across a wide, essentially unstructured range, with the most common single value occurring only 2-7% of the time, not increasing with more iterations. See `base_generalization_test.py`.

**Main finding: c = 3^p produces a clean, scale-increasing concentration on diff = −p.** Testing `diff = steps(3^p * m + x) - steps(m)` for random m (`power_of_3_test.py`, same methodology as §5a's `random_m_test.py`), across p = 1, 2, 3, 5 and several x per p:

| p (c=3^p) | m bits=16 | bits=64 | bits=256 | bits=1024 | bits=4096 |
|---|---|---|---|---|---|
| p=1, x=0 | 40.0% | 51.3% | 65.3% | 82.7% | 92.0% |
| p=1, x=2 | 55.3% | 72.0% | 80.0% | 87.3% | 94.0% |
| p=2, x=1 | 35.3% | 51.3% | 66.7% | 82.7% | 90.7% |
| p=3, x=1 | 14.0% | 40.0% | 54.0% | 75.3% | 87.3% |
| p=5, x=1 | 12.7% | 30.0% | 54.0% | 75.3% | 84.0% |

(freq(diff = −p), i.e. the fraction of trials landing exactly on the dominant value). Every (p, x) combination tested shows the same shape as §5a's original result: low frequency at small bit lengths, climbing steadily toward (but not yet reaching, in this bit-length range) 100% as m grows — strong evidence this is the same general phenomenon, not a coincidence specific to base 2.

**A trivial exact case, flagged honestly.** p=1, x=1 (i.e. `diff = steps(3m+1) - steps(m)`) measured *exactly* 100.0% at every bit length tested, with zero exceptions. This is **not a new discovery** — it's a restatement of the Collatz map's own definition: for odd m, `T(m) = 3m+1` by definition, so the trajectory of `3m+1` *is* m's own trajectory starting one step in, making `steps(3m+1) = steps(m) - 1` exactly and trivially for every odd m (all m tested here are odd, forced by the random-generation code). Included for completeness and to avoid the appearance of overclaiming.

**Mechanism (direct trajectory inspection, same style of evidence as the Addendum's coupling theory):** e.g. for c=9 (p=2), m=15: trajectory of m is `15, 46, 23, 70, 35, 106, 53, 160, 80, 40, ...` (reaches 40 after 9 steps), while the trajectory of `9*15+1=136` is `136, 68, 34, 17, 52, 26, 13, 40, ...` (reaches the *same* value 40 after only 7 steps). Once merged, the rest is identical, so `steps(136) = steps(15) - 2` exactly = −p. The general pattern: `3^p * m + x`'s trajectory literally coalesces into m's own trajectory, but reaches the shared point exactly p steps sooner — the mirror image of the original §3 mechanism (where `2^L * N_k + x`'s trajectory shares a *prefix* with N_k's, adding L steps) rather than merging early.

**Confirms the mechanisms combine additively.** Testing c = 6^p = 2^p·3^p (which bundles both the base-2 "add p steps" mechanism and the base-3 "save p steps" mechanism) gives `diff = 0` as the dominant value at every p tested (1–4) and every bit length — the two effects exactly cancel, as the additive mechanism would predict.

**A §5b-style partial proof: provable merging classes exist for this pairing too.** §5b showed that certain residues r (mod 2^k) provably force n and n+1 to merge at a fixed step-pair, for the classical pairing. The same test applies here to the pairing (m, 3^p·m+x): a residue r (mod 2^k) is a "merging class" if *every* representative m≡r merges with 3^p·m+x at the same fixed (step_from_m, step_from_n). Running this (`power_of_3_merging_classes.py`) for several (p, x):

| p, x | mod 8 | mod 64 | mod 256 | mod 1024 | diff values found |
|---|---|---|---|---|---|
| p=1, x=2 | 37.5% | 48.4% | 50.0% | 51.5% | **{−1}** only |
| p=2, x=1 | 0.0% | 4.7% | 7.4% | 8.9% | **{−2}** only |
| p=3, x=1 | 0.0% | 0.0% | 0.0% | 0.0% (first class at mod 512) | **{−3}** only |

Every single provable merging class found, at every (p, x, k) tested, gives exactly diff = −p — never any other value. This is the same style of rigorous (if partial) evidence as §5b's theorem: not a full proof for all m, but an actual growing family of residue classes where the −p result is provably exact, not just statistically likely. As with §5b, higher p needs a larger modulus before the first merging class appears (p=1 already has classes at mod 8; p=3 needs mod 512) — consistent with needing more bits of "setup" to force a larger, exact step-savings.

**Second control: other prime powers (5, 7, 11, 13) also show no effect.** To confirm the effect is specific to the primes 2 and 3 (the ones that actually appear in the Collatz map's arithmetic) rather than "any small prime power," the same test was run for c = base^p with base ∈ {5, 7, 11, 13}, p ∈ {1,2,3}, x ∈ {0,1}, 240 samples each (64 and 1024 bits): every single case showed no concentration at all — top values occurring in only 2-5% of trials, indistinguishable from noise, just like the decimal/prime/composite controls above.

**Conclusion.** The phenomenon is not a base-2 peculiarity. It generalizes specifically to multipliers built from the prime factors that appear in the Collatz map's own arithmetic (2 in `n/2`, 3 in `3n+1`) — with base 2 adding L steps and base 3 *removing* p steps — and does not appear for any multiplier unrelated to that arithmetic, including other prime powers with no special relationship to 2 or 3. This resolves the original README's §5 open question about base-generalization, in a more specific and interesting way than a simple yes/no: it isn't "any base b," it's specifically the two primes the map itself is built from.

See `power_of_3_test.py` and `base_generalization_test.py` (self-contained, pure Python).

## Acknowledgments (Addendum)

This addendum was developed collaboratively with an AI assistant (Claude); all code was executed and its output verified before being reported here, including catching and correcting a tautological definition error mid-investigation (§1) and a data-construction bug that initially produced a misleading result (§6). §10 was added in a later session after discovering the §6–§8 tooling files had been lost (never saved outside that session); its findings were independently re-derived with new code and cross-checked against the exhaustively-verified merging-class data from §5b. Please independently re-verify before citing.
