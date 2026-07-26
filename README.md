# Self-Similar Digit Blocks and Collatz Total Stopping Time

Status: amateur/independent investigation, computationally verified, not peer-reviewed.

**TL;DR.** Repeatedly prepending a fixed bit-block `x` (bit length `L`) to a number increases its Collatz total stopping time by an amount that converges — with probability approaching 1 — to exactly `L`. A specific family of blocks (`3·(x/2)+1 = 2^L`) makes this 100% exact from the first iteration. The same investigation also found that the classical "about half of n, n+1 share a total stopping time" folklore (OEIS A006577) is not a stable ~50% limit — it climbs toward 100% as n grows — and that this agreement is fully explained, at every scale tested, by the two trajectories literally coalescing to a common value before either reaches 1. A later follow-up found the whole phenomenon generalizes beyond base 2: multiplying by `3^p` produces a mirror-image effect (steps *decrease* by exactly `p`), tying the effect directly to the two primes (2 and 3) that appear in the Collatz map's own arithmetic.

This is not a claim about the Collatz conjecture itself (true/false, cycles, divergent trajectories). It is a set of concrete, checkable statements about specific number families and about the statistics of the total stopping time function.

**This file is a short summary.** The full investigation — every table, the coupling-theory Addendum, the honest record of mistakes made and corrected along the way, and the newest power-of-3 generalization — is in **[README_extra.md](README_extra.md)**, organized into 13 sections + Addendum. This file exists so a new reader isn't faced with 500+ lines up front.

## Established results (with strength of evidence)

| # | Result | Evidence |
|---|---|---|
| 1 | `steps(2^L·N_k + x) − steps(N_k)` concentrates on `L`, converging toward high frequency as k grows | 105-block sweep, 0 exceptions ([§2](README_extra.md#2-main-empirical-finding)) |
| 2 | Mechanism: shared bit-prefix forces identical parity vectors for the shared length | Proved ([§3](README_extra.md#3-why-this-happens-proved-part)) |
| 3 | Alternating blocks (`3(x/2)+1=2^L` exactly) give 100% exact convergence from k=1 | Proved + verified ([§4](README_extra.md#4-the-100-clean-family)) |
| 4 | Self-similarity of the construction isn't needed — same effect for `2^L·m+x`, any random `m` | ([§5a](README_extra.md#5a-generalization-the-self-similarity-turns-out-not-to-matter)) |
| 5 | n≡4 (mod 8) provably merges with n+1 in exactly 3 steps; a slowly-growing family of such "merging classes" exists | Proved for each class found ([§5b](README_extra.md#5b-a-proved-partial-result-for-the-l1-base-case)) |
| 6 | OEIS A006577's "~50%" is **not** a stable limit — raw agreement climbs to 98%+ by ~10⁵ digits | Sampled to n~10^27297 ([§5d](README_extra.md#5d-the-classical-50-is-not-the-limit--it-approaches-100)) |
| 7 | Agreement of `steps(n)=steps(n+1)` is fully explained by literal early trajectory coalescence — zero exceptions found anywhere tested | Exhaustive to N=10⁸ (independently re-verified twice), sampled to 65536 bits ([Addendum §2](README_extra.md#2-result-early-merging-appears-to-fully-explain-agreement-at-every-scale-tested), [§12](README_extra.md#12-follow-up-independent-re-verification-of-the-exhaustive-2-claim)) |
| 8 | The coupling time `tau_couple` is a two-component mixture: an O(1) fast bulk (small residue classes) + a heavy tail scaling ~linearly in bit length (p99 exponent ≈ 0.94) | Independently re-derived, matches original to 2-3 sig figs ([Addendum §5](README_extra.md#5-resolving-the-mismatch-the-distribution-is-a-two-component-mixture), [§11](README_extra.md#11-follow-up-independent-re-verification-of-the-quantile-scaling-mixture-model-5)) |
| 9 | The phenomenon generalizes to `c=3^p` multipliers: `diff = steps(3^p·m+x) − steps(m)` concentrates on **−p** (mirror image of result #1); combines additively with base 2 (`6^p` gives diff=0); specific to primes 2,3 — other prime powers (5,7,11,13) show no effect | New this session, control-tested, partial proof via merging classes ([§13](README_extra.md#13-new-finding-the-phenomenon-generalizes-to-powers-of-3-answering-the-original-readmes-5-open-question)) |

## Still open

- No closed-form proof that Pr(diff_k = L) → 1 in general (only the alternating family is proved exactly).
- No precise theoretical derivation of the two-component mixture's parameters (fast/slow split, the ≈bits^0.94 tail exponent) — currently fit, not derived.
- No proof that zero-exception early-merging holds for *all* n (only verified exhaustively/by sampling to large but finite scales).
- The −p mechanism for powers of 3 (result #9) is empirically and partially (via merging classes) confirmed but has no closed-form proof analogous to §3's for base 2.
- Whether other an+b Collatz-like maps (not just the standard 3n+1) show analogous effects — untested.

## Repository contents

```
README.md                    this summary
README_extra.md              full write-up (13 sections + Addendum)
LICENSE                      MIT
results/                     raw CSV data from §2 and §1's single-block run

Core (§1-4):
collatz_block_repeat.py, big_survey.py, exact_power_test.py

§5a-5d / Addendum (n, n+1 agreement statistics):
random_m_test.py, merging_residue_classes.py, merging_classes_fast.c,
general_merging_test.py, large_scale_sampling.c, coupling_experiment3/4/5.py*,
coupling_scaling.c*, fit_gamma.py*, fast_slow_test.py*, minimal_K_test.py*,
hypothesis2_test.py*

§10-13 follow-ups (this session's independent re-verification + new finding):
k_cluster_analysis.py, quantile_scaling_analysis.py, quantile_scaling_large.py,
coupling_exhaustive_verify.py, base_generalization_test.py, power_of_3_test.py,
power_of_3_merging_classes.py

Unrelated side-quest:
cycle_search.c               exhaustive Collatz-cycle search, negative result

* starred files are referenced in README_extra.md but were lost (never saved
  to disk) before this repo could capture them — see README_extra.md §10 for
  the provenance note and independent re-verification of their claims.
```

## Relation to existing theory

The base-2 construction sits inside the 2-adic extension of the Collatz map studied by Bernstein & Lagarias ("The 3x+1 conjugacy map," *Canad. J. Math.* 1996): as k→∞, `N_k` converges 2-adically to `-x/(2^L-1)`, and a trajectory's first m steps are determined by n mod 2^m — the fact underlying §3's proof. See [README_extra.md §6](README_extra.md#6-relation-to-existing-theory) and [§13](README_extra.md#13-new-finding-the-phenomenon-generalizes-to-powers-of-3-answering-the-original-readmes-5-open-question) for further discussion and open questions about what is/isn't covered by existing literature (citations there should be independently checked before citing — see the caveats below).

## Reproducing this

```
python3 collatz_block_repeat.py      # main phenomenon, single block
python3 big_survey.py                # 105-block systematic sweep
python3 power_of_3_test.py           # newest: base-3 generalization
```

No external dependencies beyond CPython 3.8+ for the Python scripts; the `.c` files need GMP + OpenMP. Full reproduction instructions for every experiment are in README_extra.md.

## Caveats

- This does not bear on the truth of the Collatz conjecture itself.
- All claims are backed by finite computation, not proofs, except where explicitly marked "proved."
- This investigation was done collaboratively with an AI assistant (Claude). All code was executed and its output verified before being reported, including several caught-and-corrected mistakes (documented in README_extra.md rather than hidden). Please independently re-verify before citing.

## Acknowledgments

Investigation prompted by tracing down and correcting an AI-hallucinated claim about a Collatz "periodic point." Full acknowledgments in README_extra.md.
