# Self-Similar Digit Blocks and Collatz Total Stopping Time

Status: amateur/independent investigation, computationally verified, not peer-reviewed.

**TL;DR.** Repeatedly prepending a fixed bit-block `x` (bit length `L`) to a number increases its Collatz total stopping time by an amount that converges — with probability approaching 1 — to exactly `L`. A specific family of blocks (`3·(x/2)+1 = 2^L`) makes this exact from the first iteration in every case tested — provably so for a quarter of all residues (§17), empirically for the rest. The same investigation also found that the classical "about half of n, n+1 share a total stopping time" folklore (OEIS A006577) is not a stable ~50% limit — it climbs toward 100% as n grows — and that this agreement is fully explained, at every scale tested, by the two trajectories literally coalescing to a common value before either reaches 1. A later follow-up found the whole phenomenon generalizes beyond base 2: multiplying by `3^p` produces a mirror-image effect (steps *decrease* by exactly `p`). A short coefficient-matching argument then proved *why*: whenever the two trajectories merge with fixed parity structure, the step difference is forced to be exactly `L − p`, and the multiplier must be 3-smooth (`2^L·3^p`) for merging to be possible at all — so the restriction to the primes 2 and 3 is a consequence of unique factorization, not a coincidence. That proof turns out never to use the number 3 at all: the same law holds for every `an+b` map with `a` odd, with the map's own multiplier playing 3's role. What remains genuinely open is not the *value* of the effect but its *frequency*.

This is not a claim about the Collatz conjecture itself (true/false, cycles, divergent trajectories). It is a set of concrete, checkable statements about specific number families and about the statistics of the total stopping time function.

**This file is a short summary.** The full investigation — every table, the coupling-theory Addendum, the honest record of mistakes made and corrected along the way, and the newest power-of-3 generalization — is in **[README_extra.md](README_extra.md)**, organized into 19 sections + Addendum. This file exists so a new reader isn't faced with 500+ lines up front.

## Established results (with strength of evidence)

| # | Result | Evidence |
|---|---|---|
| 1 | `steps(2^L·N_k + x) − steps(N_k)` concentrates on `L`, converging toward high frequency as k grows | 105-block sweep, 0 exceptions ([§2](README_extra.md#2-main-empirical-finding)) |
| 2 | Mechanism: shared bit-prefix forces identical parity vectors for the shared length | Proved ([§3](README_extra.md#3-why-this-happens-proved-part)) |
| 3 | Alternating blocks (`3(x/2)+1=2^L` exactly) give 100% exact convergence from k=1 | **Label corrected:** §4 proves only the identity `3(x/2)+1=2^L`, not that `diff=L` follows. §17 now proves `diff=L` for `m≡2 (mod 4)` (¼ of residues, every even L, one exception at m=2); the other ¾ remain empirical ([§4](README_extra.md#4-the-100-clean-family), [§17](README_extra.md#17-filling-the-gap-in-4-a-proof-that-the-alternating-family-actually-gives-diff--l)) |
| 4 | Self-similarity of the construction isn't needed — same effect for `2^L·m+x`, any random `m` | ([§5a](README_extra.md#5a-generalization-the-self-similarity-turns-out-not-to-matter)) |
| 5 | n≡4 (mod 8) provably merges with n+1 in exactly 3 steps; a slowly-growing family of such "merging classes" exists | Independently re-derived; **not new** — this is Garner (1985) / LaDue (2017) Corollary 5.1, see [§6](README_extra.md#6-relation-to-existing-theory) ([§5b](README_extra.md#5b-a-proved-partial-result-for-the-l1-base-case)) |
| 6 | OEIS A006577's "~50%" is **not** a stable limit — raw agreement climbs to 98%+ by ~10⁵ digits | Sampled to n~10^27297. **⚠ Priority unverified:** the "it increases" half may be in Gao (1993); paywalled, unchecked — see [§19](README_extra.md#19-cross-checking-5d-against-5b14-the-log-growth-law-provably-breaks-down) ([§5d](README_extra.md#5d-the-classical-50-is-not-the-limit--it-approaches-100)) |
| 7 | Agreement of `steps(n)=steps(n+1)` is fully explained by literal early trajectory coalescence — zero exceptions found anywhere tested | Exhaustive to N=10⁸ (independently re-verified twice), sampled to 65536 bits ([Addendum §2](README_extra.md#2-result-early-merging-appears-to-fully-explain-agreement-at-every-scale-tested), [§12](README_extra.md#12-follow-up-independent-re-verification-of-the-exhaustive-2-claim)) |
| 8 | The coupling time `tau_couple` is a two-component mixture: an O(1) fast bulk (small residue classes) + a heavy tail scaling ~linearly in bit length (p99 exponent ≈ 0.94) | Independently re-derived, matches original to 2-3 sig figs ([Addendum §5](README_extra.md#5-resolving-the-mismatch-the-distribution-is-a-two-component-mixture), [§11](README_extra.md#11-follow-up-independent-re-verification-of-the-quantile-scaling-mixture-model-5)) |
| 9 | The phenomenon generalizes to `c=3^p` multipliers: `diff = steps(3^p·m+x) − steps(m)` concentrates on **−p** (mirror image of result #1); combines additively with base 2 (`6^p` gives diff=0); specific to primes 2,3 — other prime powers (5,7,11,13) show no effect | New this session, control-tested; the `−p` value and the prime-2,3 specificity are now **proved** (#11, #12) ([§13](README_extra.md#13-new-finding-the-phenomenon-generalizes-to-powers-of-3-answering-the-original-readmes-5-open-question)) |
| 10 | Result #9's merging-class fraction (like #6/#5b's classical analogue) does not converge to a stable limit as modulus grows; unlike the classical case's near-constant `increment(k)×k`, here it climbs steadily with no sign of flattening in the tested range (k up to 23-25), and p=3 grows measurably faster than p=2 | New this session, pushed to modulus 2^25 (p=2) / 2^23 (p=3) ([§14](README_extra.md#14-follow-up-pushing-13s-merging-class-fraction-to-larger-moduli-does-it-show-the-same-log-type-growth-as-5b)) |
| 11 | **Whenever a merge happens with fixed parity structure, the diff value is forced to be exactly `L − p` — no other value is possible.** One coefficient-matching argument proves the observed `+L` (#1), `−p` (#9), and `0` (6^p) all at once | **Proved** (elementary; likely known folklore — see caveat); verified on 315 merging classes, 0 violations ([§15](README_extra.md#15-a-proof-that-the-diff-value-is-forced-one-argument-explaining-5a-13-and-the-6p-case)) |
| 12 | **The multiplier must be 3-smooth (`c = 2^L·3^p`)** for the merging mechanism to exist at all — any other prime factor makes it structurally impossible. This proves #9's control results (5,7,11,13, 10^k, primes, composites) rather than leaving them empirical. Also: the diff value is independent of the offset `x` | **Proved** (corollary of #11, no computation) ([§16a-b](README_extra.md#16a-corollary-the-multiplier-must-be-3-smooth-this-proves-13s-control-results)) |
| 13 | **None of this is about the number 3.** #11/#12 hold verbatim for every `an+b` map with `a` odd (`c = 2^L·a^p`, `diff = L−p`) — the map's own multiplier plays 3's role. For `5n+1` the special multipliers are 2 and 5; for `7n+1`, 2 and 7 | **Proved** (§15's proof never used a=3); verified on 5n+1 and 7n+1 ([§18](README_extra.md#18-the-whole-thing-was-never-about-3-generalization-to-every-anb-map)) |

## Still open

- No closed-form proof that Pr(diff_k = L) → 1 in general. (Even for the alternating family, only `m≡2 mod 4` is proved — see §17.) Note §15 now proves the *value* of diff is forced whenever a merge occurs — what remains open is entirely about the *frequency* of merging, not its value.
- No precise theoretical derivation of the two-component mixture's parameters (fast/slow split, the ≈bits^0.94 tail exponent) — currently fit, not derived.
- No proof that zero-exception early-merging holds for *all* n (only verified exhaustively/by sampling to large but finite scales).
- ~~The −p mechanism for powers of 3 (result #9) is empirically and partially (via merging classes) confirmed but has no closed-form proof analogous to §3's for base 2.~~ **Resolved by §15** for the *value* of the diff (it is forced to be −p); the frequency question remains open.
- ~~Whether other an+b Collatz-like maps (not just the standard 3n+1) show analogous effects — untested.~~ **Resolved by §18:** the `diff = L−p` law and the "multiplier must be `2^L·a^p`" constraint hold verbatim for every `an+b` map with `a` odd — §15's proof never used `a=3`. Verified for 5n+1 and 7n+1. (Existence/frequency of merging remains open, as always.)
- The true asymptotic growth of the power-of-3 merging-class fraction (result #10) is unknown — it doesn't match the classical case's near-logarithmic pattern in the tested range, but whether it's a genuinely different growth law or an eventual-flattening transient is unresolved. §16c reduces this to counting parity-vector pairs satisfying `2^L·γ = 3^i'·x + γ'`; note the classical analogue of this classification problem is *known* to be hard — Garner's proposed classification was disproven by Elia & Tucker (2015).
- For the *classical* merging-class fraction, §19 shows the log-growth law **must** break down between k≈31,000 and k≈91,000 (combining #6 and #5b's data via `agreement(B) ≥ f(B)`) — but where it bends, and to what, is unknown.
- **Verification task, not a research question:** obtain Gao (1993), *Discrete Math.* 112, and check whether its density table already establishes result #6's increasing trend (see §19). This gates whether #6 is novel.

## Repository contents

```
README.md                    this summary
README_extra.md              full write-up (19 sections + Addendum)
LICENSE                      MIT
results/                     raw CSV data from §2, §1's single-block run,
                              and §14's power-of-3 merging-class sweep

Core (§1-4):
collatz_block_repeat.py, big_survey.py, exact_power_test.py

§5a-5d / Addendum (n, n+1 agreement statistics):
random_m_test.py, merging_residue_classes.py, merging_classes_fast.c,
general_merging_test.py, large_scale_sampling.c, coupling_experiment3/4/5.py*,
coupling_scaling.c*, fit_gamma.py*, fast_slow_test.py*, minimal_K_test.py*,
hypothesis2_test.py*

§10-19 follow-ups (this session's independent re-verification + new findings):
k_cluster_analysis.py, quantile_scaling_analysis.py, quantile_scaling_large.py,
coupling_exhaustive_verify.py, base_generalization_test.py, power_of_3_test.py,
power_of_3_merging_classes.py, power_of_3_merging_classes_fast.py,
merging_classes_pow3_fast.c, merge_diff_theorem_verify.py,
alternating_family_theorem_verify.py, general_an_b_merge_verify.py

Unrelated side-quest:
cycle_search.c               exhaustive Collatz-cycle search, negative result

* starred files are referenced in README_extra.md but were lost (never saved
  to disk) before this repo could capture them — see README_extra.md §10 for
  the provenance note and independent re-verification of their claims.
```

## Relation to existing theory

The base-2 construction sits inside the 2-adic extension of the Collatz map studied by Bernstein & Lagarias ("The 3x+1 conjugacy map," *Canad. J. Math.* 1996): as k→∞, `N_k` converges 2-adically to `-x/(2^L-1)`.

**A literature survey (2026-07-26) found that result #5 above is not new** — it's a ~40-year-old known result (Garner 1985, formalized by LaDue 2017), and the general merging-class methodology of §5b/§5c parallels other published work (LaDue 2017, Winkler arXiv:1709.03385). It also found citation errors in an earlier, AI-assisted pass at this section (a fact wrongly attributed to Bernstein & Lagarias 1996 that actually belongs to Terras 1976). Full citation trail, corrections, and caveats: [README_extra.md §6](README_extra.md#6-relation-to-existing-theory).

**A second survey (same date) placed the merging-class work more precisely** ([§16d](README_extra.md#16d-literature-placement-updating-6)): the parity-vector-pair method used in §15/§16 *is* Garner's method, specialized there to consecutive integers; Elia & Tucker (*INTEGERS* 15, 2015, arXiv:1511.09141) prove §5b's theorem with an identical proof, and — importantly for #10 — **disprove Garner's conjectured classification of merging pairs**, showing the underlying structure is known to be complicated. What appears to remain this repo's own is the extension from `c=1` (consecutive integers, where `diff=0`) to general 3-smooth multipliers `c = 2^L·3^p`, giving the `diff = L−p` family (#11) and the 3-smoothness constraint (#12) — recorded with the standing caveat that elementary results are often already known.

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
