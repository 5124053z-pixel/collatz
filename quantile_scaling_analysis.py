"""
Independent re-verification of Addendum §5 (quantile scaling of tau_couple).

The original tooling behind §3-5 of the Addendum (coupling_scaling.c,
fit_gamma.py) was lost (see §10's provenance note). This reproduces the
core claim from scratch in pure Python: that tau_couple's quantiles scale
as power laws in bit length, with the exponent increasing smoothly from
~0 (p10, essentially bits-independent) to ~1 (p99, nearly linear) --
i.e. the distribution is a mixture of an O(1) fast bulk and an O(bits)
heavy tail.

Pure Python can't reach 65536 bits in reasonable time, but the claim is
about the *scaling exponent* across a range, which should be visible
already across a few orders of magnitude (64 to 4096 bits here).
"""
import random
import math


def tau_couple(n, max_t):
    x, y = n, n + 1
    for t in range(1, max_t + 1):
        if x == 1 or y == 1:
            return None  # reached 1 without merging -> not an agreeing pair (or censored)
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        y = y // 2 if y % 2 == 0 else 3 * y + 1
        if x == y:
            return t
    return None  # censored: didn't merge within max_t


def quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    idx = min(len(sorted_vals) - 1, int(q * len(sorted_vals)))
    return sorted_vals[idx]


def power_law_fit(xs, ys):
    """log-log linear regression, returns (exponent, R^2)."""
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    sxx = sum((v - mx) ** 2 for v in lx)
    sxy = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    slope = sxy / sxx
    intercept = my - slope * mx
    ss_res = sum((ly[i] - (slope * lx[i] + intercept)) ** 2 for i in range(n))
    ss_tot = sum((v - my) ** 2 for v in ly)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, r2


def main():
    rng = random.Random(20260726)
    bit_windows = [64, 128, 256, 512, 1024, 2048, 4096]
    samples_per_window = 400

    per_window_taus = {}
    for bits in bit_windows:
        max_t = max(500, bits * 8)  # generous cap; mean tau grows sub-linearly to linearly with bits
        taus = []
        agree, total = 0, 0
        for _ in range(samples_per_window):
            n = rng.getrandbits(bits) | (1 << (bits - 1))
            t = tau_couple(n, max_t)
            total += 1
            if t is not None:
                agree += 1
                taus.append(t)
        taus.sort()
        per_window_taus[bits] = taus
        print(f"bits={bits:5d}  agreed(merged)={agree:4d}/{total}  "
              f"({100*agree/total:.1f}%)  mean_tau={sum(taus)/len(taus):.1f}" if taus else
              f"bits={bits:5d}  agreed=0/{total}")

    print()
    print("=== quantiles of tau_couple (conditional on merging) per bit-length window ===")
    print(f"{'bits':>6} {'p10':>8} {'p50':>8} {'p90':>8} {'p99':>8} {'max':>8}")
    quantile_series = {"p10": [], "p50": [], "p90": [], "p99": []}
    valid_bits = []
    for bits in bit_windows:
        taus = per_window_taus[bits]
        if len(taus) < 20:
            print(f"{bits:6d}  (too few merged samples: {len(taus)})")
            continue
        p10 = quantile(taus, 0.10)
        p50 = quantile(taus, 0.50)
        p90 = quantile(taus, 0.90)
        p99 = quantile(taus, 0.99)
        print(f"{bits:6d} {p10:8d} {p50:8d} {p90:8d} {p99:8d} {taus[-1]:8d}")
        valid_bits.append(bits)
        quantile_series["p10"].append(p10)
        quantile_series["p50"].append(p50)
        quantile_series["p90"].append(p90)
        quantile_series["p99"].append(p99)

    print()
    print("=== power-law exponent alpha (quantile ~ bits^alpha), log-log regression ===")
    print(f"(Addendum §5 originally reported: p10~0.058, p50~0.384, p90~0.661, p99~0.935)")
    for label, ys in quantile_series.items():
        # guard against zero/negative values before log
        pairs = [(b, y) for b, y in zip(valid_bits, ys) if y and y > 0]
        if len(pairs) < 3:
            print(f"  {label}: not enough data points")
            continue
        xs = [p[0] for p in pairs]
        yv = [p[1] for p in pairs]
        alpha, r2 = power_law_fit(xs, yv)
        print(f"  {label}: alpha={alpha:.3f}  R^2={r2:.3f}  (n={len(pairs)} windows)")


if __name__ == "__main__":
    main()
