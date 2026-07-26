"""
Extended version of quantile_scaling_analysis.py (§11), pushing the bit-length
range up to 65536 bits to match the original (lost) Addendum §5 run's scale,
for a tighter comparison of the fitted power-law exponents.
"""
import random
import math
import time


def tau_couple(n, max_t):
    x, y = n, n + 1
    for t in range(1, max_t + 1):
        if x == 1 or y == 1:
            return None
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        y = y // 2 if y % 2 == 0 else 3 * y + 1
        if x == y:
            return t
    return None


def quantile(sorted_vals, q):
    if not sorted_vals:
        return None
    idx = min(len(sorted_vals) - 1, int(q * len(sorted_vals)))
    return sorted_vals[idx]


def power_law_fit(xs, ys):
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
    windows = [
        (64, 800), (256, 600), (1024, 400), (4096, 300), (16384, 200), (65536, 120),
    ]

    per_window_taus = {}
    for bits, n_samples in windows:
        max_t = max(3000, bits * 10)
        t0 = time.time()
        taus = []
        agree = 0
        for _ in range(n_samples):
            n = rng.getrandbits(bits) | (1 << (bits - 1))
            t = tau_couple(n, max_t)
            if t is not None:
                agree += 1
                taus.append(t)
        taus.sort()
        per_window_taus[bits] = taus
        elapsed = time.time() - t0
        mean_tau = sum(taus) / len(taus) if taus else float("nan")
        print(f"bits={bits:6d}  n_samples={n_samples:5d}  agreed={agree:5d} "
              f"({100*agree/n_samples:.1f}%)  mean_tau={mean_tau:.1f}  ({elapsed:.1f}s)")

    print()
    print("=== quantiles of tau_couple (conditional on merging) ===")
    print(f"{'bits':>6} {'p10':>8} {'p50':>8} {'p90':>8} {'p99':>8} {'max':>8}")
    quantile_series = {"p10": [], "p50": [], "p90": [], "p99": []}
    valid_bits = []
    for bits, _ in windows:
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
    print("=== power-law exponent alpha (quantile ~ bits^alpha) ===")
    print("(original Addendum §5: p10~0.058, p50~0.384, p90~0.661, p99~0.935)")
    print("(this repo's earlier 64-4096 bit re-run, §11: p10~0.000, p50~0.487, p90~0.860, p99~1.026)")
    for label, ys in quantile_series.items():
        pairs = [(b, y) for b, y in zip(valid_bits, ys) if y and y > 0]
        if len(pairs) < 3:
            print(f"  {label}: not enough data points")
            continue
        xs = [p[0] for p in pairs]
        yv = [p[1] for p in pairs]
        alpha, r2 = power_law_fit(xs, yv)
        print(f"  {label}: alpha={alpha:.3f}  R^2={r2:.3f}  (n={len(pairs)} windows, {valid_bits[0]}-{valid_bits[-1]} bits)")


if __name__ == "__main__":
    main()
