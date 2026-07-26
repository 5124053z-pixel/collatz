"""
Original README §5's open question: "Does the phenomenon generalize to bases
other than 2^L (i.e. repeating a block in base b for other b), or to other
an+b Collatz-like maps?" -- never previously tested.

The §1-§4 mechanism relies on N_{k+1} = 2^L * N_k + x, which (because 2^L is
a power of the Collatz map's own base) preserves N_k's *bits* unchanged as
the low-order bits of N_{k+1}: N_{k+1} = N_k (mod 2^{L*k}). That's the whole
proof mechanism in §3. If we instead build N_{k+1} = c * N_k + x for a
multiplier c that is NOT a power of 2, that congruence-preservation argument
breaks completely -- N_k's *bits* are no longer a fixed prefix of N_{k+1}'s
bits (only its base-c "digit" is preserved, which doesn't align with the
Collatz map's binary structure at all).

This tests empirically whether the diff_k = steps(N_{k+1}) - steps(N_k)
concentration phenomenon survives when c is not a power of 2, and if so,
around what value (the natural guess being round(log2(c)), i.e. the number
of bits c multiplies in on average).
"""
import math
from collections import Counter


def steps(n, cap=2_000_000):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        s += 1
        if s > cap:
            return None
    return s


def run(c, x, K, label):
    assert 0 <= x < c
    N = x if x != 0 else 1  # avoid N_1=0 degenerate case
    diffs = []
    prev_steps = steps(N)
    for k in range(1, K):
        N = c * N + x
        s = steps(N)
        if s is None or prev_steps is None:
            break
        diffs.append(s - prev_steps)
        prev_steps = s
    counts = Counter(diffs)
    top = counts.most_common(5)
    expected_bits = math.log2(c)
    print(f"{label}: c={c} x={x}  log2(c)={expected_bits:.3f}")
    print(f"  diff_k values (k=1..{len(diffs)}): top5 = {top}")
    if diffs:
        mean_diff = sum(diffs) / len(diffs)
        print(f"  mean diff_k = {mean_diff:.2f}   mode = {top[0][0]}  "
              f"mode_freq = {100*top[0][1]/len(diffs):.1f}%")
    print()


def main():
    K = 60

    print("=== control: c = power of 2 (should reproduce the original §2 finding) ===")
    run(c=2**8, x=0b10101010, K=K, label="c=2^8 (L=8 bits, alternating x)")
    print()

    print("=== c = power of 10 (decimal repeating block, matches §5's 'other base' question) ===")
    run(c=10**3, x=123, K=K, label="c=10^3 (decimal, x=123)")
    run(c=10**3, x=999, K=K, label="c=10^3 (decimal, x=999)")
    run(c=10**4, x=1234, K=K, label="c=10^4 (decimal, x=1234)")
    print()

    print("=== c = arbitrary non-power-of-2, non-power-of-10 multiplier ===")
    run(c=3**5, x=100, K=K, label="c=3^5=243 (base-3 digit block)")
    run(c=1000003, x=500001, K=K, label="c=1000003 (prime, ~2^20)")
    run(c=12345, x=6000, K=K, label="c=12345 (arbitrary composite)")


if __name__ == "__main__":
    main()
