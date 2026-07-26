"""
A new finding while testing README §5's open question ("does the phenomenon
generalize to bases other than 2^L?"): for multiplier c = 3^p (NOT a power
of 2), diff = steps(3^p * m + x) - steps(m) concentrates on exactly -p,
regardless of x, with frequency that appears to increase with m's bit
length -- the same qualitative shape as the original §5a random_m_test.py
result, but for a completely different, previously untested family.

Mechanism (seen directly by trajectory inspection): the trajectory of
3^p*m+x merges into the trajectory of m itself (literal value coincidence,
same style of evidence as the Addendum's coupling theory), but reaches the
shared point exactly p steps sooner -- so once merged, steps(3^p*m+x) =
steps(m) - p exactly.
"""
import random
from collections import Counter


def steps(n, cap=2_000_000):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        s += 1
        if s > cap:
            return None
    return s


def test(p, x, bit_lengths, trials=200, rng=None):
    c = 3 ** p
    print(f"p={p} (c=3^{p}={c}), x={x}")
    for bits in bit_lengths:
        diffs = []
        for _ in range(trials):
            m = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
            sm = steps(m)
            sn = steps(c * m + x)
            if sm is not None and sn is not None:
                diffs.append(sn - sm)
        counts = Counter(diffs)
        exact = counts.get(-p, 0)
        print(f"  m bits={bits:5d}: freq(diff=-p) = {100*exact/len(diffs):5.1f}%   "
              f"top3={counts.most_common(3)}")
    print()


def main():
    rng = random.Random(20260726)
    bit_lengths = [16, 64, 256, 1024, 4096]

    for p in [1, 2, 3, 5]:
        for x in [0, 1, 3 ** p - 1]:
            test(p, x, bit_lengths, trials=150, rng=rng)


if __name__ == "__main__":
    main()
