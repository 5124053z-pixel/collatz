"""§17: verification of the alternating-family theorem.

Theorem (see README_extra.md §17). Let L >= 2 be even and let
x = 2*(2^L - 1)/3 be the alternating block 1010...10 of bit length L --
i.e. §4's "100%-clean" family, characterized by 3*(x/2) + 1 = 2^L.
Then for EVERY m = 2 (mod 4) with m >= 6, the trajectories of m and
n = 2^L*m + x merge, and

    steps(n) - steps(m) = L    exactly.

(The sole edge case is m = 2, whose trajectory reaches 1 after a single
step and so never completes the 3 steps the argument needs -- exactly
analogous to the n = 4 exception in §5b's classical theorem.)

Explicitly, writing m = 2u with u odd and w = (3u+1)/2:
    m reaches w in 3 steps      (2u -> u -> 3u+1 -> w)
    n reaches w in L+3 steps    (n -> n/2 -> 2^(L+1)*w -> ... -> w)

This is the first proof that the §4 family actually yields diff = L;
§4 established only the algebraic identity 3*(x/2)+1 = 2^L and observed
the diff = L consequence empirically.

This script checks the two step-counts and the resulting diff directly.

Run: python alternating_family_theorem_verify.py
"""


def steps_to_one(n, cap=100000):
    """Total stopping time under the plain Collatz map."""
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        c += 1
        if c > cap:
            raise RuntimeError("no termination within cap")
    return c


def nth_iterate(n, k):
    for _ in range(k):
        n = n // 2 if n % 2 == 0 else 3 * n + 1
    return n


def check(L, num_u=400):
    """Verify the theorem for a given (even) L over many odd u."""
    assert L % 2 == 0 and L >= 2, "L must be even and >= 2 for x to be an integer"
    x = 2 * (2 ** L - 1) // 3
    assert 3 * (x // 2) + 1 == 2 ** L, "x must satisfy the §4 identity"

    bad_merge = bad_diff = 0
    for t in range(1, num_u + 1):   # t >= 1, i.e. u >= 3, m >= 6 (excludes the m=2 edge case)
        u = 2 * t + 1            # u odd
        m = 2 * u                # m = 2 (mod 4)
        n = (2 ** L) * m + x
        w = (3 * u + 1) // 2

        # claimed: m reaches w in 3 steps, n reaches w in L+3 steps
        if nth_iterate(m, 3) != w or nth_iterate(n, L + 3) != w:
            bad_merge += 1
            continue
        # claimed: the resulting total-stopping-time difference is exactly L
        if steps_to_one(n) - steps_to_one(m) != L:
            bad_diff += 1

    status = "ALL OK" if (bad_merge == 0 and bad_diff == 0) else \
             f"{bad_merge} merge failures, {bad_diff} diff failures"
    print(f"L={L:3d}  x={x:<22} m=2 (mod 4), m>=6, {num_u} values of u: {status}")


def check_edge_case(L=12):
    """The single exception: m=2 reaches 1 in one step, before the argument's
    3 steps can complete."""
    x = 2 * (2 ** L - 1) // 3
    m, n = 2, (2 ** L) * 2 + x
    print(f"    edge case m=2: steps(m)={steps_to_one(m)} (< 3), "
          f"diff={steps_to_one(n) - steps_to_one(m)} != L={L}, as expected")


def check_other_classes(L, num=200):
    """Control: the other residue classes mod 4 do NOT merge at the fixed
    (3, L+3) step pair -- i.e. the theorem's hypothesis is not vacuous."""
    x = 2 * (2 ** L - 1) // 3
    for r in (0, 1, 3):
        hits = 0
        for t in range(num):
            m = 4 * t + r
            if m < 1:
                continue
            n = (2 ** L) * m + x
            if nth_iterate(n, L + 3) == nth_iterate(m, 3):
                hits += 1
        print(f"    control m={r} (mod 4): fixed-step merge in {hits}/{num} cases")


if __name__ == "__main__":
    for L in (2, 4, 6, 12, 20):
        check(L)
    print()
    print("Edge case (documented exception):")
    check_edge_case(12)
    print()
    print("Control (the theorem is specific to m = 2 mod 4):")
    check_other_classes(12)
