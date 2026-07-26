"""§18: verification that §15/§16 generalize to every an+b Collatz-like map.

The §15 argument never uses a = 3. For the generalized map

    T(n) = n/2           if n is even
    T(n) = a*n + b       if n is odd        (a odd, b odd)

the affine form on a fixed parity vector is (a^i * m + gamma) / 2^h, so the
same coefficient matching gives, for the pairing (m, c*m + x):

    c = a^(i-i') * 2^((beta-i') - (alpha-i))

and hence (both exponents being forced non-negative, exactly as in §16a):

    uniform merging is possible ONLY for  c = 2^L * a^p,
    and then                              diff = beta - alpha = L - p.

So the "special multipliers" are always 2 and the map's own a -- for a = 3
this recovers §15/§16, but the structure is not special to 3.

This script tests a = 5 and a = 7 (and re-tests a = 3 as a control),
checking both identities on merging classes found by brute force.

Note: for a >= 5 the map has divergent orbits and several cycles, so
"total stopping time" need not be defined. The theorem is about the
relative step counts at the merge, which is well defined regardless;
that is what is checked here.

Run: python general_an_b_merge_verify.py
"""


def make_step(a, b):
    def step(n):
        return n // 2 if n % 2 == 0 else a * n + b
    return step


def first_merge(m, n, step, max_steps=250, cap=10 ** 40):
    """First shared value; returns (steps_from_m, steps_from_n) or None."""
    seen_m, seen_n = {m: 0}, {n: 0}
    cur_m, cur_n = m, n
    for s in range(1, max_steps):
        cur_m = step(cur_m)
        if cur_m > cap:
            return None
        if cur_m in seen_n:
            return (s, seen_n[cur_m])
        seen_m.setdefault(cur_m, s)
        cur_n = step(cur_n)
        if cur_n > cap:
            return None
        if cur_n in seen_m:
            return (seen_m[cur_n], s)
        seen_n.setdefault(cur_n, s)
    return None


def odd_steps(n, k, step):
    """Count odd-steps in the first k iterations."""
    c = 0
    for _ in range(k):
        if n % 2:
            c += 1
        n = step(n)
    return c


def check(a, b, L, p, x, k_max=11, num_samples=4, want=6):
    """Find merging classes for c = 2^L*a^p under the an+b map, verify
    i = i'+p  and  alpha-i = beta-i'-L,  hence diff = L-p."""
    step = make_step(a, b)
    c = (2 ** L) * (a ** p)
    predicted = L - p
    print(f"--- map {a}n+{b} | c = 2^{L}*{a}^{p} = {c}, x = {x} "
          f"| predicted diff = {predicted} ---")
    shown = checked = violations = 0
    for k in range(1, k_max + 1):
        modulus = 1 << k
        for r in range(modulus):
            info = None
            ok = True
            for j in range(2, 2 + num_samples):
                m = r + j * modulus
                res = first_merge(m, c * m + x, step)
                if res is None:
                    ok = False
                    break
                if info is None:
                    info = res
                elif res != info:
                    ok = False
                    break
            if not ok:
                continue
            alpha, beta = info
            checked += 1
            m = r + 2 * modulus
            i = odd_steps(m, alpha, step)
            i_prime = odd_steps(c * m + x, beta, step)
            id1 = (i == i_prime + p)
            id2 = (alpha - i == beta - i_prime - L)
            id3 = (beta - alpha == predicted)
            if not (id1 and id2 and id3):
                violations += 1
                print(f"  !! VIOLATION k={k} r={r}: alpha={alpha} beta={beta} "
                      f"i={i} i'={i_prime} id1={id1} id2={id2} id3={id3}")
            elif shown < 3:
                print(f"  mod 2^{k}, r={r:4d}: alpha={alpha:3d} (i={i}), "
                      f"beta={beta:3d} (i'={i_prime}) -> i-i'={i-i_prime}=p OK, "
                      f"halvings {alpha-i}={beta-i_prime-L} OK, diff={beta-alpha} OK")
                shown += 1
        if checked >= want:
            break
    print(f"  => {checked} merging classes checked, {violations} violations\n")


if __name__ == "__main__":
    # NOTE on triviality: taking x = b makes n = a*m + b = T(m) for odd m,
    # so the pair merges at (alpha, beta) = (1, 0) for a trivial reason --
    # the same degenerate case §13 flagged for p=1, x=1. The cases below use
    # x != b (or p >= 2) so the confirmations are genuinely non-trivial.

    # control: the familiar a=3 case
    check(a=3, b=1, L=0, p=1, x=2)

    # the point of this section: a != 3 behaves identically
    check(a=5, b=1, L=0, p=1, x=3, k_max=13)
    check(a=5, b=1, L=0, p=2, x=1, k_max=13)   # c = 25, predicted diff = -2
    check(a=5, b=1, L=1, p=1, x=1, k_max=13)   # c = 10, predicted diff =  0
    check(a=7, b=1, L=0, p=1, x=3, k_max=13)
