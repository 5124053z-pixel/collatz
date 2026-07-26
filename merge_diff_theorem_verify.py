"""§15: computational verification of the coefficient-matching theorem.

Theorem (see README_extra.md §15). On a merging class for the pairing
(m, c*m+x) with c = 2^L * 3^p, the merge is forced to satisfy

    i = i' + p                   (odd-step counts)
    alpha - i = beta - i' - L    (halving counts)

hence  diff = beta - alpha = L - p,  and no other value is possible.

This explains, in one argument, three separate empirical findings of this
repo: diff = +L for the base-2 pairing (result #1 / §5a), diff = -p for the
power-of-3 pairing (result #9 / §13), and diff = 0 for the mixed 6^p case.

This script verifies both identities directly, by extracting the actual
odd-step and halving counts from real merging classes found by brute force.

Run: python merge_diff_theorem_verify.py
"""
import sys
sys.path.insert(0, ".")


def trace(n, steps):
    """Return (values, n_odd_steps) after exactly `steps` Collatz steps."""
    vals = [n]
    odd = 0
    for _ in range(steps):
        if n % 2:
            odd += 1
            n = 3 * n + 1
        else:
            n = n // 2
        vals.append(n)
    return vals, odd


def first_merge(m, n, max_steps=400):
    seen_m, seen_n = {m: 0}, {n: 0}
    cur_m, cur_n = m, n
    for step in range(1, max_steps):
        if cur_m != 1:
            cur_m = cur_m // 2 if cur_m % 2 == 0 else 3 * cur_m + 1
            if cur_m in seen_n:
                return (step, seen_n[cur_m])
            seen_m.setdefault(cur_m, step)
        if cur_n != 1:
            cur_n = cur_n // 2 if cur_n % 2 == 0 else 3 * cur_n + 1
            if cur_n in seen_m:
                return (seen_m[cur_n], step)
            seen_n.setdefault(cur_n, step)
    return None


def check(L, p, x, k_max=12, num_samples=5, verbose_limit=4):
    """Find merging classes for c = 2^L*3^p and verify the two identities."""
    c = (1 << L) * (3 ** p)
    predicted = L - p
    print(f"--- c = 2^{L}*3^{p} = {c}, x = {x};  predicted diff = L-p = {predicted} ---")
    shown = 0
    total_classes = 0
    violations = 0
    for k in range(1, k_max + 1):
        modulus = 1 << k
        for r in range(modulus):
            info = None
            ok = True
            for j in range(2, 2 + num_samples):
                m = r + j * modulus
                res = first_merge(m, c * m + x)
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
            total_classes += 1
            # extract odd-step counts i (from m) and i' (from n)
            m = r + 2 * modulus
            n = c * m + x
            _, i = trace(m, alpha)
            _, i_prime = trace(n, beta)
            id1 = (i == i_prime + p)
            id2 = (alpha - i == beta - i_prime - L)
            id3 = (beta - alpha == predicted)
            if not (id1 and id2 and id3):
                violations += 1
                print(f"  !! VIOLATION k={k} r={r}: alpha={alpha} beta={beta} "
                      f"i={i} i'={i_prime}  id1={id1} id2={id2} id3={id3}")
            elif shown < verbose_limit:
                print(f"  mod 2^{k}, r={r:5d}: alpha={alpha:3d} (i={i}), "
                      f"beta={beta:3d} (i'={i_prime})  -> i-i'={i-i_prime}=p OK, "
                      f"halvings {alpha-i}={beta-i_prime-L} OK, diff={beta-alpha} OK")
                shown += 1
        if total_classes >= 40:
            break
    print(f"  => {total_classes} merging classes checked, {violations} violations\n")


if __name__ == "__main__":
    # power-of-3 cases (§13/§14)
    check(L=0, p=2, x=1)
    check(L=0, p=3, x=1)
    # base-2 case (result #1 / §5a): predicted diff = +L
    check(L=1, p=0, x=1)
    check(L=3, p=0, x=5)
    # mixed 6^p case (§13: predicted diff = 0)
    check(L=2, p=2, x=1)
