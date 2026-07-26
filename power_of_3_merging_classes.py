"""
§13's diff=-p phenomenon (for c=3^p) tested at the residue-class level,
the same style of rigor as §5b's theorem for the classical n, n+1 case.

Instead of "does n merge with n+1", here we ask "does m merge with
c*m+x" for a fixed multiplier c=3^p and offset x. A residue class
r (mod 2^k) is a "merging class" if every representative m = r (mod 2^k)
provably merges with c*m+x at the same (step_from_m, step_from_n) pair
-- which forces diff = step_from_n - step_from_m to be that class's
fixed, exact value (empirically always -p, matching §13's finding).
"""
import sys


def collatz_traj(n, max_steps=400):
    traj = [n]
    while n != 1 and len(traj) < max_steps:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        traj.append(n)
    return traj


def is_merging_class(r, modulus, c, x, num_samples=6, max_check_steps=300):
    seen = set()
    for j in range(2, 2 + num_samples):
        m = r + j * modulus
        if m < 1:
            continue
        n = c * m + x
        tm = collatz_traj(m, max_check_steps)
        tn = collatz_traj(n, max_check_steps)
        idx_n = {v: i for i, v in enumerate(tn)}
        merge = None
        for i, v in enumerate(tm):
            if v in idx_n:
                merge = (i, idx_n[v])
                break
        if merge is None:
            return False, None
        seen.add(merge)
    if len(seen) == 1:
        return True, seen.pop()
    return False, None


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    x = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    k_max = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    c = 3 ** p

    print(f"c=3^{p}={c}, x={x}")
    print(f"{'mod':>8} {'merging':>10} {'total':>10} {'fraction':>10}")
    for k in range(3, k_max + 1):
        modulus = 2 ** k
        merging = []
        for r in range(modulus):
            ok, info = is_merging_class(r, modulus, c, x)
            if ok:
                merging.append((r, info))
        diffs = {info[1] - info[0] for _, info in merging}
        print(f"{modulus:8d} {len(merging):10d} {modulus:10d} "
              f"{100*len(merging)/modulus:9.1f}%   diffs found: {sorted(diffs)}")


if __name__ == "__main__":
    main()
