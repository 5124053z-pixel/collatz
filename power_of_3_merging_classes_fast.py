"""
Faster, parallelized extension of power_of_3_merging_classes.py.

§13 found (via power_of_3_merging_classes.py, mod up to 1024) that provable
merging classes for the pairing (m, 3^p*m+x) always give diff = -p exactly,
mirroring §5b's theorem for the classical (n, n+1) pairing. §5b later pushed
its modulus up to 2^27 and found the merging-class fraction does NOT converge
to a stable limit -- it keeps climbing, with increment(k) * k roughly constant
(a logarithmic-type growth, frac(k) ~ frac(k0) + C*ln(k/k0)).

This script re-runs §13's experiment at much larger moduli to check whether
the same open-ended, non-converging growth pattern holds for the power-of-3
pairing too, or whether it behaves differently (e.g. converges, or grows at a
different rate).

Two speed changes vs. the original:
  1. Trajectories for m and c*m+x are grown one step at a time in lockstep,
     checking for a shared value as soon as it appears, instead of building
     full max_steps-length trajectories first and searching afterward.
  2. Residues are checked in parallel across all CPU cores via multiprocessing.

Usage:
    python power_of_3_merging_classes_fast.py <p> <x> <k_min> <k_max> [num_samples] [max_steps]

Example (matches §13's p=2, x=1 case, pushed further):
    python power_of_3_merging_classes_fast.py 2 1 3 20
"""
import sys
import time
from multiprocessing import Pool, cpu_count


def collatz_step(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def find_merge_incremental(m, n, max_steps):
    """Grow both trajectories one step at a time; return (step_m, step_n) at
    first shared value, or None if no merge within max_steps."""
    if m == n:
        return (0, 0)
    seen_m = {m: 0}
    seen_n = {n: 0}
    if m in seen_n:
        return (0, 0)
    cur_m, cur_n = m, n
    for step in range(1, max_steps):
        if cur_m != 1:
            cur_m = collatz_step(cur_m)
            if cur_m in seen_n:
                return (step, seen_n[cur_m])
            seen_m.setdefault(cur_m, step)
        if cur_n != 1:
            cur_n = collatz_step(cur_n)
            if cur_n in seen_m:
                return (seen_m[cur_n], step)
            seen_n.setdefault(cur_n, step)
        if cur_m == 1 and cur_n == 1:
            return None
    return None


# c, x, num_samples, max_steps set once per worker process via initializer,
# to avoid re-pickling them for every residue.
_C = _X = _NUM_SAMPLES = _MAX_STEPS = None


def _init_worker(c, x, num_samples, max_steps):
    global _C, _X, _NUM_SAMPLES, _MAX_STEPS
    _C, _X, _NUM_SAMPLES, _MAX_STEPS = c, x, num_samples, max_steps


def _check_residue(r_and_modulus):
    r, modulus = r_and_modulus
    merge_info = None
    for j in range(2, 2 + _NUM_SAMPLES):
        m = r + j * modulus
        n = _C * m + _X
        res = find_merge_incremental(m, n, _MAX_STEPS)
        if res is None:
            return None
        if merge_info is None:
            merge_info = res
        elif res != merge_info:
            return None
    return (r, merge_info)


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    x = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    k_min = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    k_max = int(sys.argv[4]) if len(sys.argv) > 4 else 18
    num_samples = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    max_steps = int(sys.argv[6]) if len(sys.argv) > 6 else 300
    c = 3 ** p

    nproc = cpu_count()
    print(f"c=3^{p}={c}, x={x}, num_samples={num_samples}, max_steps={max_steps}, workers={nproc}")
    print(f"{'k':>3} {'modulus':>12} {'merging':>10} {'fraction':>10} {'diffs':>10} {'time(s)':>8}")
    sys.stdout.flush()

    with Pool(nproc, initializer=_init_worker, initargs=(c, x, num_samples, max_steps)) as pool:
        for k in range(k_min, k_max + 1):
            modulus = 1 << k
            t0 = time.time()
            diffs = set()
            count = 0
            tasks = ((r, modulus) for r in range(modulus))
            for res in pool.imap_unordered(_check_residue, tasks, chunksize=max(1, modulus // (nproc * 8))):
                if res is not None:
                    r, (a, b) = res
                    count += 1
                    diffs.add(b - a)
            elapsed = time.time() - t0
            frac = count / modulus
            print(f"{k:3d} {modulus:12d} {count:10d} {100*frac:9.4f}% {str(sorted(diffs)):>10} {elapsed:8.1f}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
