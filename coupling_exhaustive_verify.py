"""
Independent re-verification of Addendum §2's central exhaustive claim.

The original coupling_experiment3/4/5.py scripts that produced §2's table
(N=1e6 -> 477,245 agreeing pairs, 100.0000% early-merge, 0 exceptions;
similarly for N=5e6, 15e6, 60e6) were never saved to disk (see §10's
provenance note) -- so as it stood, that specific claim had *no*
surviving verification code at all, unlike §5b/§5c/§5d whose files were
recovered. This reconstructs the check from scratch and re-runs it at
N=1,000,000 to see whether the reported numbers hold up independently.

Definition (matching the Addendum's corrected definition from §1):
  tau_couple(n) = min { t >= 1 : T^t(n) = T^t(n+1) }, searched strictly
  before either trajectory reaches 1 (t < steps(n), given steps(n)==steps(n+1)
  is a precondition for a valid "agreeing pair").
"""
import sys


def steps_all(N):
    """Memoized total stopping time for every n in [1, N], via iterative path caching."""
    cache = {1: 0}
    for start in range(2, N + 1):
        if start in cache:
            continue
        path = []
        n = start
        while n not in cache:
            path.append(n)
            n = n // 2 if n % 2 == 0 else 3 * n + 1
        s = cache[n]
        for v in reversed(path):
            s += 1
            cache[v] = s
    return cache


def early_merge_time(n, s):
    """Search t in [1, s) for T^t(n) == T^t(n+1). Returns t, or None if not found."""
    x, y = n, n + 1
    for t in range(1, s):
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        y = y // 2 if y % 2 == 0 else 3 * y + 1
        if x == y:
            return t
    return None


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    print(f"Computing steps(n) for n=1..{N} (memoized)...")
    cache = steps_all(N + 1)  # need steps(N+1) too, for the pair (N, N+1)
    # cache is a dict here; recursion via 3n+1 may also populate entries > N+1, which is fine

    print("Scanning for agreeing pairs and checking early-merge...")
    agreeing = 0
    exceptions = []
    tau_sum = 0
    for n in range(1, N):
        if cache[n] == cache[n + 1]:
            agreeing += 1
            t = early_merge_time(n, cache[n])
            if t is None:
                exceptions.append(n)
            else:
                tau_sum += t

    frac = 100.0 * (agreeing - len(exceptions)) / agreeing if agreeing else float("nan")
    print()
    print(f"N={N:,}")
    print(f"agreeing pairs: {agreeing:,}")
    print(f"early-merge fraction: {frac:.4f}%")
    print(f"exceptions: {len(exceptions)}" + (f" (first few: {exceptions[:10]})" if exceptions else ""))
    if agreeing - len(exceptions) > 0:
        print(f"mean early-merge time (of successes): {tau_sum/(agreeing-len(exceptions)):.2f}")
    print()
    print("Addendum §2 originally reported for N=1,000,000: agreeing=477,245, early-merge=100.0000%, exceptions=0")


if __name__ == "__main__":
    main()
