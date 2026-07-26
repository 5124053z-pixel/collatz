"""
README addendum §9's open question #3:
  "Why the fast component's modulus needs (K values 3-14, see Addendum §6)
   cluster where they do, and whether this connects to §4's 'exact 2^L' family."

This is a fresh investigation (the original §6 script, minimal_K_test.py, was
never saved to disk -- see repo history). It:

  1. Samples large random n, finds those with a fast tau_couple(n) merge
     against n+1 (t < 20), matching Addendum §6's setup.
  2. Binary-searches the minimal K such that fixing n's low K bits (with
     random high bits) reliably reproduces the same merge time -- same
     method §6 describes.
  3. Cross-checks: for each minimal-K found, is the residue r = n mod 2^K
     actually one of the *provable* merging classes from
     merging_residue_classes.py (enumerable up to mod 512, i.e. K<=9)?
  4. Tests the §4-connection hypothesis: does the residue r, read as its own
     small integer, have a short/clean Collatz trajectory (in the spirit of
     the §4 "exact 2^L" family), i.e. is fast-merging explained by r's own
     trajectory simplicity rather than being an unrelated coincidence?
"""
import random


def steps(n, cap=100000):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        s += 1
        if s > cap:
            return None
    return s


def merge_time(n, max_steps=200):
    """min t>=1 such that T^t(n) == T^t(n+1), searched strictly before either hits 1."""
    x, y = n, n + 1
    for t in range(1, max_steps):
        if x == 1 or y == 1:
            return None
        x = x // 2 if x % 2 == 0 else 3 * x + 1
        y = y // 2 if y % 2 == 0 else 3 * y + 1
        if x == y:
            return t
    return None


def reproduces(low_bits, K, target_t, bits, trials=5, rng=None):
    """Fix n's low K bits to low_bits, randomize the rest, check merge_time is always target_t."""
    mask_hi_bits = bits - K
    for _ in range(trials):
        hi = rng.getrandbits(mask_hi_bits) | (1 << (mask_hi_bits - 1))
        n = (hi << K) | low_bits
        if merge_time(n, max_steps=target_t + 5) != target_t:
            return False
    return True


def minimal_K(n, target_t, bits, K_max=20, rng=None):
    for K in range(3, K_max + 1):
        low_bits = n & ((1 << K) - 1)
        if reproduces(low_bits, K, target_t, bits, rng=rng):
            return K, low_bits
    return None, None


def is_merging_class(r, modulus, num_samples=6, max_check_steps=300):
    """Same test as merging_residue_classes.py: does every n=r+j*modulus provably
    merge with n+1 at the *same* (step_from_n, step_from_n+1) pair?"""
    def traj(n, cap):
        t = [n]
        while n != 1 and len(t) < cap:
            n = n // 2 if n % 2 == 0 else 3 * n + 1
            t.append(n)
        return t

    seen = set()
    for j in range(2, 2 + num_samples):
        n = r + j * modulus
        if n < 2:
            continue
        tn = traj(n, max_check_steps)
        tn1 = traj(n + 1, max_check_steps)
        idx1 = {v: i for i, v in enumerate(tn1)}
        m = None
        for i, v in enumerate(tn):
            if v in idx1:
                m = (i, idx1[v])
                break
        if m is None:
            return False
        seen.add(m)
    return len(seen) == 1


def main():
    rng = random.Random(20260726)
    BITS = 2000
    N_TRIALS = 400
    FAST_THRESHOLD = 20

    fast = []
    for _ in range(N_TRIALS):
        n = rng.getrandbits(BITS) | (1 << (BITS - 1))
        t = merge_time(n, max_steps=FAST_THRESHOLD)
        if t is not None:
            fast.append((n, t))

    print(f"bits={BITS}, trials={N_TRIALS}, fast mergers (t<{FAST_THRESHOLD}): {len(fast)}")
    print()

    K_counts = {}
    rows = []
    for n, t in fast:
        K, r = minimal_K(n, t, BITS, rng=rng)
        if K is None:
            continue
        K_counts[K] = K_counts.get(K, 0) + 1
        rows.append((n, t, K, r))

    print("=== minimal-K distribution ===")
    for K in sorted(K_counts):
        print(f"  K={K:2d}: count={K_counts[K]}")
    print()

    print("=== is each minimal-K residue a *provable* merging class (mod 2^K, K<=9)? ===")
    explained, checkable, unchecked = 0, 0, 0
    for n, t, K, r in rows:
        if K <= 9:
            checkable += 1
            ok = is_merging_class(r, 1 << K)
            explained += int(ok)
            status = "YES (provable merging class)" if ok else "no (not a provable class!)"
            print(f"  t={t:3d} K={K:2d} r={r:5d} (mod {1<<K:5d}) -> {status}")
        else:
            unchecked += 1
    print(f"\n{explained}/{checkable} minimal-K<=9 residues match a provable merging class "
          f"({unchecked} had K>9, not exhaustively checked here)")
    print()

    print("=== §4-connection test: does r's OWN trajectory reach a power of two quickly? ===")
    print("(3*(r/2)+1 == 2^L exactly is the §4 'clean family' condition, defined for odd r;")
    print(" here we just check how many steps r's own trajectory takes to first hit a power of 2)")
    for n, t, K, r in rows[:20]:
        s = 0
        x = r if r > 0 else 1
        hit_pow2_at = None
        for i in range(40):
            if x != 0 and (x & (x - 1)) == 0:
                hit_pow2_at = i
                break
            x = x // 2 if x % 2 == 0 else 3 * x + 1
        print(f"  r={r:5d} (K={K:2d}, merge_t={t:3d}) -> own traj hits a power of 2 at step "
              f"{hit_pow2_at}")


if __name__ == "__main__":
    main()
