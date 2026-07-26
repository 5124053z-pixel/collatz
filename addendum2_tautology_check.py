"""Audit of Addendum §2 / result #7: the "100% early merge, zero exceptions"
result is a tautology, not an empirical finding.

Addendum §1 documents two versions of a definitional trap that were caught
and fixed. This checks whether the third (current) definition escapes it.
It does not.

Claim: if steps(n) = steps(n+1) = s, then T^(s-1)(n) = T^(s-1)(n+1) = 2,
so a "strict early merge" at t = s-1 < s always exists.

Reason: steps(n) = s means T^s(n) = 1. The only positive integer mapping to
1 under T is 2 -- T(2) = 1, while the odd branch 3k+1 = 1 has no positive
solution. So T^(s-1)(n) = 2 is forced, and likewise for n+1. Neither has
reached 1 before step s, so t = s-1 lies in the searched range [1, s).

Hence the measured "early-merge fraction" had to be 100.0000% with zero
exceptions, for every N and every bit length, independently of any property
of the Collatz map beyond "1 has a unique positive preimage."

The same two lines give the converse, so the stated equivalence

    steps(n) = steps(n+1)  <=>  tau_couple(n) < steps(n)

is a theorem with a two-line proof, not a finding: if the trajectories
agree at any t < s they coincide from t onward and therefore finish
together; and if they finish together they agree at s-1 by the above.

What this does NOT invalidate: the *distribution* of tau_couple studied in
Addendum §3-§5. Those merges occur far earlier than s-1 (mean tau is a
small fraction of s), which is a genuine, non-forced observation. Only the
"100% / zero exceptions" headline is vacuous.

Run: python addendum2_tautology_check.py
"""


def steps_to_one(n):
    c = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        c += 1
    return c


def nth(n, k):
    for _ in range(k):
        n = n // 2 if n % 2 == 0 else 3 * n + 1
    return n


def audit(N=200000):
    agreeing = 0
    forced_at_s_minus_1 = 0
    tau_much_earlier = 0
    tau_sum = 0.0
    s_sum = 0.0

    for n in range(1, N + 1):
        s = steps_to_one(n)
        if steps_to_one(n + 1) != s or s < 1:
            continue
        agreeing += 1

        # (a) the forced merge: both must sit at 2 one step before the end
        if nth(n, s - 1) == 2 and nth(n + 1, s - 1) == 2:
            forced_at_s_minus_1 += 1

        # (b) the genuine, non-forced quantity: the FIRST merge time
        a, b = n, n + 1
        for t in range(1, s):
            a = a // 2 if a % 2 == 0 else 3 * a + 1
            b = b // 2 if b % 2 == 0 else 3 * b + 1
            if a == b:
                tau_sum += t
                s_sum += s
                if t < s - 1:
                    tau_much_earlier += 1
                break

    print(f"agreeing pairs (n <= {N}):            {agreeing}")
    print(f"  both at value 2 at step s-1:        {forced_at_s_minus_1}"
          f"  ({100.0 * forced_at_s_minus_1 / agreeing:.4f}%)")
    print(f"  => 'early merge' was FORCED in      "
          f"{100.0 * forced_at_s_minus_1 / agreeing:.4f}% of cases,")
    print(f"     so the reported 100% / 0 exceptions carries no information.")
    print()
    print(f"  first merge strictly before s-1:    {tau_much_earlier}"
          f"  ({100.0 * tau_much_earlier / agreeing:.4f}%)")
    print(f"  mean tau / mean s:                  "
          f"{tau_sum / s_sum:.4f}   <- this part is NOT forced")


if __name__ == "__main__":
    audit()
