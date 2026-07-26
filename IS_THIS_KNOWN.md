# Draft: "is this known?" question

Not part of the investigation — a draft to post as a reference request.
See the bottom of this file for venue/tagging notes.

---

## Title

Reference request: on merging Collatz trajectories, the step-count difference is forced to be `L − p` when the multiplier is `2^L a^p`

---

## Body

Let $a$ and $b$ be odd, and let $T$ be the generalized Collatz map

$$T(n) = \begin{cases} n/2 & n \text{ even} \\ an+b & n \text{ odd.}\end{cases}$$

Fix integers $c \ge 1$ and $x$, and consider pairs $(m,\; n = cm + x)$.

Call such a pair **uniformly merging on a set $S$** if $S$ is infinite and there are fixed $\alpha, \beta$ and fixed parity vectors such that, for every $m \in S$, the trajectories of $m$ and of $n = cm+x$ first meet at a common value after exactly $\alpha$ steps from $m$ and $\beta$ steps from $n$, following those same parity vectors. (Concretely: $S$ is a residue class $\bmod 2^k$ on which the merge is "the same merge" for every representative — the situation studied for consecutive integers by Garner and others.)

**Claim.** If such an $S$ exists, then necessarily

$$c = 2^L a^p \quad \text{for some } L, p \ge 0, \qquad\text{and}\qquad \beta - \alpha = L - p.$$

In particular:

* no uniformly merging pair exists at all unless $c$ is of the form $2^L a^p$ (for $a=3$: unless $c$ is 3-smooth);
* when it does exist, the difference of total stopping times is exactly $L-p$, independent of $x$ and of everything else.

**Proof.** On a fixed parity vector the map is affine (the "affine expansion lemma", going back to Terras): after $\alpha$ steps of which $i$ are odd-steps, the value reached is $(a^i m + \gamma)/2^{\alpha-i}$ with $\gamma \in \mathbb{Z}$ depending only on the parity vector. Writing the common merge value $V$ from both sides, with $i$ odd-steps from $m$ and $i'$ from $n$:

$$V = \frac{a^i m + \gamma}{2^{\alpha - i}}, \qquad V = \frac{a^{i'} n + \gamma'}{2^{\beta - i'}} = \frac{a^{i'} c\, m + a^{i'}x + \gamma'}{2^{\beta - i'}}.$$

These are affine in $m$ and agree for infinitely many $m$, so the coefficients of $m$ agree:

$$\frac{a^i}{2^{\alpha-i}} = \frac{a^{i'} c}{2^{\beta-i'}} \quad\Longrightarrow\quad c = a^{\,i-i'}\, 2^{\,(\beta-i')-(\alpha-i)}.$$

Since $\gcd(a,2)=1$ and $c$ is a positive integer, both exponents are $\ge 0$ (a negative power of $a$ would force $a \mid 2^k$; a negative power of $2$ would equate an even number with an odd one). Setting $p = i-i'$ and $L = (\beta-i')-(\alpha-i)$ gives $c = 2^L a^p$ and, subtracting, $\beta - \alpha = L - p$. $\square$

**Question.** Is this statement in the literature, in this or an equivalent form? I am looking for a reference rather than a proof.

**What I have already checked.** The $c = 1$, $x = 1$ case is the classical "consecutive integers of equal height" situation, where the claim degenerates to $\beta = \alpha$, and this is exactly the setting of

* L. E. Garner, *On heights in the Collatz $3n+1$ problem*, Discrete Math. **55** (1985) 57–64,
* B. LaTourette's senior thesis (the "corresponding stems" / "block prefix" formalism), and
* M. Elia and A. Tucker, *Consecutive integers and the Collatz conjecture*, INTEGERS **15** (2015), arXiv:1511.09141,

whose Definition 2 (block prefix, $T_b(x)+1 = T_{b'}(x+1)$) is the $c=1, x=1$ case of the constant-term condition $2^L\gamma = a^{i'}x + \gamma'$ that accompanies the claim above. So the *method* is certainly Garner's. What I have not been able to locate is any statement allowing $c \neq 1$ — i.e. allowing $\alpha \neq \beta$ — and thereby producing the $\beta - \alpha = L-p$ law and the $c = 2^L a^p$ constraint.

The argument is elementary enough that I would expect it to be known or considered immediate, which is why I am asking for a pointer rather than claiming novelty.

---

## Venue and tagging notes (not part of the post)

**Where.** MathOverflow, tagged `reference-request`, `nt.number-theory`, `collatz-conjecture`.
Fallback: if it is closed as too elementary, repost to math.stackexchange with the same tags.

**Why this framing.**
- Leads with a precise statement, then a complete short proof — so nobody has to spend effort deciding whether it is true. The question is purely "is this known", which is a well-established use of `reference-request`.
- The "what I have already checked" paragraph is doing most of the work. Collatz attracts a very large volume of low-quality amateur submissions, so demonstrating familiarity with Garner / LaTourette / Elia–Tucker up front is what separates this from that pile.
- Explicitly disclaiming novelty ("I would expect it to be known") lowers the temperature and makes it easy for someone to just drop a citation.

**Before posting.** If a MathSciNet or zbMATH subscription is available (a university account normally has one), search there rather than relying on web search — that is the actual gate on the "is it known" question, and finding the reference yourself is cheaper than asking. Suggested searches: Collatz/3x+1 combined with *coalescence*, *equal heights*, *total stopping time*, *parity vector*.

**Do not** include the other results from this repository in the same post. One question, one claim.
