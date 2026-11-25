# Performance work on lrsnash

This document explains what made lrsnash faster, biggest gains first.
"Legacy" is the stripped-down lrsnash 7.2 library this repo started from (git
tag `legacy`). The new code is the header-only library in `src/lrsnash/`. The
work came in two rounds. Round 1 (2026-11-23) changed the algorithms and the
arithmetic. Round 2 (2026-11-24) tuned the code.

The biggest wins came from counting what the program was *doing*, which was
far more work than the question needs. When one equilibrium is enough, the new
code no longer enumerates them all (sections 1 and 2). Everything else makes
the remaining exact work cheaper:

- a wider integer type only when a result really does not fit (section 3);
- pivots without division (section 4);
- less work around the enumeration's decisions (section 5);
- tuned hot loops (section 6).

None of this gives up exactness. All arithmetic is still exact, and overflow is
still reported. Each change had to pass one of two checks:

- its results are bit-identical to the reference: the same dictionaries, or
  the same complete sequence of equilibria as legacy;
- or, for the two new algorithms, every result is verified as an exact Nash
  equilibrium in GMP rationals.

## Summary

Times are µs per game on a Ryzen 7 5825U (Zen 3), one thread. Each side uses
the fastest arithmetic type that does not overflow on these games. "auto" is
the `Solver` chain of section 3: int64, then int128, then GMP. The new code's
times use the AVX2 kernel. The appendix explains how the numbers were taken,
and also gives the times without AVX2.

| case | legacy | new | speedup |
|---|---|---|---|
| 9x9 constant-sum, one equilibrium (LP) | 3,542 (int128) | 1.76 (int64) | ~2,000x |
| 15x15 constant-sum, one equilibrium (LP) | 2,272,304 (GMP) | 10.3 (auto) | ~220,000x |
| 9x9 general-sum, one equilibrium (Lemke-Howson)\* | 3,781 (int128) | 2.10 (int64) | ~1,800x |
| 15x15 general-sum, one equilibrium (Lemke-Howson)\* | 1,434,309 (GMP) | 6.70 (auto) | ~210,000x |
| 9x9 general-sum, legacy's result (`Mode::Last`) | 3,781 (int128) | 459 (int64) | ~8x |
| 12x12 general-sum, legacy's result (`Mode::Last`) | 157,347 (GMP) | 16,300 (int128) | ~10x |

The LP and Lemke-Howson are what `solve` uses by default.

\* Legacy can only run the full enumeration. Lemke-Howson returns one exact
equilibrium, which is generally not the one legacy returns (section 2).

## 1. Constant-sum games: one LP instead of an enumeration

**Why legacy was slow.** lrsnash enumerates *every* vertex of player 1's
best-response polytope. For each vertex it builds and tests a player-2
subproblem. A counting build of legacy measured the work per game:

| per game | P1 vertices | equilibria | pivots | dictionary copies | `reverse()` calls |
|---|---|---|---|---|---|
| 5x5 constant-sum | 24 | 1.0 | 71 | 72 | 123 |
| 9x9 constant-sum | 482 | 1.0 | 1,198 | 1,447 | 4,363 |
| 9x9 general-sum | 502 | 8.5 | 1,247 | 1,507 | 4,546 |
| 15x15 constant-sum | 48,530 | 1.0 | 115,243 | 145,602 | 728,152 |

On a constant-sum game exactly one of those vertices produced an equilibrium.
So the constant-sum path, the one used most, did about 500 times more work
than needed at 9x9, and about 50,000 times more at 15x15. No arithmetic trick
can compete with not doing that work.

**Why one LP is enough.** By the minimax theorem, the Nash equilibria of a
constant-sum game are exactly the pairs (optimal x, optimal y). They are
interchangeable and share one value, so any optimal pair is a correct answer.

**What is solved** (`zerosum.hpp`). Shift the payoffs to
A' = A − min(A) + 1 ≥ 1 and solve `max 1ᵀy' s.t. A'y' ≤ 1, y' ≥ 0`. The shift
makes y' = 0 feasible, so no phase 1 is needed, and it makes the LP bounded.
At the optimum z*:

- y = y'/z*;
- x = (reduced costs of the slacks)/z*, because the dual solution sits in the
  objective row;
- value = 1/z* + min(A) − 1.

So both strategies come out of one dictionary. It is lrs's integer
(fraction-free) dictionary, pivoted with the exact kernel of section 4.

- **Entering rule:** Dantzig's until the first degenerate pivot, then Bland's.
  Cycling needs degenerate pivots, and Bland's rule terminates from any basis,
  so the method always finishes.
- **Leaving rule:** minimum ratio, with ties broken by smallest index.

**Result.** The LP takes about 8 pivots at 9x9 and 14 at 15x15, instead of
about 1,200 and 115,000. In round 1, 9x9 took 1.97 µs (int64) against legacy's
3,542 µs (int128). 15x15 took 11.8 µs (auto) against legacy's 2.27 s (GMP).
Round 2 brought these to 1.76 µs and 10.3 µs.

**Trade-off.** Some degenerate games have several optimal strategies. There
the LP may return a different optimal pair than the one legacy enumerates
last. Both are optimal, and the value is the same. `Mode::Last` still gives
legacy's exact result.

## 2. General-sum games: Lemke-Howson for one equilibrium

A general-sum game has no such LP, and its equilibria can differ in their
strategies and their values. So legacy's "last equilibrium found" is one
choice among several, not *the* answer. When any one exact equilibrium will
do, `solve` now uses Lemke-Howson by default (`Mode::Any`, `lemke.hpp`).
`Mode::Last` still gives legacy's result.

**How it works.** Complementary pivoting follows one path from the artificial
equilibrium (0, 0) to a real one:

1. Drop one label.
2. Pivot alternately in the two polytopes. The variable that leaves one
   polytope duplicates a label, and the variable with that label enters the
   other polytope.
3. Stop when the dropped label is picked up again.

Degenerate games are handled by the lexicographic ratio test over the columns
of the initial slacks (= det·B⁻¹). That test never ties, so the path is unique
and ends. A 9x9 game takes about 9 pivots, against about 1,250 for the
enumeration.

**Storing only the cobasic columns** (round 2). The first version pivoted full
tableaus, with 1 + m + n columns per row. But in a fraction-free tableau, the
column of a basic variable is just det times a unit vector. Now each polytope
is an lrs dictionary that keeps only the cobasic columns.

- With D = −T on those columns, the existing lrs pivot kernel computes the
  tableau update exactly. Every pivot element is negative here. So the kernel
  negates column s, stores −det at (r, s) and leaves row r alone. That is
  exactly the leaving variable's new column.
- The lexicographic test reads the tableau column of a basic slack
  implicitly, as det or 0. So every comparison has the same sign as before.

The path and the result are the same, but the kernel does half the work: 35%
fewer instructions at 9x9 (int64). With the kernel changes of section 6, 9x9
got 29% faster and 15x15 43-45% faster.

**Result.** In round 1, with full tableaus, Lemke-Howson took 3.3 µs at 9x9
and 12.5 µs at 15x15. Now 9x9 takes 2.10 µs (int64), against 459 µs for the
full enumeration in `Mode::Last`. 15x15 takes 6.70 µs (auto).

**Trade-offs.**

- Lemke-Howson returns *an* equilibrium, generally not legacy's last one.
- Its worst-case path length is exponential (Savani–von Stengel), although
  random games give short paths.
- On constant-sum games it takes about twice as long as the LP (4.2 µs
  against 2.0 µs at 9x9), so the LP stays the constant-sum method.

**Stopping at the first equilibrium (`Mode::First`).** The other way to get
one equilibrium is to stop the enumeration at its first one, in legacy's
order. That gave only 4.4x: 138 µs against 599 µs at 9x9 in round 1, and
106 µs now. On average the first equilibrium appears about a fifth of the way
through, so stopping early cannot gain more.

## 3. Exact overflow checks: narrower integer types suffice

Each integer type is much faster than the next wider one. At 9x9 in round 1:

| case | int64 | int128 | GMP | lrsmp |
|---|---|---|---|---|
| enumeration | 618 µs | 1,117 µs | 6,203 µs | 13,069 µs |
| LP | 1.97 µs | 3.65 µs | 21.0 µs | 42.8 µs |

So using the narrowest type that works matters about as much as a faster
kernel.

**Legacy gave up early.** Upstream's `SAFE` mode checks the operands of an
operation, not its result.

- A 64-bit `qpiv`, (a·b − c·d)/e, throws if any of its four operands exceeds
  2^31 − 1 in magnitude, even when the result fits easily.
- For int128 the operand limit is 2^63 − 1.

On the benchmark games, legacy int64 overflowed on *every* 9x9 game, and
legacy int128 on every 12x12 and 15x15 game. So 9x9 needed int128, and 12x12
and up needed GMP.

**The new check is exact.** I64 and I128 throw only if a true result does not
fit. The kernel uses its wrapping fast path only where the result is proven to
fit. Every other row is computed exactly in double-width arithmetic
(section 4). Whenever neither the old nor the new check throws, the results
are the same.

- int64 handles 9x9: with payoffs 0..81 it overflowed on none of 10,000 LP
  solves and none of 1,000 enumerations.
- int128 handles 15x15.

This is why the enumeration's speedup over legacy grows with game size. In
round 1 it was 6.1x at 9x9 (int128 → int64). At 12x12 and 15x15 it was about
10x (GMP → int128).

**Promotion instead of restart** (LP only). The caller picks the type.
`Solver64`, `Solver128`, `SolverMP` and `SolverGMP` each use one type and throw
on overflow. The `Solver` alias ("auto") is an optional chain: int64, then
int128, then GMP. Its first version restarted in the next type after an
overflow. It lost to plain int128 at 15x15 and 20x20: 22 µs against 15.5 µs.
Most games overflowed int64 late, after most of the work was done.

Now the LP checks before each pivot, with `pivot_safe`, whether the current
width can provably hold the result. If not, it converts the dictionary to the
next width without loss and continues there.

- The values are exact in both widths. So every later decision, and the
  result, is the same as in a solve done entirely in the wider type.
- 15x15 then took 11.8 µs. That is *faster* than int128 alone, because the
  early pivots run in int64.

The enumeration still restarts. At 12x12, int64 overflows within the first
~14 µs, so a restart wastes little.

## 4. Pivoting without division

Every lrs pivot updates each dictionary entry as
`A[i][j] = (A[i][j]·Ars − A[i][s]·A[r][j]) / det`. In a profile of legacy
(int128, 9x9 constant-sum), `pivot()` was 47% of all instructions. libgcc's
128-bit division `__divti3` was another 17%.

**Exact division is a multiplication.** In lrs's integer dictionary every
entry is a minor of the input matrix (Bareiss / Sylvester identity), so this
division is always exact. Exact division by an odd number e equals
multiplication by e⁻¹ modulo 2^64. Four Newton steps compute that inverse.

Write det = 2^k·e', with e' odd. Once per pivot, the kernel computes e'⁻¹, and
from it R = |Ars|·e'⁻¹ and C_j = A[r][j]·e'⁻¹. Each entry is then
q = (a·R − b·C_j) >> k. That is two wrapping multiplies, a subtraction and a
shift, with no division. int128 does the same modulo 2^128.

Hardware `idiv` is decent on Zen 3, so a microbenchmark decided whether this
was worth doing. Time per exact 128/64 division:

| method | ns |
|---|---|
| libgcc `__int128 / int64` | 2.89 |
| inline `idivq` | 2.73 |
| modular inverse + range check | 0.97 |

**Used only where it is provably right.** If the true quotient does not fit,
wrapping arithmetic gives a wrong result without any error. So the fast path
is used only when |q| < 2^(63−k) is proven, with a factor-2 margin, by the
magnitude bound `(max|row i|·|Ars| + |b|·max|row r|)/det ≤ 2^(62−k)`. The
bound is evaluated in doubles. Rows that fail it are computed exactly in
128-bit arithmetic (256-bit for int128). They throw only if the true result
does not fit.

**A bound that costs nothing.** The first version found each row's maximum
with a separate pass over the row. callgrind showed that pass costing more
than the multiplies, 39M against 34M instructions, because GCC does not
vectorize it at `-O2`. Now each dictionary carries one bound:

- it is the OR of the magnitudes of all its entries, which is within a factor
  of 2 of the maximum;
- the kernel accumulates it in its main loop, at one OR per entry;
- the next pivot checks it once instead of scanning every row;
- any direct write to a dictionary marks the bound as unknown.

**Result.** The kernel alone is 2.7x faster than the reference pivot on the
LP, at about 2 TSC cycles per entry. It was checked on tens of thousands of
games: sizes up to 30x30, payoffs 0..2 (very degenerate), duplicate rows and
columns, and payoffs of ±10^6. Its dictionaries were bit-identical to those
of the reference lrs pivot, and it overflowed on exactly the same games.

## 5. The full enumeration: the same output with less work

`Mode::Last` must return exactly legacy's answer: the last equilibrium in lrs's
enumeration order. `extreme_equilibria` must return legacy's full sequence.
Player 2's warm-start dictionary carries state from one P1 vertex to the next,
so every decision must be reproduced exactly. These changes remove work
*around* the decisions, never the decisions. The complete output sequence
matched legacy on about 14,000 games up to 15x15, for all four types.

**Most player-2 subproblems fail at once.** A counter recorded how each
player-2 call ends (9x9, per game):

| outcome | count |
|---|---|
| P2 calls | 503.9 |
| `getabasis2` fails (face empty) | 0.46 |
| feasible with no pivot (equilibria) | 8.4 |
| **infeasible at `primalfeasible`'s first step** | **494.6** |
| needs `primalfeasible` pivots | 0.45 |

So 98% of calls end at the very first check of `primalfeasible`. Before that
check, legacy pays for each call with:

- a full dictionary copy (P2orig → P2);
- the warm-restart pivots;
- a second full copy (P2 → P2orig);
- about 9 `removecobasicindex` calls, each O(m).

The new driver (`nash.hpp`) produces the same values without most of this:

- **The warm-restart pivots run in place on P2orig**, and they are logged. If
  `getabasis2` fails (rare), they are undone. The integer pivot is an
  involution: pivoting again on the same entering/leaving pair restores A,
  det, B, C, Row and Col exactly. On success, P2orig is already the state that
  legacy would have copied back.
- **`primalfeasible`'s first step is evaluated read-only on P2orig**, skipping
  the linearity columns. Removing columns keeps the basis order and the order
  of the remaining cobasic columns, so the answer is the same.
- **P2 is copied and reduced only when the face might be nonempty.**

**Other work removed.**

- **`memcpy` copies.** A profile of the port showed dictionary copies at 27%
  of instructions, made element by element. For int64 and int128 a dictionary
  is one contiguous block, so a single `memcpy` replaces the loop.
- **Linear-time bookkeeping.** An O(m²) search in `lrs_getfirstbasis2`, which
  builds the inequality order, became an O(m) marking pass. Upstream sorts the
  linearity list with repeated bubble passes; that became `std::sort`.
- **No gcd reduction.** `reducearray` computed a gcd over the output of every
  P1 vertex, about 500 per game. Only reported equilibria need output values.
  The visitor gets exact numerators and a common denominator, and converts to
  floating point once.
- **Dictionaries reused across solves.** A size check that always failed made
  every solve reallocate all its dictionaries. Now a repeat solve of the same
  shape allocates nothing. That made 9x9 int64 `Mode::Last` 6% faster and
  `Mode::First` 7% faster.
- **A branch-light reverse search** (round 2). callgrind's branch simulator
  counted about 50k mispredicts on a 9x9 `Mode::Last` game. The biggest
  data-dependent source was the ratio test's candidate collection,
  `if (negative(entry)) cand[n++] = j`: the sign of a dictionary entry is
  effectively random. Three changes:
  - Candidates are collected without a branch: `cand[n] = j; n += negative(...)`.
  - For int64 and int128, `comprod` replaces upstream's sign shortcuts around
    it. The shortcuts decide some comparisons from signs alone, and `comprod`
    gives the same answer, because both denominators are negative. Two
    multiplies cost less than up to four unpredictable branches.
  - `resort` returns the element's final position. This replaces two linear
    searches after every pivot.

  Simulated mispredicts fell 18%, conditional branches 8% and instructions 1%.
  The effect on time was not measured separately.

**Result.** At 9x9, `Mode::Last` took 3,781 µs in legacy (int128). Round 1
brought it to 618 µs (int64) and 1,117 µs (int128). After round 2 it takes
459 µs and 989 µs.

## 6. Tuning the kernel code (round 2)

### 6.1 An AVX2 kernel for int64

Downstream builds use CMake's Release flags: `-O3 -DNDEBUG`, with no `-march`.
So AVX2 code must be chosen at run time. Compiled for AVX2, GCC already
vectorizes the int64 row loop:

- each 64-bit product comes from three `vpmuludq` (32x32-bit multiplies);
- the arithmetic shift comes from a logical shift plus the sign.

That is about 7 instructions per entry, against 13-14 for the scalar loop. So
the AVX2 kernel is the same C++ code, compiled a second time inside a
`target("avx2")` function.

- Each I64 object sets a flag from `__builtin_cpu_supports("avx2")` when it is
  constructed.
- `pivot` tests that flag once per pivot. Rows of 7 or more entries take the
  AVX2 kernel.
- `-DLRSNASH_NO_AVX2` leaves the AVX2 kernel out.

One detail mattered. Without `__restrict`, every row paid a runtime aliasing
check. GCC ignores `__restrict` on block-scope pointers, so the row loop gets
`__restrict` parameters in a small helper function.

Against the plain kernel, AVX2 saves:

| case | time saved |
|---|---|
| 9x9 LP | 9% |
| 9x9 Lemke-Howson | 6% |
| 9x9 `Last` and `First` | 15% |
| 15x15 Lemke-Howson (auto) | 20% |

At 9x9 the row length limits the gain. Of each 10-entry row, 8 entries go
through the vector loop and 2 through the scalar epilogue.

The cut-over was timed by row length, with one build always using AVX2 and one
never using it. The table shows the AVX2 build's change in time:

| row length (game) | LP | LH | Last |
|---|---|---|---|
| 5 (4x4) | +3.7% | +2.0% | |
| 6 (5x5) | -0.4% | +1.2% | -3.2% |
| 7 (6x6) | -4.5% | -0.9% | -7.7% |
| 8 (7x7) | -8.6% | -8.4% | -12.3% |
| 9 (8x8) | -10.7% | -9.3% | -13.0% |

Hence the threshold of 7 entries (`kAvx2MinCols`).

### 6.2 int128: keep the cold path out of the hot loop

After a refactor, int128 `Mode::Last` became 6-10% slower. The refactor
changed neither the results nor `arith.hpp`. callgrind showed the kernel doing
the same work, with one more instruction per entry. GCC had started inlining
the rare exact 256-bit row into the kernel. That doubled the kernel's size
(541 to 1,130 instructions), and the cold path's register pressure spilled into
the hot loop. The fixes are all bit-identical:

- **Cold paths out of line.** Both exact rows are `noinline, cold`, so
  inlining choices elsewhere cannot change the hot loop.
- **A shift specialized on k.**
  - k = 0 (det odd, about half of all pivots) needs no shift.
  - For 0 < k < 64, the shift is `shrd` + `sar` instead of the general 128-bit
    shift (`test` plus two `cmov`).
  - GCC emits the general shift even in the k < 64 branch. An empty
    `asm volatile` inside the loop hides k's range from it.
- **No libm call per pivot.** `ldexp` became `pow2`, which builds the power of
  two from its bits.

The per-entry loop went from 40 instructions to 30 (k = 0) or 35
(0 < k < 64). int128 paths ran 13-17% fewer instructions. Time fell 17.5% for
9x9 `Last`, 18.5% for 12x12 `Last` and about 13% for the 15x15 LP.

### 6.3 Small games

The first version of round 2 made the plain int64 kernel 4-6.5% slower on LPs
from 4x4 to 9x9, although it ran fewer instructions. Three changes removed
that:

- **No k = 0 case for int64.** For int64 it saved only one `sar` per entry,
  unlike for int128. It also added a second copy of the row loop and an
  unpredictable branch per pivot: k is 0 for about half of the pivots, in no
  pattern.
- **The plain kernel is inlined** (`always_inline` on `I64::pivot`) into its
  callers: the LP, Lemke-Howson and the lrs pivot.
  - The AVX2 kernel stays out of line.
  - Left to itself, GCC inlined it in some builds and not in others.
  - Out of line, the kernel stayed 2-3% slower on small games.
- **The AVX2 threshold moved from 8 to 7 entries** (see the table above).

Now no case is measurably slower than before the round. The worst cases are
within the 1-2% noise: 2x2 Lemke-Howson at +1.8% and the 4x4 LP at +0.4%.

## 7. Not done

Tried and set aside:

- **`-march=native`.** It gives 13% on the general-sum path, but downstream
  builds use no `-march`. The runtime-dispatched AVX2 kernel is used instead.
- **A smaller lrsmp digit capacity.** `MP<12>` gains only about 17% over
  upstream's `MP<112>`. lrsmp's cost is in its digit loops, not in memory.
- **Skipping rows where the pivot column is zero.** Only 2% of those entries
  are zero.
- **The smaller player as P1, or reverse-order enumeration for `Mode::Last`.**
  Either could cut work, but either can change which equilibrium `Last`
  returns in degenerate games.
- **Stack storage for solvers built per call.** A fresh solver costs about 2k
  instructions, 10-17% of a 9x9 LP or Lemke-Howson solve. Heap storage was kept
  (CONTEXT.md, decision 13). Instead, keep one solver per thread and reuse it:
  repeat solves of the same size allocate nothing.

Still open, biggest lever first:

- **A float-guided exact LP for constant-sum games.**
  1. Find the optimal basis with a double-precision simplex.
  2. Solve that basis exactly (one small Bareiss solve).
  3. Check optimality exactly. If the check fails, fall back to the exact
     simplex.

  The 9x9 LP spends 70% of its time in the kernel over about 8 pivots, so this
  could roughly halve it. In degenerate games it can return a different
  optimal pair than `Mode::Any` does now. That is allowed, but it changes
  behaviour.
- **A narrow path for int128 dictionaries with small entries.**
  - While every entry fits in 64 bits, a pivot entry costs three 64-bit
    multiplies instead of six.
  - This helps `Solver128`, and the auto chain's int128 restarts for
    Lemke-Howson and the enumeration.
  - It needs a bound that says "fits in int64". `hmag` (the high word) does
    not give that.
- **A fixed 256/512-bit integer type** with the same modular trick. It would
  help at 30x30 and larger, where GMP dominates: the LP takes 1.1 ms at 30x30
  and 23 ms at 50x50.
- **LP and Lemke-Howson rows padded to a multiple of 4**, so that the AVX2
  loop needs no scalar epilogue. By instruction count this is about even at
  9x9, so it needs timing first.
- **Promotion for the enumeration** (section 3).

## Appendix: measurements

All times are µs per game on a Ryzen 7 5825U (Zen 3), single-threaded. Games
come from the generator in `tools/bench.cc`, which the legacy benchmark also
used: `mt19937(seed)`, payoffs uniform in 0..den+1. den is 80 for games up to
9x9 and 100 for larger ones. `perf` was not available
(`perf_event_paranoid = 4`), so the other evidence came from:

- callgrind instruction counts;
- counters compiled into scratch copies of the code;
- `rdtsc` timing;
- microbenchmarks.

**Round 1** (2026-11-23, `-O2`). Legacy is the CMake-built legacy libraries;
"ovf" means legacy overflowed on every game.

| case | legacy int64 | legacy int128 | legacy lrsmp | legacy GMP | new |
|---|---|---|---|---|---|
| 9x9 constant-sum, LP | ovf | 3,542 | 17,277 | 11,799 | int64 1.97, int128 3.65, lrsmp 42.8, GMP 21.0, auto 2.08 |
| 15x15 constant-sum, LP | ovf | ovf | 6,659,186 | 2,272,304 | int128 15.5, auto 11.8 |
| 9x9 general-sum, `Last` | ovf | 3,781 | 18,571 | 12,463 | int64 618, int128 1,117, lrsmp 13,069, GMP 6,203 |
| 9x9 general-sum, `First` | | | | | 138 |
| 9x9 general-sum, Lemke-Howson | | | | | 3.3 |
| 12x12 general-sum, `Last` | | ovf | | 157,347 | 16,020 |
| 15x15 general-sum, `Last` | | ovf | | 1,434,309 | int128 139,297 |
| 15x15 general-sum, Lemke-Howson | | | | | 12.5 |

**Round 2** (2026-11-24, `-O3 -DNDEBUG`, pinned to one core).

- The three builds alternated, and each time is the best of 5 runs.
- "Before" is `18834e9`, the code after round 1.
- "Plain kernel" is built with `-DLRSNASH_NO_AVX2`.
- Run-to-run noise is about 1-2%.

| case | before | plain kernel | AVX2 (default) |
|---|---|---|---|
| 4x4 / 6x6 LP int64 | 0.310 / 0.763 | -2.4% / -4.2% | +0.4% / -6.2% |
| 7x7 LP int64 | 1.05 | -0.5% | -6.8% |
| 9x9 LP int64 | 1.96 | -1.4% | -10.0% |
| 9x9 LP auto | 1.99 | +0.3% | -9.9% |
| 15x15 LP auto | 12.02 | -7.7% | -14.6% |
| 15x15 LP int128 | 16.86 | -13.3% | -12.5% |
| 2x2 / 4x4 LH int64 | 0.138 / 0.471 | +1.2% / -7.3% | +1.8% / -5.9% |
| 6x6 / 9x9 LH int64 | 1.14 / 2.97 | -16.5% / -24.5% | -15.2% / -29.3% |
| 15x15 LH auto | 12.11 | -31.1% | -44.7% |
| 15x15 LH int128 | 26.41 | -42.7% | -42.6% |
| 6x6 / 9x9 Last int64 | 42.3 / 608.4 | -8.9% / -11.6% | -17.3% / -24.6% |
| 9x9 Last int128 | 1,199 | -17.5% | -17.5% |
| 9x9 First int64 | 141.2 | -11.6% | -24.9% |
| 12x12 Last int128 | 19,984 | -18.8% | -18.5% |

**How the summary was computed.** No timing was rerun for it.

- Each "new" time in the summary is round 2's "before" time, adjusted by its
  AVX2 change.
- Legacy was timed at `-O2` and the new code at `-O3`. The two flag sets timed
  the same within noise (CONTEXT.md, decision 8), so the summary's speedups
  are approximate.
- Absolute times also moved between the two rounds' runs, most for 12x12
  `Last`. Compare times within one table.

**Reproducing.** CONTEXT.md has the test commands ("How to re-verify
quickly") and the round-2 timing commands (next step 3). `tools/bench.cc` is
the benchmark.
