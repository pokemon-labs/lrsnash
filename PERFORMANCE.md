# Performance work on lrsnash (2026-11-23)

This document records what was changed to make lrsnash faster on small games,
and why each change was tried. The short version: the biggest wins came from
counting what the program was *doing*, which turned out to be far more work
than the question needs. The remaining wins came from noticing that every
division in lrs's integer pivot is exact. Numbers are per game, built with
`-O2`, on a Ryzen 7 5825U (Zen 3). "Legacy" means the current `src/*.cc`
libraries.

## Results

| case | legacy, best working type | new | speedup |
|---|---|---|---|
| 9x9 constant-sum (`solve(ConstantSumGame*)`) | 3,542 µs (int128) | 2.0 µs (int64, LP) | ~1,700x |
| 15x15 constant-sum | 2,272,304 µs (GMP) | 11.8 µs (auto, LP) | ~190,000x |
| 9x9 general-sum, exact legacy result | 3,781 µs (int128) | 618 µs (int64) | 6.1x |
| 15x15 general-sum, exact legacy result | 1,434,309 µs (GMP) | 139,297 µs (int128) | 10.3x |
| 9x9 general-sum, first equilibrium | – | 138 µs | 27x |
| 9x9 general-sum, Lemke-Howson (opt-in) | – | 3.3 µs | ~1,150x |

Per arithmetic type at 9x9:

| case | int64 | int128 | lrsmp | GMP |
|---|---|---|---|---|
| constant-sum, legacy → new LP | overflow on every game → 1.97 µs | 3,542 → 3.65 µs (970x) | 17,277 → 42.8 µs (404x) | 11,799 → 21.0 µs (563x) |
| general-sum Last, legacy → new | overflow on every game → 618 µs | 3,781 → 1,117 µs (3.4x) | 18,571 → 13,069 µs (1.4x) | 12,463 → 6,203 µs (2.0x) |

None of these changes trades away exactness. Every result is still computed
in exact integer arithmetic, and overflow is still reported.

## How changes were chosen

`perf` was not available (`perf_event_paranoid = 4`), so the evidence came
from four sources:

- **callgrind** instruction profiles, by function and by source line;
- **counters** compiled into a scratch copy of the code: P1 vertices, P2
  calls and their outcomes, pivots, copies, cache restores, `reverse()`
  calls;
- **`rdtsc`** wall-clock timing around the pivot kernel and the two halves
  of the solver;
- **microbenchmarks** for one question each, such as "which way of doing an
  exact 128/64 division is fastest on this CPU?".

The loop was the same each time:

1. profile;
2. find the biggest cost;
3. ask what the math actually requires there;
4. microbenchmark the idea if it is cheap to test;
5. implement it;
6. prove it does not change results;
7. measure;
8. profile again.

Step 6 was never skipped. Each change had to pass one of two gates. Either the
results are bit-identical to the reference path (same dictionaries, or the
same complete sequence of equilibria as legacy), or, for the two new
algorithms, every result is checked to be an exact Nash equilibrium in GMP
rationals.

## 1. Count the work first

The first profile said what you would expect. `pivot()` was 47% of
instructions, and libgcc's `__divti3` (128-bit division) another 17%. That
suggests making arithmetic faster. Before doing that, a counting build was run
to see *how much* arithmetic each game needed:

| per game | P1 vertices | equilibria | pivots | dict copies | `reverse()` calls |
|---|---|---|---|---|---|
| 5x5 constant-sum | 24 | 1.0 | 71 | 72 | 123 |
| 9x9 constant-sum | 482 | 1.0 | 1,198 | 1,447 | 4,363 |
| 9x9 general-sum | 502 | 8.5 | 1,247 | 1,507 | 4,546 |
| 15x15 constant-sum | 48,530 | 1.0 | 115,243 | 145,602 | 728,152 |

lrsnash enumerates *every* vertex of player 1's best-response polytope. For
each vertex it builds and tests a player-2 subproblem. For a constant-sum game
exactly one vertex produces an equilibrium. So the zero-sum path, the one
the user runs almost always, did about 500 times more work than the question
needed at 9x9, and about 50,000 times more at 15x15. No arithmetic trick can
compete with not doing that work. That set the order: first change the
algorithm for constant-sum games, then make the remaining arithmetic cheap.

## 2. Constant-sum games: one LP instead of an enumeration

**Why it is allowed.** By the minimax theorem, the Nash equilibria of a
constant-sum game are exactly the pairs (optimal x, optimal y). They are
interchangeable and share one value. So any optimal pair answers
the constant-sum `solve` correctly. CONTEXT.md had already recorded that "first
equilibrium" is fully correct for constant-sum games.

**What is solved** (`src/lrsnash/zerosum.hpp`). Shift the payoffs to
A' = A − min(A) + 1 ≥ 1 and solve `max 1ᵀy' s.t. A'y' ≤ 1, y' ≥ 0`. The shift
makes y' = 0 feasible, so no phase 1 is needed, and it makes the LP bounded.
At the optimum z*:

- y = y'/z*;
- x = (reduced costs of the slacks)/z*, because the dual solution sits in the
  objective row;
- value = 1/z* + min(A) − 1.

So both strategies come out of one dictionary.

**How.** The tableau is lrs's integer (fraction-free) dictionary, pivoted with
the same exact kernel as below. The entering rule is Dantzig's until the first
degenerate pivot, then Bland's. Without a degenerate pivot, cycling is
impossible, and Bland's rule terminates from any basis, so the method always
finishes. Leaving rule: minimum ratio, ties broken by smallest index.

**Result.** About 8 pivots at 9x9 and 14 at 15x15, instead of about 1,200 and
115,000. 9x9 takes 2.0 µs. The one semantic difference: in degenerate games
with several optimal strategies, the LP may return a different, equally
optimal pair than legacy's "last enumerated" one. Legacy-exact behaviour is
still available as `solve(&constant_sum_game, &out, Mode::Last)`.

## 3. Exact division by multiplying with an inverse

**Evidence.** After the algorithmic change, pivots still dominate everything
that is left, and the expensive part of a pivot is the division in
`A[i][j] = (A[i][j]·Ars − A[i][s]·A[r][j]) / det`.

**The math that makes it cheap.** In lrs's integer dictionary every entry is a
minor of the input matrix (Bareiss / Sylvester identity), so this division is
always exact. Exact division by an odd number e equals multiplication by e⁻¹
modulo 2^64, and that inverse takes four Newton steps to compute. Write
det = 2^k·e'. Then compute e'⁻¹ once per pivot. Each entry becomes
q = (a·R − b·C_j) >> k, with R and C_j precomputed by multiplying with e'⁻¹:
two wrapping multiplies, a subtract and a shift, with no division at all. For
int128 the same works modulo 2^128.

**Deciding with a microbenchmark.** Hardware `idiv` on Zen 3 is decent, so
the trick needed measuring first. Per exact 128/64 division:

| method | ns |
|---|---|
| libgcc `__int128 / int64` | 2.89 |
| inline `idivq` | 2.73 |
| modular inverse + range check | 0.97 |

**Keeping overflow detection exact.** Wrapping arithmetic fails silently: if
the true quotient does not fit, the answer is garbage. So the fast path is
used only when the quotient is *proven* to fit: |q| < 2^(63−k). The proof is a
magnitude bound,
`(max|row i|·|Ars| + |b|·max|row r|)/det ≤ 2^(62−k)`, evaluated in doubles
with a factor-2 margin. Rows that fail it use exact 128-bit arithmetic (or
exact 256-bit for int128). Those rows throw only if the true result does not
fit, so the overflow semantics are exactly as before.

**A second round, driven by the line-level profile.** The first version
computed each row's bound with a separate pass over the row. callgrind showed
that pass costing *more* than the multiplies (39M against 34M instructions),
because GCC does not vectorize it at `-O2`. The fix: keep one magnitude bound
per dictionary (`Dict::bound`). The kernel accumulates it for free inside the
main loop, and the next pivot checks it once instead of scanning every row.
Any direct write to a dictionary marks the bound unknown.

**Result.** Kernel-only, 2.7x faster than the reference pivot on the LP. The
kernel now runs at about 2 TSC cycles per entry. On the zero-sum path the
fast kernels gave dictionaries bit-identical to the reference upstream pivot,
on tens of thousands of games: sizes up to 30x30, payoffs 0..2 (very
degenerate), duplicate rows and columns, and ±10^6. They overflowed on
exactly the same games.

Plain int64 now handles 9x9 games with benchmark-like payoffs (0..81): it
overflowed on none of 10,000 LP solves and none of 1,000 enumerations.
Legacy's SAFE mode checks operand sizes conservatively and overflowed on
*every* such game. Much larger payoffs can still overflow int64; the caller
then picks a wider type.

## 4. General-sum games: same answers, less work

Here the output must stay identical to legacy: the "last equilibrium found",
in lrs's enumeration order. Player 2's warm-start dictionary carries state
from one P1 vertex to the next, so every decision must be reproduced exactly.
The changes therefore remove work around the decisions, never the decisions.

The header-only port was profiled after the kernels went in (13.7M
instructions per game, against 46M in legacy):

- **Dictionary copies were 27%.** They copied element by element. For int64
  and int128 a dictionary is one contiguous block, so a single `memcpy`
  replaces the loop.
- **An O(m²) search** in `lrs_getfirstbasis2`, which builds the inequality
  order, became an O(m) marking pass. Upstream also sorts the linearity list
  by calling a bubble pass m−1 times; that is just a sort of distinct row
  indices, so it became `std::sort`.
- **`reducearray` was skipped.** It computes a gcd over every vertex's output
  for all ~500 P1 vertices, but only the equilibrium that gets reported needs
  output values. The visitor now gets exact numerators and a common
  denominator, and converts to floating point once.

Then a counter was added for how each player-2 subproblem ends (9x9, per
game):

| outcome | count |
|---|---|
| P2 calls | 503.9 |
| `getabasis2` fails (face empty) | 0.46 |
| feasible with no pivot (equilibria) | 8.4 |
| **infeasible at `primalfeasible`'s first step** | **494.6** |
| needs `primalfeasible` pivots | 0.45 |

So 98% of calls end at the very first check of `primalfeasible`. Legacy pays
for each of those calls with:

- a full dictionary copy (P2orig → P2);
- the warm-restart pivots;
- a second full copy (P2 → P2orig);
- about 9 `removecobasicindex` calls, each O(m);
- and only then the check that fails.

The restructured `nash2_main` produces the same values with less work:

- The warm-restart pivots run **in place on P2orig**. They are logged, and if
  `getabasis2` fails (rare) they are undone. The integer pivot is an
  involution: pivoting again on the same entering/leaving pair restores A,
  det, B, C, Row and Col exactly. On success, P2orig is already the state
  legacy would have copied back.
- `primalfeasible`'s first step is evaluated **read-only** on P2orig, with the
  linearity columns skipped. Removing columns preserves the basis order and
  the order of the remaining cobasic columns, so the answer is the same.
- P2 is copied and reduced only when the face might be nonempty.

The undo path was exercised 1,959 times in the sequence tests, with 1,326
pivots undone, and the sequences still matched.

**Result.** 9x9 general-sum went from 3,781 µs (legacy) to 1,117 µs (first
port) to 836 µs (copies, bound, O(m)) to 618 µs. The complete equilibrium
sequence (every output, in order) matched legacy on about 14,000 games up to
15x15 for all four arithmetic types. The speedup grows with size: 4.9x at
9x9, 10x at 12x12 and 15x15. Legacy is forced into GMP at those sizes, while
exact overflow detection keeps the new code in int128.

**What limits it now** (`rdtsc` split, 9x9):

- P1's reverse search: 334 µs. About 4,500 `reverse()`/`lrs_ratio` tests,
  branchy code that the CPU mispredicts.
- P2: 284 µs, mostly about 700 warm-restart pivots.
- The pivot kernel itself: 41% of wall time.

The next steps here are branch-light ratio tests and SIMD. Neither changes
the algorithm.

## 5. Picking one equilibrium faster

For general-sum games the remaining cost is the enumeration itself. Two
options were measured:

- **Stop at the first equilibrium (`Mode::First`).** 138 µs against 599 µs
  for the full enumeration, only 4.4x. Instrumenting showed the first
  equilibrium appears on average about a fifth of the way through, so early
  exit cannot give more.
- **Lemke-Howson (`Mode::Any`, `src/lrsnash/lemke.hpp`).** Complementary
  pivoting follows one path from the artificial equilibrium (0,0) to a real
  one, usually in a handful of pivots. Implementation:
  - first version: both tableaus stored all variables, with column c holding
    the variable with label c, so "enter the complement" meant "the same
    column in the other tableau";
  - the full-tableau Bareiss pivot used the same fast exact kernel. Since
    the second round each polytope is an lrs dictionary, which does half the
    work for the same path (section 10.2);
  - degeneracy is handled by the lexicographic ratio test over the slack
    columns (= det·B⁻¹), which never ties.

  Result: 3.3 µs at 9x9 (about 9 pivots) and 12.5 µs at 15x15. For every
  possible dropped label on 2,750 test games, every result was an exact
  equilibrium with exactly the reported payoffs. It is **opt-in only**: it
  returns *an* equilibrium, generally not legacy's "last" one, and its
  worst-case path length is exponential (Savani–von Stengel), although random
  games are short.

Lemke-Howson was also tried on constant-sum games. It loses to the LP there:
4.2 µs against 2.0 µs at 9x9, and 39 µs against 18 µs at 15x15. Its paths are
longer on zero-sum games, so the LP stays the constant-sum method.

## 6. Precision: exact detection, explicit choice, optional promotion

Exact overflow detection matters more than raw kernel speed for the user's
sizes. With payoffs of about 0..100, it lets int64 handle 9x9 and int128
handle 15x15. Legacy SAFE needed int128 for 9x9 and GMP for 12x12 and 15x15.

The arithmetic type is the caller's choice. `Solver64`, `Solver128`,
`SolverMP` and `SolverGMP` each use exactly one type, and on overflow they
throw `overflow_error` (user preference). The `Solver` alias is an optional
chain: int64, then int128, then GMP. A scaling study showed the first version
of that chain *losing* at 15x15 and 20x20: 22 µs against 15.5 µs for int128
directly. Most games overflowed int64 late, after most of the work was done,
and then restarted.

Fixing that with **promotion** came from the same exactness argument. Before
each LP pivot, `pivot_safe` checks whether the current width can provably
hold the result. If not, the dictionary is converted losslessly to the next
width and the solve continues. The values are exact in both widths, so every
later decision, and the result, is the same as in a solve done entirely in
the wider type. That was confirmed on 750 games from 12x12 to 30x30 against
pure GMP. After promotion, 15x15 runs in 11.8 µs, *faster* than int128 alone,
because the early pivots run in int64. In the enumeration, int64 overflows
within the first ~14 µs at 12x12, so restarting wastes nothing and promotion
was not needed there.

## 7. Measured and set aside

- **Compiler flags.** `-O3 -march=native` gives 13% on the general-sum path.
  The user's build flags are unknown, so the code was tuned at plain `-O2`.
- **lrsmp digit capacity.** `MP<12>` against the upstream default `MP<112>`
  gains only about 17%. lrsmp's cost is in its digit loops, not memory, so
  the upstream default stays.
- **Skipping rows where the pivot column is zero.** Only 2% of those entries
  are zero.
- **Smaller player as P1, or reverse-order enumeration for Last mode.** Both
  could reduce work, but either can change which equilibrium Last returns in
  degenerate games, which would break the exact-legacy guarantee.
- **Lemke-Howson for constant-sum games.** Slower than the LP (section 5).

## 8. Remaining opportunities

The branch-light reverse search and the AVX2 kernel from the first version
of this list were done in the second round (section 10). Still open:

- **A narrow path for int128 dictionaries with small entries.**
  - While every entry fits in 64 bits, one pivot entry costs three 64-bit
    multiplies instead of six (or 128-bit storage with the I64 arithmetic).
  - This helps `Solver128`, and the int128 restarts of the auto chain for
    Lemke-Howson and the enumeration.
  - It needs a bound that says "fits in int64": `hmag` (the high word) does
    not.
- **A float-guided exact LP for constant-sum games.**
  - Find the optimal basis with a double-precision simplex, solve that basis
    exactly (one small Bareiss solve), and check optimality exactly. Fall
    back to the exact simplex if the check fails.
  - The 9x9 LP is 70% kernel over about 8 pivots, so this could roughly
    halve it.
  - It can return a different optimal pair than today's `Mode::Any` in
    degenerate games, which is allowed but is a behaviour change.
- **Padding LP/LH rows to a multiple of 4** so that the AVX2 loop needs no
  scalar epilogue. By instruction count this is about even at 9x9 (3 vector
  iterations against 2 plus 2 scalar entries), so it waits for timing.
- **A fixed 256/512-bit integer type** with the same modular trick, for
  30x30 and larger, where GMP dominates: the LP takes 1.1 ms at 30x30 and
  23 ms at 50x50.
- **Inline storage for solvers that are constructed per call**
  (CONTEXT.md Q6), and the "promotion" idea for the enumeration path.

## 9. Bug found along the way

Legacy returns a wrong game **value** on about 0.4% of 9x9 games with
payoffs 0..1000, for example 0.0054 instead of 0.4374. `mptoi` truncates the
payoff numerator (about payoff × denominator) to `long`, while the strategy
numerators still fit. The new code converts exact values to double directly.

## 10. Second round (2026-11-24)

The user was running their own benchmarks during most of this round, so
after a first set of native measurements the evidence was callgrind
instruction counts. These are deterministic and do not depend on machine
load, and they are backed by reading the generated code.

The final state, after the fixes in section 10.6, was then timed natively:
`-O3 -DNDEBUG`, pinned to one core, best of 5, the builds alternating. The
"plain kernel" build has `-DLRSNASH_NO_AVX2`. Times are in µs per game; the
percentages are against `18834e9`:

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

Run-to-run noise is about 1-2%, so nothing is measurably slower than
`18834e9` any more. The commit before this fix, `b37a781`, was still 4-5%
slower on 2x2 and 4x4 games (section 10.6).

Section 10.4 predicted that AVX2's real gain would exceed its
instruction-count reduction, because callgrind counts each vector
instruction as one. The timings do not bear that out. The gains roughly
track the counts, against the plain kernel:

| case | time | instructions |
|---|---|---|
| 9x9 LP | -14% | -14% |
| 9x9 LH | -11% | -14% |
| 9x9 Last | -16% | -15% |
| 15x15 LH auto | -21% | -27% |

### 10.1 A regression that was codegen, not code

Timing `fff2737` against the session-5 refactor (`18834e9`, A/B
interleaved) showed the dictionary-pool fix working: 9x9 Last int64 -6%,
First -7%. But int128 Last was 6-10% *slower*, although the refactor was
bit-identical and `arith.hpp` had not changed.

Callgrind showed the kernel doing exactly the same work (the same 4,817,184
entry updates), but with one more instruction per entry. The disassembly
explained it:

- GCC had started inlining the rare exact 256-bit row (`pivot_row_exact`)
  into the kernel. The function doubled in size (541 to 1,130
  instructions), and the register pressure of the cold path spilled into
  the hot loop.
- Both versions reloaded the shift count from the stack on every entry and
  used the general 128-bit shift (`shrd`, `sar`, `test`, two `cmov`).

The fixes, all bit-identical:

- **Cold paths out of line.** Make both exact rows `noinline, cold`, so
  inlining choices elsewhere cannot change the hot loop.
- **Specialize on the shift k.** k = 0 (det odd, about half of all pivots)
  needs no shift at all. For int128 with 0 < k < 64 the shift is just `shrd`
  + `sar`. (int64 got the k = 0 case too at first; it was taken out again,
  see section 10.6.) Getting GCC to emit that took a trick:
  - inside the `k < 64` branch it knows the range, so it folds `k & 63` to
    `k`;
  - it then expands the general shift anyway, because the range is not
    visible at RTL expansion time;
  - an empty `asm volatile("" : "+r"(k))` inside the loop hides the range,
    and the `& 63` survives. Clang emits the short shift without help.
- **No libm call per pivot.** `ldexp` became `pow2`, which builds the power
  of two from its bits.

Result: int128 paths -13% to -17% in instructions. The per-entry loop is
30 instructions for k = 0 and 35 for 0 < k < 64, against 40 before.

### 10.2 Lemke-Howson without the identity columns

Lemke-Howson pivoted full tableaus: 1 + m + n columns per row. But n (or m)
of those columns belong to basic variables, and in a fraction-free tableau
such a column is just det times a unit vector. The same polytope as an lrs
dictionary keeps only the 1 + m cobasic columns.

With `D = -T` on the cobasic columns, the existing lrs pivot kernel
reproduces the tableau update exactly. Every pivot element is negative
(`T[r][s] > 0`), so the kernel negates column s, stores -det at (r, s) and
leaves row r alone, which is exactly the leaving variable's new column.

The lexicographic ratio test compares tableau columns of the initial
slacks. For a basic slack, that column is read implicitly as det or 0. Each
comparison is the same exact sign as before, so the path is identical.

This was proved by a new gate against the frozen tableau version: every
missing label, all four types, identical pivots, exact solutions and
overflow. The result is half the kernel work at 9x9: LH instructions -35%
(9x9 int64), -53% (15x15 int128), -46% (15x15 auto). The full-tableau
pivot then had no users and was deleted.

### 10.3 Branch-light reverse search

The callgrind branch simulator put about 50k mispredicts on a 9x9 Last
game. Leaving aside the kernel's short-loop exits, which a real predictor
handles better than the simulator, the biggest data-dependent source was
the candidate collection of the ratio test:
`if (negative(entry)) cand[n++] = j`. The sign of a dictionary entry is
effectively random.

Changes:

- Collect candidates without a branch (`cand[n] = j; n += negative(...)`).
- For int64/int128, replace upstream's sign shortcuts around `comprod` with
  `comprod` itself. The shortcuts only skip the product, and give the same
  sign, because both denominators are negative. Two multiplies cost less
  than up to four unpredictable branches.
- Let `resort` return the element's final position. This replaces two linear
  searches after every pivot.

Simulated mispredicts -18%, conditional branches -8%, instructions -1%.

### 10.4 AVX2 without intrinsics

The instruction counts below are from the first AVX2 version. The dispatch
later moved to one test at the top of `pivot`, for rows of 8 or more
entries; CONTEXT.md session 6 has the counts and timings of the final
version.

Compiled for AVX2, GCC already vectorizes the I64 fast row loop:

- each 64-bit product from three `vpmuludq` (32x32-bit);
- the arithmetic shift from a logical shift plus the sign;
- the magnitude bound from `vpcmpgtq`.

That is 28 instructions per 4 entries, against about 13 per entry scalar.

So the AVX2 kernel is the same C++ row loop instantiated a second time
inside a `target("avx2")` function, chosen per pivot by a flag that each
I64 object sets from `__builtin_cpu_supports` when it is constructed.

One detail mattered. The loop needed `__restrict` parameters, in a small
helper function: GCC ignores `__restrict` on block-scope pointers, and
without it every row paid a runtime aliasing check.

`test_kernel` checks both kernels against GMP, and every gate passes with
AVX2 on and with `-DLRSNASH_NO_AVX2`.

| instructions per game | `18834e9` | no AVX2 | AVX2 |
|---|---|---|---|
| 9x9 LP int64 | 20,082 | -1.2% | -9.4% |
| 9x9 LH int64 | 33,660 | -33.9% | -38.9% |
| 15x15 LP auto | 138,735 | -7.6% | -16.2% |
| 15x15 LH auto | 191,764 | -44.8% | -56.5% |
| 9x9 Last int64 | 5,500,058 | -1.8% | -11.6% |

The 9x9 gain is limited by the row length. Of each 10-entry row, 8 entries
go through the vector loop and 2 through the scalar epilogue, and every row
pays about 12 instructions of setup.

I expected the real speed-up to be larger than the instruction counts
suggested:

- the scalar loop is bound by its two `imul`s per entry, about 2 cycles on
  Zen 3;
- the vector loop needs 1.5 `vpmuludq` per entry, and Zen 3 issues two per
  cycle.

Timed afterwards, it roughly matches the counts (see the table at the start
of section 10). The AVX2 default stays on (CONTEXT.md Q8).

### 10.5 A kernel test that is known to catch bugs

All of this rests on the fixed-width kernels being exact. The old gates
only reached the paths that natural games reach. Before this round, no test
had produced a k >= 64 pivot for int128, and none had come near the
fast-path limit.

`tools/test_kernel.cc` checks the kernel contract directly against GMP:

- exact values;
- overflow if and only if a true result does not fit;
- a sound bound.

It runs on random Bareiss sequences, and on constructed inputs whose
results straddle the fast-path limit. With det e = 2^k e', a pivot element
of +-e*m and pivot-row entries that are multiples of e, every division is
exact, so any magnitude can be placed exactly where it is wanted.

It was mutation-tested: 6 of 6 injected bugs fail it. These include a fast
path loosened by a factor of 3.8. The random part of the test alone missed
that one, which is why the constructed cases exist.

### 10.6 The small-LP regression

The first timing of this round found one regression. The plain int64
kernel was 4-6.5% slower on constant-sum LPs from 4x4 to 9x9, although it
executed 1.2% fewer instructions. On AVX2 machines, dictionaries under 8
columns take that path too.

These are a few percent of sub-microsecond solves, so the investigation
first stopped after three cheap steps (A/B interleaved, best of 5):

1. **Drop the int64 k = 0 specialization.** It saved one `sar` per entry,
   but it added a second copy of the row loop and an unpredictable branch
   per pivot: k is 0 for about half of the pivots, in no pattern. Without
   it, the inner loop is instruction for instruction the one in `18834e9`.
   - 4x4 went from +6.0% to +3.5%, and 6x6 from +4.6% to +1.4%.
   - 9x9 did not move (+5.3%, then +5.2%).
   - The AVX2 path did not change.
2. **Restore `18834e9`'s scalar kernel verbatim** for the plain path. The
   function compiled to the same size, to the byte, and 9x9 was still
   +4.5%. So the rest of the gap is not in the kernel's code.
3. **Build both with the same code alignment**
   (`-falign-functions=64 -falign-loops=64`). The gap then closes at 4x4 and
   6x6 (0.0% and -0.1%) and shrinks to +2.3% at 9x9. Most of what remains
   is where the code lands in the binary.

Step 1 is kept, since it also removes code. Step 2 is not kept.

Step 2 left one lead. In the default build, its variant brought 4x4 and
6x6 LPs back level with `18834e9` (-0.5% and -0.1%, one run). In that
build `pivot` had become a small dispatcher, and GCC had inlined the whole
plain kernel into its callers. That was committed as `b37a781` with the
rest noted as a partial regression. The user then asked for two
follow-ups:

4. **Inline the plain kernel into its callers** (`always_inline` on
   `I64::pivot`; the AVX2 kernel stays out of line). Three variants were
   timed: GCC's choice, forced inline and forced out of line. Only forced
   inline is reliable, because GCC inlined the kernel in one build and not
   in the other. Against `b37a781`:
   - default build: 4x4 LP +4.0% -> +0.1%, 2x2 LH +5.2% -> +1.9%;
   - plain build: 4x4 LP +2.6% -> -2.0%, 9x9 LP +4.1% -> -1.3%;
   - Lemke-Howson at 9x9, Last and First did not change.

   Forced out of line, the kernel stayed 2-3% slower on small games. It is
   not known whether the cost was code placement or the call and frame
   setup around the kernel; inlining removes it either way.
5. **Time the AVX2 cut-over.** Two builds were compared at every row
   length: one always takes the AVX2 kernel, the other never does. AVX2 at
   each length, against the plain kernel:

   | row length (game) | LP | LH | Last |
   |---|---|---|---|
   | 5 (4x4) | +3.7% | +2.0% | |
   | 6 (5x5) | -0.4% | +1.2% | -3.2% |
   | 7 (6x6) | -4.5% | -0.9% | -7.7% |
   | 8 (7x7) | -8.6% | -8.4% | -12.3% |
   | 9 (8x8) | -10.7% | -9.3% | -13.0% |

   The cut-over moved from 8 to 7 entries (`kAvx2MinCols`), so 6x6 games
   now take AVX2 too. The 6x6 LP went from +2.0% to -6.2% against
   `18834e9`.

## Reproducing

See "How to re-verify quickly" in CONTEXT.md for the test commands, and
`tools/bench.cc` for the benchmark. The legacy benchmark numbers come from the
same game generator (`mt19937(seed)`, payoffs uniform in [0, den+1]) linked
against the CMake-built legacy libraries.
