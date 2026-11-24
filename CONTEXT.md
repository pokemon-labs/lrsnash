# lrsnash (stripped) — working context

This file holds the background, goals, constraints and plan for the ongoing
verification, optimization and C++ rewrite work on this repo. Read it before
making changes. `PERFORMANCE.md` explains the performance work in detail: what
was measured, why each change was tried, and what it gained.

**Branches.** All work since the initial stripped port lives on the `rewrite`
branch: this file, `PERFORMANCE.md`, `src/lrsnash/` and the new `tools/`.
`main` holds only the legacy code. Downstream users still build the legacy
libraries until `lib.h` is switched over (see "Decisions").

## What this repo is

A stripped-down copy of **lrsnash from lrslib 7.2** (David Avis), cut down to the
part needed to compute Nash equilibria of two-player games and exposed to C++
callers through a small API (`src/lib.h`, namespace `LRSNash`).

- Original upstream source: `/home/user/lrslib-072/` (reference only, do not edit).
  - Relevant files: `lrsnashlib.c/.h` (the Nash driver, derived from `nash.c` in
    lrslib-060 by Terje Lensberg, 2015), `lrslib.c/.h` (reverse-search core),
    `lrsmp.c/.h`, `lrslong.c/.h`, `lrsgmp.c/.h` (arithmetic back ends),
    `lrsnash.c` (CLI driver), `setupnash*.c`, `2nash.c`, `games/` (sample inputs),
    `man/man1/lrsnash.1`.
  - Upstream builds: `lrsnash` (gmp), `lrsnash1` (64-bit, `-DLRSLONG -DSAFE`),
    `lrsnash2` (128-bit, `-DLRSLONG -DSAFE -DB128`), plus an `-DMP` build.
- This repo:
  - `src/lib.h`, `src/lib.cc` — public API and the Nash driver (port of
    `lrsnashlib.c`).
  - `src/lrslib.{h,cc}` — stripped reverse-search core.
  - `src/lrsmp.*`, `src/lrslong.*`, `src/lrsgmp.*` — arithmetic back ends.
  - `test.cc` — random constant-sum games, checks value and exploitability.
  - `benchmark.cc` — 10000 random 9x9 constant-sum games, fixed seed.
  - `CMakeLists.txt` — builds four static libraries. The test and benchmark
    targets are commented out and are only for local testing; downstream code
    includes `src/lib.h`.
  - `src/lrsnash/` — the new header-only implementation (branch `rewrite`).
    `lib.h` does not use it yet.
    - `arith.hpp` — the four arithmetic policies (`I64`, `I128`, `MP<W>`,
      `GMP`) with the fast exact pivot kernels.
    - `lrs.hpp` — the lrs reverse-search core (`Dict`, `Problem`).
    - `nash.hpp` — the lrsnash driver (exact upstream equilibrium order).
    - `zerosum.hpp` — the exact LP for constant-sum games.
    - `lemke.hpp` — exact Lemke-Howson.
    - `solve.hpp` — the `SolverT<Ar...>` entry points (`solve` and
      `extreme_equilibria`, overloaded on `ConstantSumGame` and
      `GeneralSumGame`; decision 9).
  - `tools/` — benchmarks, exactness tests and differential tests (see
    "Tooling" below).

## The algorithm and the paper

**Primary paper:** D. Avis, G. D. Rosenberg, R. Savani, B. von Stengel,
*"Enumeration of Nash equilibria for two-player games"*, Economic Theory
42(1):9–37, 2010 (published online 2009; the lrsnash man page cites it as
"Economic Theory 42 (2009) 9-37").

lrsnash is the "lrs" method described there. For an m×n bimatrix game (A, B),
it uses the best-response polyhedra of both players, each with a payoff
variable. It enumerates the vertices of player 1's polyhedron with lexicographic
reverse search. For each vertex x it takes the labels that are tight at x and
imposes the complementary constraints as **linearities** (equalities) on player 2's
polyhedron. It then enumerates the vertices y of that face; every (x, y) found is
an extreme equilibrium. Player 2's dictionary is saved after the first
successful `getabasis2`, and later calls restart from it (`linindex`) instead of
building from scratch, which saves many pivots.

Supporting references:
- D. Avis, K. Fukuda, *"A pivoting algorithm for convex hulls and vertex
  enumeration of arrangements and polyhedra"*, Discrete Comput. Geom. 8 (1992)
  295–313. This is reverse search vertex enumeration.
- D. Avis, *"lrs: A revised implementation of the reverse search vertex
  enumeration algorithm"*, in Polytopes — Combinatorics and Computation (2000).
  This is the lrs implementation, lexicographic pivoting and the integer
  (fraction-free) dictionary.
- B. von Stengel, *"Computing equilibria for two-person games"*, Handbook of
  Game Theory vol. 3 (2002). This covers the best-response polytope background.

## Usage profile (drives optimization priorities)

- **Small matrices:** usually 9x9, possibly up to about 15x15, and never close
  to 100x100. So per-call overhead (allocation, setup, `lrs_init`, dictionary
  copies) probably dominates over asymptotic behavior. Fixed-size or stack
  storage and compile-time dimension bounds are reasonable options.
- **Almost always zero-sum or constant-sum** (`solve_fast`), but not always.
  MatrixUCB adds exploration rewards to the payoffs, which makes the game
  general-sum (`solve_full`). Both paths must stay correct. The zero-sum path
  is the one to optimize hardest.
- **Many independent solves run in parallel by the caller.**
- `src/lib.h` is the actual import that downstream programs use, so its API is
  the contract. The CMake test and benchmark targets exist only for local
  testing.

## Current public API (`src/lib.h`)

- `solve_fast(const FastInput*, FloatOneSumOutput*)` — **constant-sum** game.
  `data` is the row player's payoff (row-major, rows×cols). The column player's
  payoff is `den - data`. The output value is divided by `den`.
- `solve_full(const Input*, FloatOneSumOutput*)` — general bimatrix game.
  `row_data` and `col_data` are row-major payoff matrices. Added in commit e6f5246.
- Output: `row_strategy[rows]`, `col_strategy[cols]` and `value`, all as floats.
  Only **one** equilibrium is returned: the last one enumerated overwrites the
  earlier ones. All equilibria are still enumerated, so there is no early exit.

Note: despite the intent, the code is **not header-only yet**. It is `.cc` files
compiled into four static libraries (`lrslib`=MP, `lrslib64`, `lrslib128`,
`lrslibgmp`). The arithmetic type is chosen by preprocessor macros (`MP`,
`LRSLONG [+B128]`, `GMP`, plus `SAFE` for overflow-checked 64/128).

## Planned API direction

- Let the caller choose between **one equilibrium** and **all equilibria**.
  - For constant-sum games any equilibrium is fully correct. The value is
    unique and equilibria are interchangeable, so solving one exact LP (or
    stopping early) is a pure speedup.
  - For general-sum games (`solve_full`), equilibria are **not** unique in value
    or strategy. Returning only one (legacy keeps the last one) is a policy
    choice, not "the" answer. An "all" mode returns every extreme equilibrium
    (x, y) together with both players' payoffs.
- **Implemented on `rewrite`** in `src/lrsnash/solve.hpp` (not yet exposed
  through `lib.h`):
  - `SolverT<Ar...>`, with aliases `Solver64`, `Solver128`, `SolverMP`,
    `SolverGMP` (explicit type; overflow throws) and `Solver` (optional auto
    chain).
  - Two independent choices (decision 9): the game class is the input
    type, and the method says what is returned.

    |                   | one equilibrium (`Mode`)  | every extreme equilibrium |
    |---|---|---|
    | `ConstantSumGame` | `solve(&g, &out, mode)`   | `extreme_equilibria(&g)`  |
    | `GeneralSumGame`  | `solve(&g, &out, mode, k)` | `extreme_equilibria(&g)` |

    `Mode {Last, First, Any}`. `solve` writes a `FloatEquilibrium`;
    `extreme_equilibria` returns `std::vector<Equilibrium>`.

  The defaults are decided; see "Decisions".
- The caller chooses the arithmetic type. Dynamic switching stays opt-in (user
  preference, 2026-11-23).

## Hard constraints (user requirements)

1. **Stay a header library.** The end state is header-only, and new work must not
   move away from that.
2. **No hacky globals.** Remove reliance on global or `__thread` mutable state:
   `FirstTime`, `lrs_global_list`/`lrs_global_count`, `dict_count`/`dict_limit`/
   cache counters, `lrs_digits`/`lrs_record_digits`, `MAXD*` for 128-bit, and
   `lrs_init()` side effects. Move that state into explicit context or solver
   objects.
3. **Single threaded.** Do not add internal parallelism, because callers will run
   many solves in parallel. Each solve must therefore be **reentrant and
   thread-safe by construction**, with no shared mutable state, rather than made
   safe through `__thread`.
4. **Support all four arithmetic types:** 64-bit (`lrslong`), 128-bit
   (`lrslong` + `B128`), `lrsmp` (built-in multiprecision) and GMP. These need to
   coexist in one header library, so the natural direction is a template
   parameter or traits type instead of macros. Keep the `SAFE` overflow checking
   semantics for the fixed-width types.

## Work plan

1. **Identify the paper.** Done (see above).
2. **Correctness of the port.** Done: the static per-function diff and the
   differential runs vs the upstream binaries passed (see the results log).
3. **Correctness of the original.** Done: exact NE checks and support
   enumeration passed. Still not checked: completeness of the extreme
   equilibrium set on degenerate games, which needs an independent
   vertex-pair enumerator.
4. **Algorithmic (math) optimizations.** Done for this round:
   - the exact LP for constant-sum games;
   - Lemke-Howson as an opt-in single-equilibrium mode;
   - first-equilibrium early exit;
   - the P2 warm restart in place, with a read-only feasibility pre-check.

   Considered but not done: choosing the smaller player as P1, and
   reverse-order enumeration for Last mode. Both can change which equilibrium
   Last returns in degenerate games. Details are in `PERFORMANCE.md`.
5. **Programming optimizations.** Done for this round:
   - division-free exact pivot kernels with a tracked magnitude bound;
   - `memcpy` dictionaries and pooled storage;
   - no gcd reduction of vertex output;
   - O(m) linearity bookkeeping;
   - lossless precision promotion for the LP.

   Optional next steps (see "Next steps"):
   - a branch-light `lrs_ratio`/`reverse`;
   - an AVX2 kernel;
   - a fixed 256/512-bit type.
6. **C++ rewrite (next major step).** Turn the header-only port into the only
   implementation, written as idiomatic modern C++ (C++17). It must keep the
   hard constraints: header-only, no globals, reentrant, and all four
   arithmetic types.
   1. **Public API.** Switch `src/lib.h` to the header-only solvers. The
      structs keep their fields, under the new names `ConstantSumGame`,
      `GeneralSumGame` and `FloatEquilibrium` (aliases of the legacy structs
      until then; at the switch they become the definitions and the legacy
      names go, decision 9). Expose explicit arithmetic choice, `Mode`,
      `extreme_equilibria`, and a reusable per-thread solver object.
      Deferred until after 6.4 (decision 3); the type is chosen in code
      (decision 4).
   2. **Build.** Replace the four static libraries with one CMake INTERFACE
      target, `lrsnash` (decision 4). Turn `tools/test_*` and `bench` into
      CTest/optional targets.
   3. **Remove legacy.** Delete `src/lib.cc`, `src/lrslib.*`, `src/lrsmp.*`,
      `src/lrslong.*` and `src/lrsgmp.*` once the switch is done. Before that,
      capture any legacy-vs-new regression data still needed.
      `tools/instrument_legacy.py` needs the legacy sources, so archive a
      legacy build or tag `main` first.
   4. **Idiomatic internals.** Done for `lrs.hpp` and `nash.hpp` (session
      5, see the progress notes). They used to mirror upstream C line by
      line, with upstream names, `long*` out-parameters and magic row
      indices. Behaviour stays **bit-identical**. Gates for every further
      refactor:
      - `tools/test_refactor.cc` (against the frozen pre-refactor headers);
      - the full-sequence differential test (`tools/test_api.cc`, against the
        archived legacy build);
      - `tools/test_lp.cc` and `tools/test_lh.cc`;
      - a static review against the paper and upstream (user preference).
   5. **Errors and output.** Exceptions only (no printf, exit or abort);
      already true for the new code. Exactly `rows`/`cols` output entries.
      Exact (rational) output for `extreme_equilibria` if wanted (Q2).
   6. **Optional:** restore rational `num/den` payoffs (scale each player's
      matrix by the lcm of its denominators), and add compile-time dimension
      bounds or stack storage for tiny games.

## Results log

### Step 2: fidelity of the stripped port (done, 2026-11-23): PASSED
- **Static check:** each function was diffed against upstream after stripping
  comments and running clang-format. The scripts are in the session
  scratchpad. Every difference falls into one of these groups:
  - removed debug, verbose and print code;
  - removed options that nash never sets (hull, nonnegative, voronoi,
    maximize/minimize, lponly, givenstart, restart, extract, bound/prune,
    maxdepth/mindepth, maxoutput);
  - `exit()` replaced with `throw`.
  Upstream sets `Q->nash = TRUE`, so upstream `copy_dict` uses the elementwise
  copy. The stripped version always copies elementwise, so behavior is the same.
  Odd but harmless: `divint` in `lrsmp.cc` repeats the line
  `c[...] = qh;` three times.
- **The game matrix maps correctly.**
  - P1 rows are `-B^T` (fast: `data - den`, full: `-col_data`).
  - P2 rows are `-A`.
  - Integer `den = 1` in `lrs_set_row_constraint` is equivalent to upstream
    `lrs_set_row`.
- **Dynamic check:** 13,200 random games, compared against upstream binaries
  built with the same arithmetic.
  - Dimensions 1..9 with payoffs 0..80, and 1..15 with payoffs 0..100.
  - A third of the games use payoffs 0..2 to force degeneracy, and some include
    duplicate rows.
  - Covers both `solve_fast` and `solve_full` with MP, 64, 128 and GMP.
  - The stripped "last equilibrium" matched upstream's last printed equilibrium
    in **every** case: 0 mismatches.
  - Upstream crashes (stack smash) on long input file paths. Use short relative
    paths when driving it.
- **Overflow behavior matches upstream.** Without `-DMA`, upstream also just
  exits on overflow, and we throw instead. It is frequent, though:
  - 64-bit overflows on about 13% of 9x9 games with payoffs 0..80, and about
    30% at 15x15 with payoffs 0..100.
  - 128-bit overflows on about 7% at 15x15 with payoffs 0..100.
  - MP and GMP never overflowed in these tests.

### Step 3: correctness of upstream lrsnash (done, 2026-11-23): PASSED
- **Exact checker:** an independent Python `Fraction` checker
  (`verify.py`, in the scratchpad) ran on 1,600 random games against upstream
  GMP.
  - Sizes 1..6 x 1..6.
  - Half the games are degenerate (payoffs -3..3, some duplicate rows).
  - Half are generic (payoffs ±10^6).
- **Results:**
  - Every reported (x, y) is an exact NE.
  - The reported payoffs are exact.
  - No duplicate equilibria were reported, and no game had zero equilibria.
  - For generic games, the reported set **equals** the set found by support
    enumeration (800 games).
- **Not checked:** completeness of the extreme-equilibrium set for degenerate
  games. That needs an independent vertex-pair enumerator.

### Static review: code vs. paper (2026-11-23)
The user has already done extensive testing by running the program, so
correctness work should be **static analysis**: check that the code matches
the math of the paper. The randomized runs above are supporting evidence only.

**Paper's method (ARSvS 2010, the lrs-based method).** For an m×n game (A, B):
- P̄ = {(x,v) : x ≥ 0, 1ᵀx = 1, Bᵀx ≤ 1v}. Label i ∈ M if x_i = 0; label
  j ∈ N if (Bᵀx)_j = v.
- Q̄ = {(y,u) : Ay ≤ 1u, y ≥ 0, 1ᵀy = 1}. Label i ∈ M if (Ay)_i = u; label
  j ∈ N if y_j = 0.
- (x,y) is an extreme equilibrium iff x and y are vertices and
  L(x) ∪ L(y) = M ∪ N.
- Algorithm: enumerate the vertices x of P̄ with lrs. For each x, enumerate
  the vertices of the face of Q̄ where every constraint whose label is
  **not** in L(x) is tight.

**Mapping to the code.** lrs rows are `b + a·z ≥ 0`, with column 0 = b,
columns 1..k = strategy, last column = payoff variable.
- Q1 (P̄, `BuildRepP1Is1`):
  - rows 1..m are `x_i ≥ 0` (label i);
  - rows m+1..m+n are `v − (Bᵀx)_j ≥ 0` (label m+j), where
    fast uses B = den − A, so the coefficients are `data − den`;
  - row m+n+1 is `Σx = 1` (EQ).
- Q2 (Q̄, `BuildRepP1Is0`):
  - rows 1..m are `u − (Ay)_i ≥ 0` (label i);
  - rows m+1..m+n are `y_j ≥ 0` (label m+j);
  - row m+n+1 is `Σy = 1` (EQ).
- **The same row number means the same label in Q1 and Q2.** That is the
  invariant the coupling relies on. ✔
- Row 0 (`FillFirstRow`) is overwritten by `getfirstbasis`, so it is
  irrelevant. ✔
- The payoff variables v and u have no sign row. They are free decision
  variables and always basic, so **negative payoffs are fine**. Positive
  payoffs are needed only by the normalized-polytope form, which this code
  does not use. ✔

**Coupling (`nash2_main`).** For every basic slack in P1 with
`A[Row[i]][0] ≠ 0` (slack > 0, so the constraint is not tight and its label
is ∉ L(x)), the code maps the index through `Q1->inequality` to an input row,
and that row becomes a linearity of Q2. `Σy = 1` (row m+n+1) is then
appended.
- Cobasic slacks and degenerate basic slacks (value 0) are tight, so they
  are correctly treated as labels present at x.
- The tightness test depends only on x, not on which lexmin basis
  represents x. ✔
- The repeated `reorder` call (nlin − 1 rounds of forward and backward
  bubble passes) is a full sort. It is needed because `getfirstbasis`
  requires linearity in sorted order. ✔
- `j < Q1->linearity[0]` only guards against the Σx row. That row is removed
  from the dictionary, so it never appears among the basics. ✔

**Completeness and soundness.**
- The vertices of a face of Q̄ are vertices of Q̄.
- Any y with L(y) ⊇ (M∪N) \ L(x) lies in that face.
- So the output set is exactly the set of completely labeled vertex pairs,
  which is the extreme equilibria. ✔
- The face is bounded: some i ∈ M is missing from L(x), because x_i = 0 for
  all i is impossible. That makes (Ay)_i = u tight and pins u. ✔
- P̄ is pointed, so lrs's lexmin rule outputs each vertex once. Each
  (x, y) pair is therefore emitted once. ✔
- lrs reverse search is correct on degenerate input because of the
  lexicographic ratio test (Avis 2000). That is the basis for completeness on
  degenerate games.

**Warm restart (`getabasis2`, `!FirstTime` branch).** Starting from player
2's dictionary saved after the previous successful call:
- each basic slack that is now a linearity is pivoted out against a cobasic
  non-linearity column with a nonzero entry;
- `i--` re-examines the position after `reorder1` shifts entries left. This
  happens when the entering index is larger than the leaving index. ✔
- If no such column exists, the row is constant on the face:
  - a zero constant means the linearity is redundant, and the slack stays
    basic and identically 0. ✔
  - a nonzero constant means the face is empty, the function returns FALSE,
    and P2orig is left unchanged. ✔
- Primal infeasibility created by these pivots is repaired by
  `primalfeasible`. ✔
- `FirstTime` is cleared only after the first *successful* build. Until
  then, every call restarts from the original dictionary. ✔
- `nredundcol` and `lastdv` from the first call stay valid, because they
  depend only on the column rank, not on the linearities. ✔

**Latent issue, harmless.** In the `!FirstTime` branch, redundant
linearities stay in `Q2->linearity`. So `lrs_getfirstbasis2` compacts
`Q2->inequality` as if they had been removed, and the mapping becomes off
by one. Q2's `inequality` is never used after that, because player 2's
output uses only decision variables, so results are unaffected. Do not start
relying on `Q2->inequality` in the rewrite without fixing this.

**Stripped vs. upstream.** This was a function-by-function static diff (see
step 2). Every difference is removed I/O, a removed option that nash never
sets, or `exit` replaced by `throw`.

### Confirmed bugs in the current stripped code (fix in the rewrite)
- **Leak on overflow:** every exception thrown mid-solve leaks P1, P2orig, P2
  with its cache ring of up to 50 dictionaries, Q1, Q2, and the output
  vectors. Nothing is RAII. Legacy 64-bit overflows on every full-size 9x9
  game (payoffs 0..81); the 13% figure was for random dimensions 1..9.
- **Wrong `value` from `mptoi` truncation** (confirmed 2026-11-23). The
  payoff numerator (about payoff × denominator) overflows `long` while the
  strategy numerators still fit. About 0.4% of 9x9 games with payoffs
  0..1000 return a wrong value, for example 0.0054 instead of 0.4374; the
  strategies are correct. Fixed in the new code, which converts exact values
  to double.
- **Data races (formally UB):** plain globals are written by every call to
  `lrs_init`:
  - `lrs_digits` and `lrs_record_digits` (MP and GMP);
  - `MAXDm`, `MAXDl` and `MAXDa` (128-bit);
  - in lrsmp, `gcd` has `static` locals that are initialized lazily.

  These races are benign in practice (every thread writes the same values),
  but they are real for parallel callers. lrsmp's `lrs_record_digits` is a
  process-wide high-water mark that every thread both reads and writes.
- **Buffer overrun:** `lrs_nashoutput` writes `Q->n - 1 = rows + 1` (or
  `cols + 1`) entries. So `row_strategy[rows]` and `col_strategy[cols]` receive
  the payoff value. Callers must over-allocate, which `benchmark.cc` does with
  `+2`. The harness confirmed the write past the end in 100% of runs.
- **Truncating conversion:** `mptoi` truncates. For 128-bit it casts
  `__int128` to `long`. For MP it calls `notimpl()`, which calls `exit(1)`,
  when the value is above about 10^18. The library should never call `exit`.
  Several `exit(1)` calls remain in `lrsmp.cc`, `lrslong.cc` and `lrsgmp.cc`.
- **Overflow messages:** `lrs_overflow` and `lrs_exit` `printf` to stdout before
  they throw.

## Initial observations (to verify, not yet confirmed bugs)

- `MAXCOL = 1000` fixed stack arrays in the `Fill*` functions overflow if
  rows or cols is greater than about 998. There is no bounds check. This is
  irrelevant at our sizes, but the arrays are oversized (stack usage) and a
  bounds assert would still be cheap.
- `lrs_nashoutput` converts through `mptoi` to `long` and then to `float`. For
  MP and GMP this can overflow or truncate large numerators and denominators.
  It also bypasses the exactness that the bignum types provide. **Confirmed**:
  see "Wrong `value` from `mptoi` truncation" above.
- `solve_full` divides `value` by `den`, the same as `solve_fast`, but its
  constraint rows do not use `den`. Confirm what convention callers expect and
  which player's payoff `value` represents. It is currently read from player 2's
  polytope, so it is player 1's payoff.
- The returned (row, col) pair is "last found". Confirm that it is always a
  consistent equilibrium pair. The row strategy is written only when that P1
  vertex produced at least one P2 match, so the pairing looks consistent, but
  this should be proven.
- `lrs_set_row_constraint` in `lib.cc` is a local replacement for upstream
  `lrs_set_row`/`lrs_set_row_mp`. Its denominator handling is effectively a
  no-op, which is fine for integer input, but it should be checked against
  upstream.
- In `nash2_main`, `for (i = 1; i < nlinearity; i++) reorder(...)` calls
  `reorder` repeatedly. Check whether this matches upstream and whether it is
  needed.
- Error handling prints with `printf` and then throws `std::exception()`.
- The test and benchmark targets are commented out in CMake, and `test.cc`
  only exercises `solve_fast`.
- `lrsgmp.cc` and `lrsmp.cc` declare `lrs_digits`/`lrs_record_digits` as plain
  globals, while `lrslong.cc` uses `__thread`. This is inconsistent, and MP and
  GMP are not thread-safe as written.

## Globals audit of the current stripped code (2026-11-23)
The current code is **not** global-free. Most state was made `__thread`. That
is safe for one solve per thread, but it is still shared state, and it breaks
when an exception escapes mid-solve. Some state was missed entirely and is
still a plain process-wide global:
- `MAXDm`, `MAXDl`, `MAXDa` (`lrslong.cc:19`, 128-bit), rewritten by every
  `lrs_init()`.
- `lrs_digits` and `lrs_record_digits` in `lrsmp.cc:27-28`.
  `lrs_record_digits` is read and written by every normalize.
- `lrs_digits` and `lrs_record_digits` in `lrsgmp.cc:16-17` (written but never
  read).
- `static` locals in lrsmp `gcd()` (`lrsmp.cc:285-292`).

The user's callers allocate `rows+2` / `cols+2` output arrays, which is why the
write past the end never crashed. The new code must write exactly `rows` and
`cols` entries, so the existing `+2` buffers keep working.

## Performance work (2026-11-23, session 3): results

**Where the time went (legacy, 9x9 constant-sum).** P1 has about 482
vertices; each one triggers a P2 subproblem, and only one of them is the
equilibrium. That is about 1,200 pivots, 1,450 dictionary copies and 4,400
`reverse()` tests per game. 65% of instructions were the per-entry `qpiv`
including libgcc `__divti3`. At 15x15 there are about 48,500 P1 vertices per
game.

**Tricks implemented (all exact):**
1. **Constant-sum games: one LP** (`zerosum.hpp`) instead of the
   enumeration: about 8 pivots at 9x9 and 14 at 15x15. Both strategies come
   from one optimal dictionary: y from the primal, x from the reduced costs
   of the slacks. Payoffs are shifted to >= 1, so y' = 0 is feasible and the
   LP is bounded. Dantzig's rule until the first degenerate pivot, then
   Bland's rule, so it always terminates. It returns *an* optimal pair, which
   in degenerate games can differ from legacy's "last" equilibrium; the
   value is unique.
2. **Division-free exact pivots** (`arith.hpp`, `I64::pivot` /
   `I128::pivot`). Every lrs division is exact (Bareiss), so dividing by
   det = 2^k e' becomes multiplying by e'^-1 mod 2^64 (or 2^128) plus an
   arithmetic shift: two wrapping multiplies per entry. The fast path needs
   |q| < 2^(63-k). This is proven per pivot from a magnitude bound kept with
   each dictionary (`Dict::bound`, updated inside the kernel loop), or else
   per row. Rows that fail use exact 128-bit (I64) or 256-bit (I128)
   arithmetic, so overflow is still reported iff a true result does not fit.
   Microbenchmark (Zen 3): 0.97 ns per division, vs 2.9 ns (libgcc) and
   2.7 ns (`idiv`).
3. **Lemke-Howson** (`lemke.hpp`): one equilibrium of a general game. It is
   exact and lexicographic (safe on degenerate games), about 9 pivots at 9x9.
   **Opt-in only** (`Mode::Any`), because the equilibrium it finds differs
   from legacy's.
4. **Exact general-sum port** (`nash.hpp`). Decisions are identical to
   upstream, with less work per P1 vertex:
   - P2orig gets its warm-restart pivots in place, with an undo log for the
     rare `getabasis2` failure (an integer pivot is an involution).
   - `primalfeasible`'s first step is emulated read-only, skipping the
     linearity columns. 98% of P2 calls end there with an empty face. P2 is
     copied and reduced only when the face may be nonempty.
   - Vertex outputs are not gcd-reduced (`reducearray`).
   - Dictionaries are copied with `memcpy`.
   - The linearity ordering is O(m), and `std::sort` replaces the repeated
     bubble passes.
5. **Precision promotion** (only for the `Solver` auto alias; LP only).
   Before each pivot, `pivot_safe` checks whether the current width can
   provably hold the result. If not, the dictionary is converted losslessly
   to the next width and the solve continues, so nothing is redone. The
   enumeration restarts on overflow instead; there int64 fails within the
   first microseconds, so nothing is wasted.

**Verification (exact checks plus differential vs legacy):**
- The fast kernels give dictionaries bit-identical to the reference lrs
  pivot (`dense_pivot` / `dense_tableau_pivot`) and across
  I64/I128/lrsmp/GMP. This was checked on tens of thousands of games: sizes
  up to 30x30, degenerate payoffs 0..2, duplicate rows and columns,
  payoffs of ±10^6. The fast and reference kernels overflow on exactly the
  same games.
- LP: every result is an exact equilibrium (GMP rationals). The value equals
  legacy's on about 30k games, except 2 where legacy is wrong (see below).
- LH: for every possible dropped label, every result is an exact NE with
  exact payoffs (about 33k runs).
- `nash.hpp` / `solve_all`: the **complete output sequence** (every
  `lrs_nashoutput`, in order) is identical to the legacy library. This holds
  for all 4 types and auto, on about 14k games up to 15x15. The undo path was
  exercised: 1,959 undo events, 1,326 pivots undone.
- Promotion: auto equals pure GMP on 750 games, 12x12 to 30x30.

**Timings** (µs per game, `-O2`, Ryzen 7 5825U; legacy = the current
`src/*.cc` libraries; payoffs 0..81 at 9x9, 0..100 otherwise; "ovf" = legacy
overflows on every game):

| case | legacy i64 | legacy i128 | legacy MP | legacy GMP | new |
|---|---|---|---|---|---|
| 9x9 const-sum, LP | ovf | 3,542 | 17,277 | 11,799 | i64 1.97, i128 3.65, MP 42.8, GMP 21.0, auto 2.08 |
| 15x15 const-sum, LP | ovf | ovf | 6,659,186 | 2,272,304 | i128 15.5, auto 11.8 |
| 9x9 general, Last (exact legacy) | ovf | 3,781 | 18,571 | 12,463 | i64 618, i128 1,117, MP 13,069, GMP 6,203 |
| 9x9 general, First | | | | | 138 |
| 9x9 general, Any (LH) | | | | | 3.3 |
| 12x12 general, Last | | ovf | | 157,347 | 16,020 |
| 15x15 general, Last | | ovf | | 1,434,309 | 139,297 |
| 15x15 general, Any (LH) | | | | | 12.5 |

Scaling: the enumeration speedup grows with size (4.9x at 9x9, 10x at 12x12
and 15x15), because legacy is forced into GMP while the new code stays in
exact int128. The LP at 30x30 or larger is GMP-bound: 1.1 ms (30x30) and
23 ms (50x50).

**New confirmed legacy bug: wrong `value`.** `mptoi` truncates the payoff
numerator to `long`. On 9x9 games with payoffs 0..1000, about 0.4% return a
wrong value, for example 0.0054 instead of 0.4374. The payoff numerator
(about payoff × denominator) overflows while the strategy numerators still
fit. The new code converts exact values to double directly.

**Correction:** legacy `lrslib64` overflows on 100% of full-size 9x9 games
(payoffs 0..81), and `lrslib128` on 100% of 15x15 games (payoffs 0..100).
The earlier 13% and 7% figures were for random dimensions 1..9 and 1..15.

**Precision policy (user, 2026-11-23):** the caller should choose the data
type, and on overflow it simply throws. The fixed-type solvers
`Solver64/128/MP/GMP` do exactly that. The `Solver` auto alias is an
optional extra; keep it optional.

## Decisions (answered by the user, 2026-11-23)

Decisions 1-8 use the API names of the time. In the names of decision 9:
`solve_fast` = `solve(const ConstantSumGame*)`, `solve_full` =
`solve(const GeneralSumGame*)`, and `solve_all` =
`extreme_equilibria(const GeneralSumGame*)`.

1. **`solve_fast` (constant-sum) default: `Mode::Any`**, the exact LP. The
   value always equals legacy's. In degenerate games the strategies can
   differ from legacy's "last" equilibrium; that is accepted.
2. **`solve_full` (general-sum) default: `Mode::Any`**, Lemke-Howson. This
   changes the provisional `Last`. It returns one exact equilibrium, usually
   not the one legacy returns. `Last` (legacy-exact) and `First` stay
   available as explicit modes.
3. **Switching `lib.h`/CMake and deleting the legacy `.cc` files: not yet.**
   Do the idiomatic-internals refactor (step 6.4) first, while `lib.h` stays
   on legacy.
4. **Type selection after the switch: one target, type in code.** A single
   CMake INTERFACE target `lrsnash`. Callers pick the type in code
   (`LRSNash::Solver64 s; s.solve(&game, &out)`). Downstream call sites and
   CMake will change when the switch happens; the four legacy target names
   are not kept.
5. **API extras: all three wanted.**
   - a caller-chosen Lemke-Howson missing label (default label 1);
   - exact rational output from `solve_all`;
   - rational `num/den` payoff input (lcm scaling per player).
6. **`solve_full` value convention: keep it**, as player 1's payoff / den.
7. **Performance work: all three wanted.**
   - per-call overhead (stack storage or size bounds) for the LP and LH
     defaults;
   - a branch-light `lrs_ratio`/`reverse` for Last/First;
   - an AVX2 kernel, which must use runtime dispatch (see 8).
8. **Downstream build flags: CMake Release**, which is `-O3 -DNDEBUG` with
   no `-march`. Benchmark at `-O3 -DNDEBUG`. Code that needs AVX2/BMI2 must
   use a target attribute plus runtime dispatch. The 2026-11-23 measurement:
   `-O3` equals `-O2` within noise; `-march=native` gives 11% (LP), 30% (LH)
   and 13% (Last) at 9x9 int64, and nothing at 15x15 int128.

9. **API names: renamed** (the user, 2026-11-23: "just do the rename").
   `solve_fast`/`solve_full` read as "first"/"all" equilibria, and the old
   API mixed two independent axes. Now:
   - `solve(const ConstantSumGame*, FloatEquilibrium*, Mode = Any)`;
   - `solve(const GeneralSumGame*, FloatEquilibrium*, Mode = Any,
     long missing_label = 0)`;
   - `extreme_equilibria(const ConstantSumGame*)` (new: this combination
     did not exist before);
   - `extreme_equilibria(const GeneralSumGame*)` (was `solve_all`).

   `ConstantSumGame`, `GeneralSumGame` and `FloatEquilibrium` are aliases of
   the legacy `lib.h` structs `FastInput`, `Input` and `FloatOneSumOutput`.
   Legacy `lib.h` keeps the old names until the switch, because downstream
   still builds it.

   The user approved the rename without reviewing the sub-choices. Claude
   applied them as it had recommended; they can be revisited:
   - `constant_sum`, not `zero_sum`: player 2 gets `den - data`, not
     `-data`;
   - `extreme_equilibria`, not `solve_all`: only the extreme equilibria are
     returned, and degenerate games can have infinitely many equilibria;
   - the name `FloatEquilibrium` for the output struct;
   - `Equilibrium` payoffs are now **divided by `den`**, like `value`
     (`solve_all` reported matrix units; only `tools/` used it);
   - pointers, not references, as in legacy;
   - `Mode` is unchanged.

Earlier decision: the caller chooses the arithmetic type and overflow
throws. The `Solver` auto alias stays as an optional extra (user,
2026-11-23).

## Open questions (for the user)

The user decides these in batches (2026-11-23). Every open design or
interface question goes here, as do decisions Claude took alone that the user
may want to revisit. Work continues around them. When one is answered, move it
to "Decisions" above.

**Q1. Public API names**: answered, see decision 9. The sub-choices it
lists were applied without review and stay revisitable.

**Q2. Exact rational output for "all"** (decision 5 says wanted; its form is
open). *Blocks:* that item of next step 4.
- Inside, the results are already exact: numerators over a common
  denominator, in the solver's type `T`.
- Options:
  - (a) a template visitor that gets `const T*` exactly as `nash.hpp` produces
    them: no copies, works for all 4 types, but the caller deals with `T`;
  - (b) return GMP rationals (puts GMP in the API);
  - (c) fixed-width num/den pairs that throw on overflow;
  - (d) strings: no (lean headers).
- For the auto chain `Solver<I64, I128, GMP>`, `T` depends on which type
  succeeded: either convert to the widest type, or call a generic-lambda
  visitor with whatever type succeeded.
- Recommend (a) with a generic lambda, and keep the double-valued convenience
  call.

**Q3. Rational payoff input** (decision 5 says wanted; its form is open).
*Blocks:* that item of next step 4.
- Inside: scale each player's matrix by the lcm of its denominators (best
  responses do not change), then divide the reported payoffs back. If a
  scaled entry does not fit in `long`, throw `overflow_error`.
- Open:
  - the input struct: `num`/`den` arrays per entry, or one denominator per
    player?
  - how it combines with today's scalar `den`, which both defines
    constant-sum games and scales `value`.

**Q4. Rejecting 0-size games** (decided by Claude, 2026-11-23; revisit if it
matters). Every `solve_*` throws `invalid_argument` unless rows, cols >= 1.
- This was not validated before. With `cols = 0` the old lrs code went down a
  rank-deficient path (it computed rays).
- Legacy `lib.h` never defined what 0 sizes do.
- Revisit if downstream can pass empty games.

**Q5. `den <= 0`** (noticed by Claude, 2026-11-23; nothing changed).
- `value` and the `Equilibrium` payoffs are divided by `den`, so `den = 0`
  gives inf/NaN. For a `ConstantSumGame` that is exactly the zero-sum case
  (B = -A). A negative `den` flips the sign. Legacy behaves the same.
- Options:
  - (a) reject `den <= 0` with `invalid_argument`;
  - (b) document "`den > 0`";
  - (c) for `den = 0`, report unscaled payoffs.
- Recommend (a), because legacy's normalisation assumes `den > 0`.

## Current status / next steps (handoff, 2026-11-23)
**User preferences.**
- Verify correctness by **static analysis against the paper and upstream
  source**. Randomized differential runs are secondary evidence.
- The caller picks the arithmetic type explicitly.
- Library headers stay lean: no `std::string` or similar (exception
  messages are `const char*` literals).
- Benchmark only on mains power with the `performance` power profile
  (`powerprofilesctl get`); otherwise do no timing and keep test trial
  counts low (a few hundred games).

**State.**
- Sessions 4 and 5 are committed on `rewrite` in the commit after
  `fff2737`. It is local only, not pushed; the user said "proceed". Later
  work stays uncommitted until the user asks.
- All tests pass:
  - `tools/test_lp.cc`, `tools/test_lh.cc`, `tools/test_api.cc`;
  - full-sequence equality with legacy for all 4 types;
  - exact NE checks;
  - `tools/test_refactor.cc`: bit-identical to the pre-refactor headers.
- The headers compile cleanly with g++ 15 and clang++ 21 at
  `-Wall -Wextra`, with and without `-DNDEBUG` (checked by building
  `tools/bench.cc`, which instantiates every solver path).
- `lib.h`/CMake still build only the legacy code.

**Header-only code (`src/lrsnash/`).**
- `arith.hpp`: the four policies `I64`, `I128`, `MP<W=112>` and `GMP`, with
  the fast kernels (`pivot`, `tableau_pivot`, `pivot_safe`,
  `magnitude_bound`, `from_i64/i128`).
- `lrs.hpp`: the lrs core (`Dict` and `Problem`), using the kernels and
  `memcpy` copies. Its header documents the dictionary layout, maps every
  function to upstream, and lists the removed dead code with the reason
  each piece is dead.
- `nash.hpp`: the lrsnash driver with a visitor (First, Last and All). Its
  header explains the two polyhedra and the label numbering.
- `zerosum.hpp`: the LP, with `promote_from` for the auto chain.
- `lemke.hpp`: Lemke-Howson.
- `solve.hpp`: `SolverT<Ar...>` and its aliases. It provides:
  - `solve(const ConstantSumGame*, out, Mode = Any)`: `Any` is the LP, and
    `Last` is exactly legacy `solve_fast`;
  - `solve(const GeneralSumGame*, out, Mode = Any, missing_label = 0)`:
    `Any` is Lemke-Howson, and `Last` is exactly legacy `solve_full`;
  - `extreme_equilibria(const ConstantSumGame*)` and
    `extreme_equilibria(const GeneralSumGame*)`: every extreme equilibrium,
    in upstream order, with both payoffs divided by `den`.

  It writes exactly `rows`/`cols` entries and has no globals. Every
  `solve_*` throws `std::invalid_argument` unless `rows, cols >= 1`.

**Functionality vs upstream 7.2.**
- Kept: all-equilibria enumeration in upstream order, and all 4 arithmetic
  types.
- Added: first-equilibrium, LP and LH modes.
- Removed when the code was originally stripped (before this work):
  rational num/den payoffs (could be restored by lcm-scaling each player's
  matrix), CLI and I/O options, the `MAXSTRAT=200` cap, and `2nash`.

**Next steps (in order, per the 2026-11-23 decisions).**
1. Done (session 4): small API changes. `solve_full` defaults to
   `Mode::Any`, and the Lemke-Howson missing-label argument is added.
2. Done (session 5): idiomatic internals (step 6.4) for `lrs.hpp` and
   `nash.hpp`. `lib.h` stays on legacy meanwhile.
3. Performance work (decision 7), benchmarked at `-O3 -DNDEBUG`.
   **Blocked until timing is possible** (mains power and the `performance`
   profile; the user said not to benchmark in session 5). First measure the
   session-5 pool fix for Last/First, which was not timed:
   - per-call overhead of the LP and LH (the defaults);
   - a branch-light `lrs_ratio`/`reverse` for the fixed-width types (P1
     reverse search is about half of Last-mode time);
   - an AVX2 kernel behind a target attribute with runtime dispatch.
4. API extras: exact rational output of `extreme_equilibria`, rational
   payoff input.
   **Blocked on the interface questions Q2 and Q3.**
5. Later, when the user says so: switch `lib.h` to one INTERFACE target
   `lrsnash` with the type chosen in code (6.1, 6.2), archive legacy, then
   delete it (6.3).

**Progress (2026-11-23, session 4).**
- Step 1 done: `solve_full` defaults to `Mode::Any`, with the
  `missing_label` argument (0-based, `invalid_argument` if out of range).
  `tools/test_lh.cc` checks it through the public API (0 mismatches).
- `tools/test_api.cc` now checks `Last` against legacy's exact sequence.
  It counts legacy's own float errors (`mptoi`) separately.
- Step 2 (refactor) is designed but **not started in code**. Plan:
  - one `pivot()` that also updates B/C positions;
  - descriptive names, with an upstream mapping table;
  - shared `start_basis`/`finish_basis` for player 1 and player 2;
  - drop dead code: the ray/linearity-space output (nredundcol is always 0
    by full column rank), `ismin` (dead with col 0), FillFirstRow, the
    origin check, counters and `minratio[m]`.

  Gate: compare exact outputs and pivot counts against a frozen copy of the
  pre-refactor headers (namespace `LRSNashRef`; rebuild it from git if the
  scratchpad copy is gone).
- Found a perf bug for the refactor to fix: `Problem::init` clears the
  dictionary pool on every solve, because `pool[0]->d_orig` was changed by
  `resize`. So Last/First reallocate every dictionary each solve.
- 50x50 constant-sum LP, payoffs 0..100, `-O3`:
  - int64 and int128 overflow on every game;
  - GMP 21.6 ms, auto 19.8 ms, lrsmp 133 ms;
  - 30x30: auto 1.17 ms, GMP 2.09 ms;
  - 50x50 Lemke-Howson: GMP 32 ms, auto 15.6 ms.

**Progress (2026-11-23, session 5).**
- Step 2 done. `lrs.hpp` went from 773 to 651 lines and `nash.hpp` from 455
  to 272.
  - `pivot()` now includes upstream's `update()`, and takes `long&`
    instead of `long*`.
  - Descriptive names, with an upstream mapping table in the header of
    `lrs.hpp`.
  - Player 1 and player 2 share `start_order`, `find_basis`,
    `remove_linearity_columns` and `finish_basis`. P2 differs only in
    saving `P2base` (was `P2orig`) between `find_basis` and the column
    removal, and in the warm restart.
  - Dead code is removed, each piece with its static reason in the header
    of `lrs.hpp`: rays and redundant columns, `ismin`, the
    `set_row`/`reducearray` gcd reductions, `FillFirstRow`, the origin
    checks, `depth`/`MAXD`, `count[]` except pivots, `minratio[m]`, `m_A`,
    `d_orig`, `Gcd`/`Lcm`, `temparray` and `reorder`.
  - Kept on purpose: 1-based basis positions, because `B[0]` is the
    objective row. That is the dictionary's structure, not a C habit.
- Facts the removal relies on, which the code now enforces:
  - No redundant column. Both systems have an identity block plus payoff
    rows whose payoff-variable coefficient is 1, so the column rank is
    `nvars` when `rows, cols >= 1`. `find_basis` throws `logic_error`
    otherwise. The new `solve_*` dimension check makes that unreachable;
    before, `cols = 0` was unvalidated and the old code computed a ray.
  - Every input row has gcd at most 1 (upstream divides only by a gcd
    above 1). This is a debug-only `assert` in `set_row`.
- The pool bug is fixed: `init` now compares only `m` and `stride`. A
  repeat solve of the same shape allocates no dictionary. This is checked
  through `Nash::dict_allocations()` in `test_refactor`; a mutant with the
  old behaviour fails 50 of 50 games. The speed gain is **not measured**
  (no benchmarking this session).
- One bug was made and fixed during the refactor, caught by the gate, not
  by static review. The lookup of a linearity column must use the
  *current* `d` (`r + P->d`): every removal shifts the higher indices down
  by one, and upstream relies on that. Watch for this in later reviews of
  index arithmetic.
- New gate `tools/test_refactor.cc`. It checks the exact visitor event
  sequence, early stop after 1-3 events, pivot count and overflow
  behaviour against `fff2737` (namespace `LRSNashRef`), for all 4 types,
  plus the repeat-solve check. Results:
  - before the user asked for low counts: 25,000 games (maxdim 4, 9 and
    14), ~7.8M pivots, 0 mismatches;
  - after the final edits: 300-400 games, 0 mismatches.
- Gate-script lesson: one build failed (`Nash::T_ar` was dropped, and
  `solve.hpp` uses it), but the tests still "passed" by running the stale
  binaries. The script must delete the old binaries and check every build's
  exit status (see "How to re-verify").
- A second static pass compared every function line by line with the
  frozen reference. The points it needed:
  - `ratio_test` now keeps `Nmin`/`Dmin` as pointers into `A`. This is
    safe because every `comprod` takes `const&`, writes only member
    temporaries, and each pass sets them at its first candidate.
  - The warm restart skips upstream's zeroing of row 0. `P2base` is saved
    with a zero row 0, and pivots keep a zero row zero.
  - P1's "linearity not cobasic → fail" branch is unreachable.
    `find_basis` brings only decision variables in, so after compaction
    every remaining linearity is cobasic.

**Progress (2026-11-23, session 5, after commit `1ae7810`).**
- The API rename was applied (decision 9) as its own commit.
  - `tools/test_api.cc` now also checks that, on constant-sum games, the
    `ConstantSumGame` overloads (`extreme_equilibria`, `solve` Last/First)
    give exactly the `GeneralSumGame` results on `(a, den - a)`.
  - `tools/test_lh.cc` checks that every overload rejects 0-size games.
  - `tools/bench.cc` modes are now `{cs,gs}-{any,first,last,all}`, with the
    old names kept as aliases. Games are drawn as before, so results stay
    comparable.
  - All gates passed at low counts (`test_lp` 300, `test_lh` 150, `test_api`
    200, `test_refactor` 300 games): 0 mismatches.
- Work-plan step 3 gap: the completeness check for degenerate games is
  written: `tools/dump_eq.cc` and `tools/check_complete.py` (brute-force
  vertex pairs in exact rationals). They are in the working tree, not
  committed. Results are below once it has run.

Not planned: a fixed 256/512-bit type (only helps at 30x30 and larger; it
is the item to revisit if 50x50 matters).
Caveat: Lemke-Howson has exponential worst-case paths (Savani–von
Stengel), though random games are short.

**How to re-verify quickly** (repo root; legacy libraries from
`cmake -S . -B build && cmake --build build`). Delete old test binaries
first and stop on any build error, so that stale binaries cannot pass:
```
g++ -O2 -std=c++17 -Isrc tools/test_lp.cc build/liblrslib.a -lgmp -o test_lp && ./test_lp 3000 9 1 legacy
g++ -O2 -std=c++17 -Isrc tools/test_lh.cc -lgmp -o test_lh && ./test_lh 1000 7 1
python3 tools/instrument_legacy.py /tmp/legacy && (cd /tmp/legacy && for f in lib lrslib lrsgmp; do g++ -O2 -c -I. $f.cc -DGMP -w -fpermissive; done)
g++ -O2 -std=c++17 -Isrc tools/test_api.cc /tmp/legacy/*.o -lgmp -o test_api && ./test_api 1000 8 1
mkdir -p /tmp/ref/lrsnash_ref && for f in arith lrs nash; do git show fff2737:src/lrsnash/$f.hpp | sed 's/\bLRSNash\b/LRSNashRef/g' > /tmp/ref/lrsnash_ref/$f.hpp; done
g++ -O2 -std=c++17 -Isrc -I/tmp/ref tools/test_refactor.cc -lgmp -o test_refactor && ./test_refactor 3000 9 1
g++ -O2 -std=c++17 -Isrc tools/bench.cc -lgmp -o bench && ./bench 9 9 80 10000 100000001 cs-any auto
```
Every mismatch count should be 0. The exception is the
`fast_test`-style comparison of `solve_fast` against the legacy *value*: it
shows rare mismatches on 9x9 games with payoffs 0..1000. Those are the legacy
`mptoi` value bug, not a regression.

**Tooling in `tools/`.** The build and run commands are in each file's header.
- `bench.cc`: benchmark of the header-only solvers, run as
  `./bench rows cols den trials seed {cs,gs}-{any,first,last,all} auto|i64|i128|mp|gmp`
  (`cs` constant-sum, `gs` general-sum; `all` = `extreme_equilibria`). The
  old mode names still work: `fast` = `cs-any`, `full` = `gs-last`, and
  `first`, `any`, `all` = `gs-…`.
- `test_lp.cc`, `test_lh.cc`: exactness tests for the LP and for
  Lemke-Howson (kernel bit-identity, exact NE checks in GMP).
- `test_api.cc` + `instrument_legacy.py`: solve.hpp vs the legacy library,
  checking full output sequences for all 4 types.
- `test_refactor.cc`: `nash.hpp`/`lrs.hpp` vs the frozen pre-refactor
  headers (`fff2737`). It checks exact events, pivots, early stop and
  overflow, and that a repeat solve allocates nothing. Use it as the gate
  for any further internal change.
- Older tools (copies; the originals are in the session scratchpad):
- `extract.py`: per-function normalized diff, run as
  `python3 extract.py upstream.c stripped.cc outdir`.
- `macros.py`: header `#define` diff.
- `diff.cc`: differential test against an upstream binary. Build it against
  one of the libraries, then run
  `./diff upstream_bin fast|full trials maxdim maxval seed`.
- `verify.py`: exact Fraction NE checker plus support enumeration, run as
  `python3 -I verify.py upstream_gmp_bin trials seed`.

To build the upstream binaries, run this in `~/lrslib-072`:
```
gcc -O2 -std=gnu99 -w -fpermissive -o nash_gmp lrsnash.c lrsnashlib.c lrslib.c lrsgmp.c lrsdriver.c -lgmp -DGMP
gcc ... lrslong.c -DLRSLONG -DSAFE [-DB128]
gcc ... -DMP -DLRS_QUIET lrsmp.c
```
Upstream crashes on long input paths, so use short relative paths and
`ulimit -c 0`.
