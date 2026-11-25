# lrsnash (stripped) — working context

This file holds the background, goals, constraints and plan for the ongoing
verification, optimization and C++ rewrite work on this repo. Read it before
making changes. `PERFORMANCE.md` explains the performance changes, biggest
first: what each one does, why it keeps results exact, and what it gained.

**Branches.** All work since the initial stripped port lives on the `rewrite`
branch: this file, `PERFORMANCE.md`, `src/lrsnash/` and the new `tools/`.
`main` holds only the legacy code. Since session 7, `lib.h` on `rewrite` is the
header-only API and the legacy sources are deleted there; they are kept at
the git tag `legacy` (`e6f5246`, = `main`).

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
- This repo (branch `rewrite`, header-only since session 7):
  - `src/lib.h` — the public header that downstream includes. It documents
    the API and includes `lrsnash/solve.hpp`.
  - `src/lrsnash/` — the implementation:
    - `arith.hpp` — the four arithmetic policies (`I64`, `I128`, `MP<W>`,
      `GMP`) with the fast exact pivot kernels.
    - `lrs.hpp` — the lrs reverse-search core (`Dict`, `Problem`).
    - `nash.hpp` — the lrsnash driver (exact upstream equilibrium order).
    - `zerosum.hpp` — the exact LP for constant-sum games.
    - `lemke.hpp` — exact Lemke-Howson.
    - `solve.hpp` — the game and output structs and the `SolverT<Ar...>`
      entry points.
  - `CMakeLists.txt` — one INTERFACE target `lrsnash` (links GMP if found).
    `-DLRSNASH_BUILD_TOOLS=ON` adds the benchmark and three CTest tests.
  - `tools/` — benchmarks, exactness tests and differential tests (see
    "Tooling" below).
  - The legacy library (`src/lib.cc`, `src/lrslib.*`, `src/lrsmp.*`,
    `src/lrslong.*`, `src/lrsgmp.*`, the root `test.cc` and `benchmark.cc`)
    is at the tag `legacy`. `tools/instrument_legacy.py` extracts it for
    the differential tests.

## Current public API (`src/lib.h`, since session 7)

- Types (`solve.hpp`): `ConstantSumGame{rows, cols, data, den}` (player 2
  gets `den - data`), `GeneralSumGame{rows, cols, den, row_data, col_data}`,
  `RationalGame{rows, cols, row_num, row_den, col_num, col_den}`, and the
  output `FloatEquilibrium{row_strategy, col_strategy, value}`. The integer
  structs have the legacy fields; their input pointers are now `const`.
- `SolverT<Ar...>`, with the aliases `Solver64`, `Solver128`, `SolverMP`,
  `SolverGMP` (one type; overflow throws) and `Solver` (int64, int128,
  then GMP). Keep one per thread and reuse it.

  |                   | one equilibrium (`Mode`)   | every extreme equilibrium |
  |---|---|---|
  | `ConstantSumGame` | `solve(&g, &out, mode)`    | `extreme_equilibria(&g[, f])` |
  | `GeneralSumGame`  | `solve(&g, &out, mode, k)` | `extreme_equilibria(&g[, f])` |
  | `RationalGame`    | `solve(&g, &out, mode, k)` | `extreme_equilibria(&g[, f])` |

  - `Mode {Last, First, Any}`. The default `Any` is the LP for constant-sum
    games and Lemke-Howson (missing label k) otherwise; `Last` is exactly
    legacy `solve_fast` / `solve_full`.
  - `extreme_equilibria(&g)` returns `std::vector<Equilibrium>` (doubles).
    With a callback f it calls `f(const ExactEquilibrium<Ar> &)` once per
    equilibrium, with the exact integers (decision 11).
- `value` and the payoffs are divided by `den` (by the lcm of each player's
  denominators for a `RationalGame`, i.e. true units). Exactly `rows` and
  `cols` entries are written. Errors are exceptions: `invalid_argument` for
  bad input (rows or cols < 1, den < 1, a denominator < 1, a bad missing
  label), `LRSNash::overflow_error` when a result does not fit.

The legacy API (`solve_fast`, `solve_full`, `FastInput`, `Input`,
`FloatOneSumOutput`, four static libraries chosen by macros) is gone; `lib.h`
lists the replacements. It was not header-only, and it returned only the
last equilibrium.

## Planned API direction

- Let the caller choose between **one equilibrium** and **all equilibria**.
  - For constant-sum games any equilibrium is fully correct. The value is
    unique and equilibria are interchangeable, so solving one exact LP (or
    stopping early) is a pure speedup.
  - For general-sum games (`solve_full`), equilibria are **not** unique in value
    or strategy. Returning only one (legacy keeps the last one) is a policy
    choice, not "the" answer. An "all" mode returns every extreme equilibrium
    (x, y) together with both players' payoffs.
- **Implemented** in `src/lrsnash/solve.hpp`, exposed through `lib.h` since
  session 7 (see "Current public API"):
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
   enumeration passed. Completeness of the extreme-equilibrium set on
   degenerate games is checked too (session 5, `tools/check_complete.py`;
   see the results log).
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

   Second round (session 6; timed, finished and committed in session 7):
   - Lemke-Howson on dictionaries;
   - int128 kernel fixes;
   - a branch-light `lrs_ratio`/`reverse`;
   - an AVX2 kernel with runtime dispatch.

   Not done: a fixed 256/512-bit type. More ideas are in PERFORMANCE.md
   section 7.
6. **C++ rewrite.** Done in session 7 for 6.1-6.6, except as noted. Turn the
   header-only port into the only
   implementation, written as idiomatic modern C++ (C++17). It must keep the
   hard constraints: header-only, no globals, reentrant, and all four
   arithmetic types.
   1. **Public API.** Done (session 7, decision 10). Switch `src/lib.h` to
      the header-only solvers. The
      structs keep their fields, under the new names `ConstantSumGame`,
      `GeneralSumGame` and `FloatEquilibrium` (aliases of the legacy structs
      until then; at the switch they become the definitions and the legacy
      names go, decision 9). Expose explicit arithmetic choice, `Mode`,
      `extreme_equilibria`, and a reusable per-thread solver object.
      Deferred until after 6.4 (decision 3); the type is chosen in code
      (decision 4).
   2. **Build.** Done (session 7). Replace the four static libraries with one
      CMake INTERFACE
      target, `lrsnash` (decision 4). Turn `tools/test_*` and `bench` into
      CTest/optional targets.
   3. **Remove legacy.** Done (session 7; tag `legacy`). Delete `src/lib.cc`,
      `src/lrslib.*`, `src/lrsmp.*`,
      `src/lrslong.*` and `src/lrsgmp.*` once the switch is done. Before that,
      capture any legacy-vs-new regression data still needed.
      `tools/instrument_legacy.py` needs the legacy sources, so archive a
      legacy build or tag `main` first.
   4. **Idiomatic internals.** Session 7 removed the dead lrs operations
      from `arith.hpp` (`gcd`, `lcm`, `exactdiv`, `mul`, `one`, including
      lrsmp's `goto`-based gcd). Done for `lrs.hpp` and `nash.hpp` (session
      5, see the progress notes). They used to mirror upstream C line by
      line, with upstream names, `long*` out-parameters and magic row
      indices. Behaviour stays **bit-identical**. Gates for every further
      refactor:
      - `tools/test_refactor.cc` (against the frozen pre-refactor headers);
      - the full-sequence differential test (`tools/test_api.cc`, against the
        archived legacy build);
      - `tools/test_lp.cc` and `tools/test_lh.cc`;
      - a static review against the paper and upstream (user preference).
   5. **Errors and output.** Done (session 7: the exact visitor, decision
      11). Exceptions only (no printf, exit or abort);
      already true for the new code. Exactly `rows`/`cols` output entries.
      Exact (rational) output for `extreme_equilibria` if wanted (Q2).
   6. **Optional:** rational payoffs done (session 7, `RationalGame`,
      general-sum only, Q9); stack storage not wanted (decision 13). Restore
      rational `num/den` payoffs (scale each player's
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
- **Completeness on degenerate games (session 5): PASSED.**
  - Static argument: lrsnash reports (x, y) for each vertex x of P1 (once,
    by lexmin) and each vertex y of the face F(x) of P2 on which the labels
    missing at x are tight. A face's vertices are exactly P2's vertices
    lying in it, so the output is exactly the set of completely labelled
    vertex pairs, i.e. the extreme equilibria.
  - The subtle point is degenerate x. Its labels must be read from the
    point, not from the basis: a basic slack with value 0 is tight. nash.hpp
    does this, because only basic slacks with a nonzero value become face
    equations.
  - Independent check: `tools/check_complete.py` enumerates both polyhedra's
    vertices by brute force (every (k-1)-subset of constraints, exact
    `Fraction`s, in Python) and compares the completely labelled pairs with
    what nash.hpp reports (`tools/dump_eq.cc`, GMP).
  - Results: 300 games up to 5x5 (190 degenerate, including all-zero
    games), 824 extreme equilibria; 0 missing, 0 extra, 0 duplicates.
  - The checker is shown to catch errors: hiding every third reported
    equilibrium fails 130 of 200 games, and reintroducing the session-5
    linearity-index bug fails 118 of 200.

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

**Latent issue, harmless** (the new code reads only player 1's
`slack_row` and says so in `nash.hpp`). In the `!FirstTime` branch, redundant
linearities stay in `Q2->linearity`. So `lrs_getfirstbasis2` compacts
`Q2->inequality` as if they had been removed, and the mapping becomes off
by one. Q2's `inequality` is never used after that, because player 2's
output uses only decision variables, so results are unaffected. Do not start
relying on `Q2->inequality` in the rewrite without fixing this.

**Stripped vs. upstream.** This was a function-by-function static diff (see
step 2). Every difference is removed I/O, a removed option that nash never
sets, or `exit` replaced by `throw`.

### Confirmed bugs in the current stripped code (fix in the rewrite)
**Resolved (session 7):** all of these are bugs of the legacy code, which is
deleted (tag `legacy`). The header-only code was checked statically for each:
no `printf`/`exit`/`abort`, no globals or `thread_local`, every allocation
owned (`std::vector`, `std::unique_ptr`), exactly `rows`/`cols` entries
written, and exact values converted to double directly.
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
**Resolved (session 7):** these are about the deleted legacy code. Each was
confirmed and fixed in the header-only code, proven there (the last-pair
consistency), or no longer applies (`MAXCOL`, `lrs_set_row_constraint`,
`reorder`, the CMake test targets).

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
  this should be proven. **Proven for the new code (session 5):** nash.hpp
  calls `on_eq1(x)` right after the `on_eq2(y)` calls of x's face whenever
  there was at least one, including when an `on_eq2` stopped the search. So
  `FloatEq` ends on the last pair reported (for First, the first), and
  `value` comes from that pair's y.
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
**Resolved (session 7):** the legacy code is deleted; the header-only code
has no global or `thread_local` state (hard constraint 2).
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

**Decisions 10-17 (answered by the user, 2026-11-24).** They answer the
`lib.h` switch (next step 5) and Q2-Q8. The user took every recommended
option except for the performance follow-ups, where they chose "both".

10. **Switch `lib.h` and CMake now, and delete legacy.** This replaces
    decision 3's "not yet".
    - First tag `main` as the legacy archive. The differential tests keep
      using a legacy build made from that tag.
    - Then switch `lib.h` and CMake to one INTERFACE target (decision 4).
    - Delete `src/lib.cc`, `src/lrslib.*`, `src/lrsmp.*`, `src/lrslong.*`
      and `src/lrsgmp.*` in the same change.
11. **Exact output of `extreme_equilibria` (was Q2): a visitor with raw `T`
    values.**
    - An overload calls the caller's callback (e.g. a generic lambda) with
      `const T*` numerators and denominators, exactly as computed, for all
      4 types, without copies.
    - For the auto chain, the callback gets whichever type succeeded.
    - The double-valued version stays.
12. **Rational payoff input (was Q3): per-entry num/den arrays.**
    - New game structs with `row_num`/`row_den`/`col_num`/`col_den` arrays,
      rows x cols each.
    - Each player's matrix is scaled by the lcm of its denominators. A
      scaled entry that does not fit in `long` throws `overflow_error`.
    - Payoffs and `value` are reported in true units, so there is no
      separate `den`.
13. **Solver storage (was Q6): heap storage, as now.** Document "keep one
    solver per thread and reuse it".
14. **0-size games (was Q4): keep throwing `invalid_argument`** unless
    rows, cols >= 1.
15. **`den <= 0` (was Q5): reject it with `invalid_argument`** in every
    overload. A zero-sum game can be passed with any `den > 0`: the
    strategies are the same, and `value` is scaled by 1/den.
16. **AVX2 detection (was Q7): `__builtin_cpu_supports`** when each I64
    object is constructed, as now.
17. **Performance follow-ups (Q8's open threshold and the thin-dispatcher
    lead): do both, and keep whatever wins.**
    - Time the AVX2 cut-over (`kAvx2MinCols`, now 8 columns).
    - Try a thin `I64::pivot` dispatcher for the small-game residual
      (PERFORMANCE.md 6.3).

## Open questions (for the user)

The user decides these in batches (2026-11-23). Every open design or
interface question goes here, as do decisions Claude took alone that the user
may want to revisit. Work continues around them. When one is answered, move it
to "Decisions" above.

Answered so far:
- Q1, public API names: decision 9. The sub-choices it lists were applied
  without review and stay revisitable.
- Q2 (exact output), Q3 (rational input), Q4 (0-size games), Q5
  (`den <= 0`), Q6 (inline storage), Q7 (AVX2 detection) and Q8 (AVX2 by
  default and its threshold): decisions 11-17, 2026-11-24.

**Q9. A constant-sum `RationalGame`** (raised by Claude, 2026-11-24).
- `RationalGame` is general-sum only, as decision 12 describes it (row and
  column arrays). A constant-sum game with rational payoffs can be passed
  as one, with `col = -row`, but then `Mode::Any` is Lemke-Howson, not the
  faster LP.
- Options: (a) leave it; (b) add a `RationalConstantSumGame{rows, cols, num,
  den}` (zero-sum: player 2 gets `-num/den`) that uses the LP.
- Recommend (a) until rational constant-sum games actually occur.

**Choices Claude made alone in session 7** (revisitable):
- The integer game structs take `const int *` (legacy: `int *`). Callers'
  code still compiles.
- The legacy root `test.cc` and `benchmark.cc` were deleted with the legacy
  sources (at the tag), not ported: `tools/bench.cc` and the tests cover
  them, and `test.cc` never called its own check.
- `RationalGame` requires |scaled payoff| <= 2^62, not just "fits in
  `long`": the LP and Lemke-Howson shift payoffs by `min - 1`, and the
  enumeration negates them, which overflow near the limits of `long`.
- `ExactEquilibrium` carries `scale1`/`scale2` (den, or the lcm for a
  `RationalGame`), so exact payoffs can be put in the double API's units.

## Current status / next steps (handoff, updated 2026-11-24)
**User preferences.**
- Verify correctness by **static analysis against the paper and upstream
  source**. Randomized differential runs are secondary evidence.
- The caller picks the arithmetic type explicitly.
- Library headers stay lean: no `std::string` or similar (exception
  messages are `const char*` literals).
- Benchmark only on mains power with the `performance` power profile
  (`powerprofilesctl get`); otherwise do no timing and keep test trial
  counts low (a few hundred games).
- Do not time, or run CPU-heavy tests, while the user runs their own
  benchmarks (2026-11-24). Short single-threaded runs under `nice -n 19`
  are fine; callgrind instruction counts stand in for timing.

**State.**
- Committed on `rewrite`, local only, not pushed:
  - sessions 4 and 5 in `1ae7810` and `18834e9` (the user said "proceed");
  - session 6, the session 5 completeness check and session 7's first
    fixes in `b37a781`. The user approved that commit ("commit sure",
    2026-11-24). Its message notes a partial regression, which the
    follow-ups below then fixed.
- The rest of session 7 is in the commit after `b37a781` (the user: "ok
  commit"). Later work stays uncommitted until the user asks.
  - Decisions 10-17 and Q9 in this file; the PERFORMANCE.md update.
  - Decision 17's follow-ups: in `arith.hpp`, `always_inline` on
    `I64::pivot` and `kAvx2MinCols` 8 -> 7.
  - The `lib.h` switch (decision 10): new `lib.h` and `CMakeLists.txt`, the
    structs defined in `solve.hpp`, the legacy sources and the root
    `test.cc`/`benchmark.cc` deleted (`git rm`; tag `legacy` made).
  - `den < 1` rejected (decision 15), the exact visitor (decision 11),
    `RationalGame` (decision 12), and the dead `arith.hpp` operations
    removed.
  - tools: `instrument_legacy.py` reads the tag; `test_lp.cc`
    (`-DLRSNASH_TEST_LEGACY`), `test_api.cc` and `diff.cc` include
    `legacy/lib.h`; new checks in `test_lh.cc` (den, `RationalGame`) and
    `test_api.cc` (exact visitor).
- All gates pass on the working tree (session 7, after the switch, the
  visitor and `RationalGame`; `gates.sh` scale 10, at `-O2` with asserts and
  at `-O3 -DNDEBUG`; the kernel was also gated with and without
  `-DLRSNASH_NO_AVX2`). CTest passes (3 tests):
  - `tools/test_kernel.cc`: the kernel contract against GMP;
  - `tools/test_lemke_ref.cc`: Lemke-Howson against the frozen `18834e9`
    copy;
  - `tools/test_lp.cc`, `tools/test_lh.cc`: exact NE checks;
  - `tools/test_api.cc`: full-sequence equality with legacy for all 4 types;
  - `tools/test_refactor.cc`: bit-identical to the pre-refactor headers;
  - `tools/check_complete.py`: completeness against brute force.
- The headers compile without warnings with g++ 15 and clang++ 21 at
  `-Wall -Wextra`, with and without `-DNDEBUG` and `-DLRSNASH_NO_AVX2`.
  This is checked by building `tools/bench.cc` (every solve path) and
  `tools/test_lh.cc` (`RationalGame` and the exact visitor).
- `lib.h` is the header-only API; CMake has one INTERFACE target.

**Header-only code (`src/lrsnash/`).**
- `arith.hpp`: the four policies `I64`, `I128`, `MP<W=112>` and `GMP`, with
  the fast kernels (`pivot`, `pivot_safe`, `magnitude_bound`,
  `from_i64/i128`). Since session 6 the I64 `pivot` also has an AVX2
  version, chosen at run time, and `tableau_pivot` is gone.
- `lrs.hpp`: the lrs core (`Dict` and `Problem`), using the kernels and
  `memcpy` copies. Its header documents the dictionary layout, maps every
  function to upstream, and lists the removed dead code with the reason
  each piece is dead.
- `nash.hpp`: the lrsnash driver with a visitor (First, Last and All). Its
  header explains the two polyhedra and the label numbering.
- `zerosum.hpp`: the LP, with `promote_from` for the auto chain.
- `lemke.hpp`: Lemke-Howson, on lrs dictionaries (cobasic columns only)
  since session 6.
- `solve.hpp`: the structs, `SolverT<Ar...>` and its aliases; see
  "Current public API". It writes exactly `rows`/`cols` entries and has no
  globals. The exact visitor reads player 1's vertex from the enumeration
  (`Nash::vertex1()`, computed before its face), so pairs are passed on
  without copies; only types before the last one of a chain buffer them.

**Functionality vs upstream 7.2.**
- Kept: all-equilibria enumeration in upstream order, and all 4 arithmetic
  types.
- Added: first-equilibrium, LP and LH modes.
- Removed when the code was originally stripped (before this work):
  rational num/den payoffs (could be restored by lcm-scaling each player's
  matrix), CLI and I/O options, the `MAXSTRAT=200` cap, and `2nash`.

**Next steps (in order, per the 2026-11-23 and 2026-11-24 decisions).**
1. Done (session 4): small API changes. `solve_full` defaults to
   `Mode::Any`, and the Lemke-Howson missing-label argument is added.
2. Done (session 5): idiomatic internals (step 6.4) for `lrs.hpp` and
   `nash.hpp`. `lib.h` stays on legacy meanwhile.
3. Performance work (decision 7), benchmarked at `-O3 -DNDEBUG`. **Done
   (sessions 6 and 7).** Rounds 1-2 are committed in `b37a781`; decision
   17's two follow-ups are uncommitted. It stops here unless the user asks
   for more speed: the hot paths are in single-digit-% territory (the 9x9
   LP takes about 1.77 µs with AVX2).
   - The pool fix was measured: Last -6%, First -7%.
   - Per-call overhead of the defaults:
     - Lemke-Howson now pivots dictionaries instead of full tableaus
       (-35% instructions at 9x9, -53% at 15x15 int128);
     - the LP's own setup is only about 8% of its cost;
     - inline storage for fresh solvers: not wanted (decision 13).
   - Branch-light reverse search: -18% simulated mispredicts in Last mode.
   - An AVX2 kernel with runtime dispatch (I64), on by default, for rows
     of 7 or more entries (decision 17).
   - Kernel fixes: the int128 regression from session 5, and the int128
     shift specializations.
   - **Final timing (session 7)**, against `18834e9`. The table is in
     the PERFORMANCE.md appendix and in the session 7 notes.
     - Default (AVX2) build: LP 6x6 -6%, 9x9 -10%, 15x15 -12% to -15%; LH
       6x6 -15%, 9x9 -29% to -31%, 15x15 -43% to -45%; Last/First -17% to
       -25%.
     - Nothing is measurably slower any more. The worst cases are 2x2 LH
       +1.8% (about 2.5 ns) and 4x4 LP +0.4%, within the 1-2% noise.
     - The partial regression in `b37a781`'s commit message (games up to
       6x6 1-6% slower; the plain-kernel 9x9 LP +4.5%) is fixed by decision
       17's follow-ups: the plain kernel is inlined into its callers, and
       the AVX2 cut-over is 7 entries.
   - The timing commands. Session 7 used `tcmp.sh` in its scratchpad, which
     does the same:

     ```
     O=/tmp/old && mkdir -p $O/src/lrsnash && for f in arith lrs nash zerosum lemke solve; do git show 18834e9:src/lrsnash/$f.hpp > $O/src/lrsnash/$f.hpp; done && git show 18834e9:src/lib.h > $O/src/lib.h
     # the same (current) bench.cc for both header trees: it prints 4 decimals
     g++ -O3 -DNDEBUG -std=c++17 -I$O/src tools/bench.cc -lgmp -o bench_old
     g++ -O3 -DNDEBUG -std=c++17 -Isrc tools/bench.cc -lgmp -o bench_new
     g++ -O3 -DNDEBUG -DLRSNASH_NO_AVX2 -std=c++17 -Isrc tools/bench.cc -lgmp -o bench_noavx2
     # per case, alternate the three binaries 5 times on one idle core, keep the best:
     for r in 1 2 3 4 5; do for b in bench_old bench_noavx2 bench_new; do taskset -c 4 ./$b 9 9 80 100000 100000001 cs-any i64; done; done
     ```
     Cases (rows cols den trials seed mode type):
     - `4 4 80 300000 . cs-any i64`, `6 6 80 200000 . cs-any i64`,
       `7 7 80 150000 . cs-any i64`, `9 9 80 100000 . cs-any i64|auto` and
       `15 15 100 20000 . cs-any auto|i128`;
     - `2 2 80 1000000 . gs-any i64`, `4 4 80 300000 . gs-any i64`,
       `6 6 80 150000 . gs-any i64`, `9 9 80 60000 . gs-any i64|auto` and
       `15 15 100 20000 . gs-any auto|i128`;
     - `6 6 80 2000 . gs-last i64`, `9 9 80 400 . gs-last i64`,
       `9 9 80 300 . gs-last i128`, `9 9 80 2000 . gs-first i64`, and
       `12 12 100 15 . gs-last i128`.
4. Done (session 7): reject `den <= 0` with `invalid_argument` in every
   overload (decision 15).
5. Done (session 7): `lib.h` and CMake switched to one INTERFACE target
   `lrsnash`, with the type chosen in code (6.1, 6.2; decision 10).
   `main` is tagged `legacy`, the differential tests build legacy from the
   tag, the legacy sources are deleted (6.3), and `lib.h` says "keep one
   solver per thread and reuse it" (decision 13).
6. Done (session 7): the exact visitor of `extreme_equilibria` (decision 11)
   and `RationalGame` (decision 12).
7. Committed (session 7). Downstream must now move to the new `lib.h` (its
   header lists the replacements for `solve_fast` / `solve_full`). Open: Q9.
   Optional: PERFORMANCE.md section 7 (the float-guided exact LP is the
   biggest lever left), and an idiomatic rewrite of the lrsmp port
   (`MP<W>`, used only by `SolverMP`).

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
- Work-plan step 3 gap closed: completeness on degenerate games PASSED
  (300 games, 190 degenerate, 0 mismatches; details in the results log).
  `tools/dump_eq.cc` and `tools/check_complete.py` are in the working tree,
  **not committed**.

**Progress (2026-11-24, session 6): performance work, next step 3.**
- **Timing.** The machine was on mains with the `performance` profile, but
  the user was running their own benchmarks: "refrain from tests that hog
  cpu". After the first native measurements below, the evidence is callgrind
  instruction counts (deterministic, single core, under `nice`) plus a static
  reading of the generated code. Later the user allowed full CPU use, and
  the final state was timed: see "Wall-clock results" at the end of this
  section, including **a known partial regression**.
- **First native measurements.** These were taken under background load,
  A/B interleaved, pinned to CPU 4, best of 5. Only the ratios are
  meaningful.
  - The session-5 pool fix (`fff2737` vs `18834e9`): 9x9 Last int64 672 ->
    634 µs (-6%), First 158 -> 148 µs (-7%).
  - **Regression:** int128 Last was 6-10% slower after the session-5
    refactor (9x9: 1,118 -> 1,185 µs; another round 1,200 -> 1,319 µs) and
    5% slower at 12x12, although `arith.hpp` did not change.
    - Cause: GCC inlined the rare exact 256-bit row (`pivot_row_exact`) into
      the kernel. The kernel grew from 541 to 1,130 instructions, and the hot
      loop gained an instruction per entry.
    - Fixed below.
  - The defaults, current code:
    - 9x9 LP: 2.05 µs (int64), 2.09 µs (auto);
    - 9x9 LH: 3.14 µs (int64), 3.21 µs (auto);
    - 15x15 LP: 12.55 µs (auto), 17.6 µs (int128);
    - 15x15 LH: 12.63 µs (auto), 28.35 µs (int128).
- **Kernel changes (`arith.hpp`).** All are bit-identical.
  - The rare exact rows (I64 `exact_row`, I128 `pivot_row_exact`) are
    `noinline, cold`, so they cannot take registers from the fast loop.
  - The row loop is specialized on the shift k:
    - k = 0 (det odd, about half of all pivots) needs no shift. Session 7
      removed this case again for I64 (see its notes);
    - for I128 with 0 < k < 64, the shift is `shrd` + `sar` instead of the
      general 128-bit shift (`test` + two `cmov`). GCC folds `k & 63` away
      when it knows k < 64 and then expands the general shift, so an empty
      `asm volatile` inside the loop hides k's range from it.
  - `ldexp`, a libm call on every I128 pivot, is replaced by `pow2`, an exact
    power of two built from its bits.
  - `tableau_pivot`, `dense_tableau_pivot` and the `LRS` template parameter
    are removed. They became dead once Lemke-Howson moved to `pivot`.
- **Lemke-Howson on dictionaries (`lemke.hpp`).**
  - Each polytope keeps only its cobasic columns, as an lrs dictionary,
    instead of the full tableau. The basic columns of a tableau are just det
    times unit vectors.
  - The lrs `pivot` kernel performs exactly the tableau update on the stored
    columns, because the pivot element is always negative here.
  - The lexicographic ratio test reads the tableau entries of the initial
    slacks from the dictionary: det or 0 for a basic slack.
  - The path and the result are the same. Kernel work per pivot drops from
    rows x (1+m+n) to rows x (1+m).
- **Instructions per game** (callgrind, `-O3 -DNDEBUG`, `18834e9` -> now):

  | case | before | now | change |
  |---|---|---|---|
  | 9x9 LP int64 | 20,300 | 19,539 | -3.7% |
  | 9x9 LP auto | 20,963 | 20,440 | -2.5% |
  | 9x9 LH int64 | 33,677 | 21,905 | -35.0% |
  | 15x15 LP auto | 141,358 | 125,777 | -11.0% |
  | 15x15 LP int128 | 210,053 | 175,114 | -16.6% |
  | 15x15 LH int128 | 436,459 | 204,642 | -53.1% |
  | 15x15 LH auto | 191,985 | 104,350 | -45.6% |
  | 9x9 Last int64 | 5,491,724 | 5,377,756 | -2.1% |
  | 9x9 Last int128 | 12,559,023 | 10,544,218 | -16.0% |
  | 9x9 First int64 | 1,195,589 | 1,173,981 | -1.8% |

- **Where the 9x9 LP's instructions go now.**
  - The kernel is 71%: an inner loop of 13.2 instructions per entry (~70%
    of the kernel), about 15 instructions per row and about 150 per pivot.
  - Pivot selection and ratio tests are 11%, setup 8%, float output 1.4%.
  - Unrolling or hoisting per-row work might save about 15% of kernel
    instructions. Whether that saves time depends on whether the loop is
    already bound by its two `imul`s per entry, so it needs wall-clock
    timing.
- **A fresh solver per game costs about 2k instructions more** (10-17% of a
  9x9 LP or LH solve): construction, `malloc`/`free` of about 20 vectors,
  and the destructor. A reused solver allocates nothing. See Q6.
- **New gates.**
  - `tools/test_kernel.cc` checks the kernel contract against GMP:
    - exact values;
    - overflow thrown if and only if a true result does not fit;
    - a sound returned bound.

    It runs on random Bareiss pivot sequences, and on constructed inputs
    whose results straddle the fast-path limit, for k = 0, 0 < k < 64 and
    k >= 64.
    - Mutation-checked: 6 of 6 injected kernel bugs fail it (a wrong shift
      for k >= 64, a missing bound accumulation, a loosened fast-path limit
      for I64 and for I128, a loosened all-rows test, a missing I128
      overflow check).
    - The first version missed the loosened limit, because random
      dictionaries almost never reach it; the constructed boundary cases were
      added for that.
  - `tools/test_lemke_ref.cc` compares Lemke-Howson with the frozen
    `18834e9` `lemke.hpp`, for every missing label and all 4 types: same
    pivots, exact solution, floats and overflow.
  - `test_lh`'s reference now overrides `pivot` (it was `tableau_pivot`).
  - All gates pass at `-O2` and at `-O3 -DNDEBUG`.
- **Branch-light reverse search (`lrs.hpp`).** All of these are
  bit-identical; `test_refactor` passed on 600 games up to 12x12 (233k
  pivots, all 4 types).
  - `ratio_test` collects its candidates without a branch:
    `cand[n] = j; n += negative(...)`.
  - For the fixed-width types (the new `Ar::fixed_width` constant),
    `ratio_test` and `is_reverse_pivot` compare with `comprod` directly.
    They no longer use upstream's sign shortcuts, which only avoid the
    product and give its sign: checked case by case, since both
    denominators are negative. `is_reverse_pivot` combines its test without
    short-circuit branches. lrsmp and GMP keep the shortcuts, because their
    products are expensive.
  - `resort` returns the element's new place. This replaces the two linear
    searches after every pivot; the lists are sorted and distinct, so the
    position is the same.
  - Measured on 9x9 int64 Last with callgrind's branch simulator:
    instructions -1.0%, conditional branches -7.8%, simulated mispredicts
    -17.9% (49.6k -> 40.7k per game).
  - The biggest remaining source (18k per game) is the kernel's inner-loop
    exit on 10-entry rows. Zen 3's predictor probably handles that better
    than the simulator does.
- **AVX2 kernel (I64, `arith.hpp`).**
  - `rows_pass_avx2` is the same row loop compiled with
    `__attribute__((target("avx2")))`. GCC and Clang vectorize it: 6
    `vpmuludq` per 4 entries for the two 64-bit products, and the
    arithmetic shift from a logical shift plus the sign. That is 7
    instructions per entry against 13-14 for the scalar loop.
  - The fast loop is a separate function `fast_row` with `__restrict`
    parameters. Without them GCC adds a runtime aliasing check to every row,
    and it ignores `__restrict` on block-scope pointers.
  - Dispatch: each I64 object sets `avx2 = __builtin_cpu_supports("avx2")`
    when it is constructed (Q7), and `pivot()` branches on that flag.
    `-DLRSNASH_NO_AVX2` compiles the AVX2 version out.
  - `test_kernel` runs the I64 contract with both kernels. All gates pass
    with AVX2 on and with `-DLRSNASH_NO_AVX2`, and the kernel, LH and
    refactor gates also pass under clang `-O3`.
  - Wall-clock speed is **not measured** (Q8).
  - Instructions per game with AVX2, against `18834e9` and against the same
    code with `-DLRSNASH_NO_AVX2` (callgrind, `-O3 -DNDEBUG`):

    | case | `18834e9` | no AVX2 | AVX2 |
    |---|---|---|---|
    | 9x9 LP int64 | 20,082 | -1.2% | -9.4% |
    | 9x9 LH int64 | 33,660 | -33.9% | -38.9% |
    | 15x15 LP auto | 138,735 | -7.6% | -16.2% |
    | 15x15 LH auto | 191,764 | -44.8% | -56.5% |
    | 9x9 Last int64 | 5,500,058 | -1.8% | -11.6% |
    | 9x9 First int64 | 1,197,914 | -1.6% | -11.3% |

    The `18834e9` column differs by 1-3% from the first table because the
    profiling driver changed in between. Compare within one table only.
  - On 9x9 the vector loop covers 8 of each row's 10 entries; the other 2
    take the scalar epilogue, and each row pays about 12 instructions of
    vector setup.
  - 15x15 LP int64 shows +2.2% without AVX2. Most of those games overflow,
    and the out-of-line `exact_row` adds a frame for every throw to unwind
    (the unwinder is 26% of that case). This is the error path, not the
    kernel.
- **Later changes from the small-game check (before committing):**
  - The first AVX2 version made small games slower in instruction count:
    LP +7.8% at 2x2, falling to +1.4% at 6x6. The dictionary LH was +14.5%
    at 2x2. Three changes followed:
  - **AVX2 dispatch** is now one test at the top of `I64::pivot`:
    `if (avx2 && cols >= kAvx2MinCols) pivot_avx2(...)`.
    - `pivot_avx2` is the whole `pivot_impl` compiled for AVX2, so the
      AVX2 path costs no extra call, and the per-pivot loops vectorize too.
    - `kAvx2MinCols = 8` (7 since session 7), chosen from instruction
      counts: shorter rows get at most one vector iteration.
    - `rows_pass_avx2` is gone. `test_kernel`'s boundary cases now go up to
      13 columns, so they reach the AVX2 path.
  - **`lemke.hpp` setup:** `col`/`row` use `resize` instead of `assign`
    (every label is written by `basic`/`cobasic`). Both sides start with an
    unknown bound, as `18834e9` did, instead of a `magnitude_bound` pass.
    **This last edit has not been gated yet.**
  - Instruction counts after these changes (`18834e9` -> no AVX2 / AVX2):

    | case | no AVX2 | AVX2 |
    |---|---|---|
    | LP 2x2 / 4x4 / 6x6 int64 | -1.4% / -1.3% / -1.2% | -0.6% each |
    | LP 9x9 int64 | -1.2% | -14.8% |
    | LP 15x15 auto | -7.6% | -18.1% |
    | LH 2x2 / 4x4 / 9x9 int64 | -4.5% / -19.3% / -34.8% | -3.6% / -19.0% / -44.1% |
    | LH 15x15 auto | -44.9% | -60.0% |
    | Last 6x6 / 9x9 int64 | -1.8% / -1.8% | -12.2% / -16.3% |
    | First 9x9 int64 | -1.6% | -16.1% |
    | Last 9x9 int128 | -13.4% | -13.4% |
    | Last 4x4 GMP | -0.2% | -0.2% |
    | LP 15x15 int64 (mostly overflow) | +2.2% | -18.8% |

- **Wall-clock results** (2026-11-24, after the user allowed full CPU use).
  - Setup: mains, `performance` profile, light load, `-O3 -DNDEBUG`, pinned
    to CPU 4.
  - The three builds alternated, best of 5, µs per game: `18834e9` (head),
    the working tree with `-DLRSNASH_NO_AVX2`, and the working tree.
  - `lemke.hpp` was already in its final state for this run.

  | case | head | no AVX2 | vs head | AVX2 | vs head |
  |---|---|---|---|---|---|
  | 4x4 LP int64 | 0.31 | 0.33 | **+6.5%** | 0.33 | **+6.5%** |
  | 6x6 LP int64 | 0.77 | 0.80 | **+3.9%** | 0.80 | **+3.9%** |
  | 9x9 LP int64 | 1.97 | 2.06 | **+4.6%** | 1.78 | -9.6% |
  | 9x9 LP auto | 1.98 | 2.07 | **+4.5%** | 1.80 | -9.1% |
  | 15x15 LP auto | 11.95 | 11.34 | -5.1% | 10.29 | -13.9% |
  | 15x15 LP int128 | 17.08 | 14.91 | -12.7% | 14.88 | -12.9% |
  | 2x2 LH int64 | 0.14 | 0.14 | 0 | 0.15 | **+7.1%** (at the timer's 0.01 µs resolution) |
  | 4x4 LH int64 | 0.48 | 0.45 | -6.2% | 0.45 | -6.2% |
  | 9x9 LH int64 | 3.02 | 2.37 | -21.5% | 2.11 | -30.1% |
  | 9x9 LH auto | 3.06 | 2.37 | -22.5% | 2.08 | -32.0% |
  | 15x15 LH auto | 12.08 | 8.37 | -30.7% | 6.64 | -45.0% |
  | 15x15 LH int128 | 26.51 | 15.17 | -42.8% | 15.16 | -42.8% |
  | 9x9 Last int64 | 610.9 | 546.9 | -10.5% | 458.6 | -24.9% |
  | 9x9 Last int128 | 1,208 | 996.6 | -17.5% | 1,001 | -17.1% |
  | 9x9 First int64 | 141.6 | 127.6 | -9.8% | 106.6 | -24.7% |
  | 12x12 Last int128 | 19,736 | 16,238 | -17.7% | 16,274 | -17.5% |

  - **AVX2 (Q8) is confirmed** where it runs. Against the plain kernel:
    9x9 LP -14% (1.78 vs 2.06 µs), 9x9 LH -11%, 15x15 LH -21%, 9x9 Last and
    First -16%. Keep it on.
  - **Known partial regression** (as session 6 left it; session 7 fixed
    part of it and traced most of the rest to code placement, see its
    notes): the plain int64 kernel is
    4-6.5% slower on constant-sum LPs at 4x4-9x9, although it executes 1.2%
    *fewer* instructions.
    - On AVX2 machines, games with 8 or more columns (LP 7x7 and up) use the
      AVX2 kernel and are faster.
    - Smaller LPs take the plain path on every machine and are 4-6.5%
      slower. LH, Last, First and all int128 paths are faster despite it.
  - Hypotheses, to check by bisecting. Time 4x4/6x6/9x9 cs-any with
    `-DLRSNASH_NO_AVX2` against variants of `arith.hpp`:
    1. The I64 `k == 0` specialization. It saves only one `sar` per entry,
       but adds a data-dependent branch per pivot (k is 0 about half the
       time, in no pattern) plus a second copy of the row loop. For I128 it
       saves 5 instructions per entry, so keep it there. Likely the main
       cause.
    2. The `Step` struct passed by reference, if not fully scalarized: 13
       stores and reloads per pivot.
    3. The `fast_row` helper or the out-of-line `exact_row` changing layout
       or registers in the row loop.
    4. Code alignment.

    Fix whatever is found, re-time, and if anything is still slower, record
    it here and in the commit message (user's request). Session 7 did
    this: hypothesis 1 was part of it, and most of the rest is
    hypothesis 4.
  - `tools/bench.cc` printed µs with 2 decimals, too coarse for 2x2 (0.14
    µs). Since session 7 it prints 4.

**Progress (2026-11-24, session 7): the regression, then the commit.**
- The final `lemke.hpp` edit is gated. All gates pass at `gates.sh` scale
  10 in four configurations: `-O2` with asserts and `-O3 -DNDEBUG`, each
  with and without `-DLRSNASH_NO_AVX2`.
- **The small-LP regression** (session 6 "Wall-clock results"). It was
  investigated A/B interleaved, best of 5, against `18834e9`, and stopped
  after three cheap steps, given the user's cost concern:
  1. **Drop the I64 `k == 0` specialization.** Kept.
     - Plain kernel: 4x4 +6.0% -> +3.5%, 6x6 +4.6% -> +1.4%, 9x9 +5.3% ->
       +5.2%.
     - The AVX2 path did not change.
     - The scalar inner loop is now instruction for instruction the one in
       `18834e9`.
  2. **Restore `18834e9`'s scalar kernel verbatim** for the plain path.
     The function compiled to the same size, to the byte, and 9x9 was still
     +4.5%. Not kept.
  3. **Build both with `-falign-functions=64 -falign-loops=64`.** 4x4 then
     shows 0.0%, 6x6 -0.1% and 9x9 +2.3%. Most of the gap is code
     placement.
  - `perf` cannot be used (`perf_event_paranoid` = 4), so there are no
    counter data.
- **Clang warnings fixed.** Clang 21 gave 6 `-Wbitwise-instead-of-logical`
  warnings on session 6's branch-free `is_reverse_pivot`. Each predicate
  now goes into a named `bool` first. GCC's code changed only in register
  allocation in `next_basis`, and Last/First timing moved by -0.2% to -1.2%,
  which is noise.
- **Timing of `b37a781`** (µs per game; percentages against `18834e9`;
  "plain" = `-DLRSNASH_NO_AVX2`). The follow-ups below changed this; the
  final table is in the PERFORMANCE.md appendix.

  | case | `18834e9` | plain | AVX2 (default) |
  |---|---|---|---|
  | 4x4 LP int64 | 0.3068 | **+2.4%** | **+4.8%** |
  | 6x6 LP int64 | 0.7538 | **+1.2%** | **+1.9%** |
  | 9x9 LP int64 | 1.9041 | **+4.5%** | -8.5% |
  | 9x9 LP auto | 1.9493 | **+2.8%** | -9.3% |
  | 15x15 LP auto | 11.77 | -6.1% | -13.8% |
  | 15x15 LP int128 | 16.46 | -13.0% | -12.4% |
  | 2x2 LH int64 | 0.1331 | **+2.4%** | **+6.3%** |
  | 4x4 LH int64 | 0.4579 | -7.0% | -2.6% |
  | 9x9 LH int64 | 2.927 | -22.3% | -29.7% |
  | 9x9 LH auto | 2.997 | -24.6% | -32.2% |
  | 15x15 LH auto | 11.64 | -30.8% | -44.5% |
  | 15x15 LH int128 | 25.45 | -42.7% | -42.6% |
  | 9x9 Last int64 | 591.1 | -11.6% | -25.5% |
  | 9x9 Last int128 | 1,162 | -17.3% | -17.5% |
  | 9x9 First int64 | 138.1 | -11.8% | -25.6% |
  | 12x12 Last int128 | 19,394 | -17.7% | -17.5% |

  Run-to-run noise is about 1-2%.
  - In the default build, games under 8 dictionary columns use the plain
    kernel too, yet they are slower there than in the plain build (4x4 LP
    +4.8% against +2.4%). Only the code around the kernel differs.
  - A lead that was not followed: in step 2's variant, `pivot` had become a
    small dispatcher, and GCC inlined the plain kernel into its callers. In
    the default build that variant had 4x4 and 6x6 LPs level with
    `18834e9` (-0.5% and -0.1%, one run). Decision 17 followed it up (below).
- **Decisions 10-17** were answered by the user. That closes every open
  question; see "Decisions".
- **Decision 17's follow-ups, both kept** (see also PERFORMANCE.md 6.1 and 6.3):
  - `always_inline` on `I64::pivot`, so the plain kernel is part of each
    caller and the AVX2 kernel stays out of line. Of three variants (GCC's
    choice, forced inline, forced out of line), only forced inline is both
    reliable and fast.
    - Default build: 4x4 LP +4.0% -> +0.1% against `18834e9`, 2x2 LH +5.2%
      -> +1.9%.
    - Plain build: 9x9 LP +4.1% -> -1.3%.
  - The AVX2 cut-over moves from 8 to 7 entries. Always-AVX2 against
    never-AVX2: at 7 entries LP -4.5%, LH -0.9%, Last -7.7%; at 6 mixed;
    at 5 slower. The 6x6 LP went from +2.0% to -6.2% against `18834e9`.
  - Final state: nothing is measurably slower than `18834e9` (worst: 2x2
    LH +1.8%, 4x4 LP +0.4%), and the 9x9 LP is -10.0%.
  - Gates pass at scale 10 in all 4 build configurations, and g++/clang++
    give no warnings.

Not planned: a fixed 256/512-bit type (only helps at 30x30 and larger; it
is the item to revisit if 50x50 matters).
Caveat: Lemke-Howson has exponential worst-case paths (Savani–von
Stengel), though random games are short.

**How to re-verify quickly** (repo root). The legacy library comes from the
`legacy` tag. Delete old test binaries first and stop on any build error, so
that stale binaries cannot pass:
```
python3 tools/instrument_legacy.py /tmp/leg && (cd /tmp/leg/legacy && for f in lib lrslib lrsgmp; do g++ -O2 -c -I. $f.cc -DGMP -w -fpermissive -o gmp_$f.o; done && for f in lib lrslib lrsmp; do g++ -O2 -c -I. $f.cc -DMP -DLRS_QUIET -w -fpermissive -o mp_$f.o; done)
g++ -O2 -std=c++17 -DLRSNASH_TEST_LEGACY -Isrc -I/tmp/leg tools/test_lp.cc /tmp/leg/legacy/mp_*.o -lgmp -o test_lp && ./test_lp 3000 9 1 legacy
g++ -O2 -std=c++17 -Isrc tools/test_lh.cc -lgmp -o test_lh && ./test_lh 1000 7 1
g++ -O2 -std=c++17 -Isrc -I/tmp/leg tools/test_api.cc /tmp/leg/legacy/gmp_*.o -lgmp -o test_api && ./test_api 1000 8 1
mkdir -p /tmp/ref/lrsnash_ref && for f in arith lrs nash; do git show fff2737:src/lrsnash/$f.hpp | sed 's/\bLRSNash\b/LRSNashRef/g' > /tmp/ref/lrsnash_ref/$f.hpp; done
g++ -O2 -std=c++17 -Isrc -I/tmp/ref tools/test_refactor.cc -lgmp -o test_refactor && ./test_refactor 3000 9 1
g++ -O2 -std=c++17 -Isrc tools/dump_eq.cc -lgmp -o dump_eq && python3 -I tools/check_complete.py ./dump_eq 200 4 1
g++ -O2 -std=c++17 -Isrc tools/test_kernel.cc -lgmp -o test_kernel && ./test_kernel 3000 1
git show 18834e9:src/lrsnash/lemke.hpp | sed 's/\bLRSNash\b/LRSNashRef/g' > /tmp/ref/lrsnash_ref/lemke.hpp
g++ -O2 -std=c++17 -Isrc -I/tmp/ref tools/test_lemke_ref.cc -lgmp -o test_lemke_ref && ./test_lemke_ref 400 12 1
g++ -O2 -std=c++17 -Isrc tools/bench.cc -lgmp -o bench && ./bench 9 9 80 10000 100000001 cs-any auto
```
Run the gates once more with `-DLRSNASH_NO_AVX2` added to every build, so
that the plain I64 kernel is covered too on an AVX2 machine. Sessions 6 and
7 used a script that does all of the above under `nice -n 19` and greps for
nonzero counts (`gates.sh` in the session scratchpad; recreate it from these
lines if it is gone). The quick subset without legacy or frozen references
is `cmake -S . -B b -DLRSNASH_BUILD_TOOLS=ON && cmake --build b && (cd b &&
ctest)`.
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
- `test_api.cc` + `instrument_legacy.py`: solve.hpp vs the legacy library
  (extracted from the `legacy` tag), checking full output sequences for all
  4 types, and the exact visitor against `extreme_equilibria`.
- `dump_eq.cc` + `check_complete.py`: completeness of the extreme
  equilibria against an independent brute-force vertex-pair enumeration in
  exact rationals, run as
  `python3 -I tools/check_complete.py ./dump_eq trials maxdim seed`.
  Keep maxdim at 5 or less: the brute force tries every (k-1)-subset of the
  constraints.
- `test_refactor.cc`: `nash.hpp`/`lrs.hpp` vs the frozen pre-refactor
  headers (`fff2737`). It checks exact events, pivots, early stop and
  overflow, and that a repeat solve allocates nothing. Use it as the gate
  for any further internal change.
- `test_kernel.cc`: the contract of the fixed-width pivot kernels (I64
  plain and AVX2, I128) against GMP:
  - exact values;
  - overflow if and only if a true result does not fit;
  - a sound returned bound.

  It uses random Bareiss pivot sequences and constructed inputs at the
  fast-path limit. Run it for any kernel change. Mutation-checked in
  session 6.
- `test_lemke_ref.cc`: Lemke-Howson vs a frozen copy (`18834e9`, namespace
  `LRSNashRef`), every missing label, all 4 types. Use it as the gate for
  any change to `lemke.hpp`.
- Older tools (copies; the originals are in the session scratchpad):
- `extract.py`: per-function normalized diff, run as
  `python3 extract.py upstream.c stripped.cc outdir`.
- `macros.py`: header `#define` diff.
- `diff.cc`: differential test of the legacy library against an upstream
  binary. Build it against the legacy objects from `instrument_legacy.py`
  (see its header), then run
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
