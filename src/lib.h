// lrsnash: exact Nash equilibria of two-player (bimatrix) games, header-only.
// Derived from lrsnash in lrslib 7.2 (David Avis; GPL v2, see COPYING).
//
//   #include "lib.h"
//   LRSNash::Solver64 solver;  // keep one per thread and reuse it
//   LRSNash::ConstantSumGame g{rows, cols, data, den};
//   LRSNash::FloatEquilibrium out{row_strategy, col_strategy, 0.f};
//   solver.solve(&g, &out);
//
// Solver64, Solver128, SolverMP and SolverGMP use one arithmetic type and
// throw LRSNash::overflow_error if a result does not fit; Solver tries
// int64, then int128, then GMP. A reused solver allocates nothing after
// warm-up; constructing one per game costs 10-17% of a 9x9 solve. The full
// API (Mode, missing labels, extreme_equilibria) is in lrsnash/solve.hpp.
//
// From the legacy API: solve_fast(&g, &out) is solver.solve(&g, &out) with
// g a ConstantSumGame (the exact LP; Mode::Last gives legacy's equilibrium),
// and solve_full(&g, &out) is solver.solve(&g, &out, LRSNash::Mode::Last)
// with g a GeneralSumGame (the default, Mode::Any, is Lemke-Howson). The
// output arrays need exactly rows and cols entries.
#pragma once

#include "lrsnash/solve.hpp"
