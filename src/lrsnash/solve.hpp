// Header-only solver entry points.
//
// SolverT<Ar...> owns all workspaces (no globals, no thread_local), so one
// solver per thread gives allocation-free, reentrant solves after warm-up.
// The arithmetic types are tried in order: a solve runs with the first one,
// and if a result does not fit (overflow_error, detected exactly) it reruns
// with the next. With a single type, that type is used and overflow_error
// propagates. All four lrs arithmetic types are supported:
//
//   Solver64  = SolverT<I64>            64-bit  (upstream lrsnash1)
//   Solver128 = SolverT<I128>           128-bit (upstream lrsnash2)
//   SolverMP  = SolverT<MP<>>           lrsmp   (upstream -DMP build)
//   SolverGMP = SolverT<GMP>            GMP     (upstream lrsnash)
//   Solver    = SolverT<I64, I128, GMP> automatic (lrsmp if GMP is absent)
//
// The exact algorithms are deterministic, so the result does not depend on
// which type finished.
//
// The API has two independent choices: the game class is the input type, and
// the method says what is returned.
//
//                      one equilibrium (Mode)      every extreme equilibrium
//   ConstantSumGame    solve(&g, &out, mode)       extreme_equilibria(&g)
//   GeneralSumGame     solve(&g, &out, mode, k)    extreme_equilibria(&g)
//
// solve(const ConstantSumGame*), by Mode:
//   Any (default) - one exact LP (zerosum.hpp). It returns an optimal
//           strategy pair, which for a constant-sum game is a Nash
//           equilibrium with the unique game value. In degenerate games with
//           several optimal strategies it can differ from the equilibrium
//           legacy lrsnash lists last.
//   Last / First - the lrsnash enumeration on (data, den - data). Last is
//           exactly the legacy solve_fast result.
// solve(const GeneralSumGame*), by Mode:
//   Any (default) - one equilibrium by Lemke-Howson (lemke.hpp). The caller
//           may choose the missing label; different labels can find
//           different equilibria. Fastest, but the equilibrium is in general
//           not the one legacy lrsnash returns.
//   Last  - legacy behaviour exactly: lrsnash enumeration order, the last
//           equilibrium found (nash.hpp). Exactly the legacy solve_full
//           result.
//   First - the first equilibrium in lrsnash order (early exit).
// extreme_equilibria: every extreme equilibrium, in upstream lrsnash's order,
//   with both players' payoffs. For a general-sum game this is upstream
//   lrsnash's output. (A degenerate game can have infinitely many
//   equilibria; each one is a convex combination of extreme ones.)
//
// Output of solve: row_strategy[0..rows-1], col_strategy[0..cols-1] (exactly
// that many entries are written) and value = player 1's payoff / den, as in
// legacy. Every method throws std::invalid_argument unless rows >= 1 and
// cols >= 1.
#pragma once

#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "../lib.h"
#include "arith.hpp"
#include "lemke.hpp"
#include "nash.hpp"
#include "zerosum.hpp"

namespace LRSNash {

// The game and output types. Until lib.h is switched over to this API they
// are the legacy lib.h structs under new names; then the definitions move
// here and the legacy names go.
//   ConstantSumGame{rows, cols, data, den}: player 1's payoff is data
//       (row-major rows x cols), player 2's is den - data.
//   GeneralSumGame{rows, cols, den, row_data, col_data}: player 1's and
//       player 2's payoffs (row-major); den only scales the reported payoffs.
//   FloatEquilibrium{row_strategy, col_strategy, value}: caller-owned
//       arrays of rows and cols floats, and value (see above).
using ConstantSumGame = FastInput;
using GeneralSumGame = Input;
using FloatEquilibrium = FloatOneSumOutput;

enum class Mode { Last, First, Any };

// One extreme equilibrium. Payoffs are divided by den, like
// FloatEquilibrium::value.
struct Equilibrium {
  std::vector<double> x, y;
  double payoff1 = 0, payoff2 = 0;
};

namespace detail {

// Visitor that keeps the last (or first) equilibrium as floats.
template <class Ar> struct FloatEq {
  float *x, *y;
  long rows, cols;
  bool first;
  double value = 0;
  bool on_eq2(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    for (long j = 0; j < cols; j++) y[j] = (float)(Ar::to_double(o[j + 1]) / d);
    value = Ar::to_double(o[cols + 1]) / d;
    return !first;  // First: stop at the first (x, y) pair
  }
  bool on_eq1(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    for (long i = 0; i < rows; i++) x[i] = (float)(Ar::to_double(o[i + 1]) / d);
    return !first;
  }
};

// Visitor that collects every equilibrium. For each player-1 vertex x,
// upstream reports the player-2 vertices y of its face first, then x.
template <class Ar> struct AllEq {
  std::vector<Equilibrium> *out;
  long rows, cols;
  double den;
  size_t pending = 0;  // first equilibrium of the current x
  bool on_eq2(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    Equilibrium e;
    e.y.resize(cols);
    for (long j = 0; j < cols; j++) e.y[j] = Ar::to_double(o[j + 1]) / d;
    e.payoff1 = Ar::to_double(o[cols + 1]) / d / den;
    out->push_back(std::move(e));
    return true;
  }
  bool on_eq1(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    std::vector<double> x(rows);
    for (long i = 0; i < rows; i++) x[i] = Ar::to_double(o[i + 1]) / d;
    const double v = Ar::to_double(o[rows + 1]) / d / den;
    for (size_t k = pending; k < out->size(); k++) {
      (*out)[k].x = x;
      (*out)[k].payoff2 = v;
    }
    pending = out->size();
    return true;
  }
};

template <class Ar> struct Workspace {
  ZeroSumLP<Ar> lp;
  LemkeHowson<Ar> lh;
  Nash<Ar> en;
};

}  // namespace detail

template <class... Ars> class SolverT {
  static_assert(sizeof...(Ars) > 0, "need at least one arithmetic type");

 public:
  // One equilibrium of a constant-sum game.
  //   Mode::Any (default) - one exact LP: an optimal strategy pair.
  //   Mode::Last / First  - legacy enumeration on (data, den - data); Last is
  //                         exactly the legacy solve_fast result.
  void solve(const ConstantSumGame *g, FloatEquilibrium *out, Mode mode = Mode::Any) {
    const long m = g->rows, n = g->cols;
    check_dims(m, n, "solve: rows and cols must be >= 1");
    const int *a = g->data;
    const long den = g->den;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    if (mode != Mode::Any) {
      auto B = [a, n, den](long i, long j) { return den - (long)a[i * n + j]; };
      enumerate(m, n, A, B, mode, out);
      out->value = (float)(out_value / g->den);
      return;
    }
    std::get<0>(ws).lp.setup(m, n, A);
    const double v = lp_chain<0>(m, n, A, out);
    out->value = (float)(v / g->den);
  }

  // One equilibrium of a general-sum game.
  //   Mode::Any (default) - one equilibrium by Lemke-Howson, dropping
  //                         missing_label (0-based): 0..rows-1 is row
  //                         strategy i = label, rows..rows+cols-1 is column
  //                         strategy j = label - rows. Ignored by Last/First.
  //   Mode::Last / First  - legacy enumeration; Last is exactly the legacy
  //                         solve_full result.
  void solve(const GeneralSumGame *g, FloatEquilibrium *out, Mode mode = Mode::Any,
             long missing_label = 0) {
    const long m = g->rows, n = g->cols;
    check_dims(m, n, "solve: rows and cols must be >= 1");
    const int *a = g->row_data, *b = g->col_data;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [b, n](long i, long j) { return (long)b[i * n + j]; };
    double v = 0;
    if (mode == Mode::Any) {
      if (missing_label < 0 || missing_label >= m + n)
        throw std::invalid_argument("solve: missing_label must be in [0, rows + cols)");
      attempt([&](auto &w) {
        w.lh.setup(m, n, A, B);
        w.lh.run(missing_label + 1);  // lemke.hpp labels are 1-based
        v = w.lh.solve_float(out->row_strategy, out->col_strategy);
      });
    } else {
      enumerate(m, n, A, B, mode, out);
      v = out_value;
    }
    out->value = (float)(v / g->den);
  }

  // Every extreme equilibrium of a constant-sum game: the enumeration on
  // (data, den - data).
  std::vector<Equilibrium> extreme_equilibria(const ConstantSumGame *g) {
    const long m = g->rows, n = g->cols;
    check_dims(m, n, "extreme_equilibria: rows and cols must be >= 1");
    const int *a = g->data;
    const long den = g->den;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [a, n, den](long i, long j) { return den - (long)a[i * n + j]; };
    return collect_all(m, n, A, B, den);
  }

  // Every extreme equilibrium of a general-sum game.
  std::vector<Equilibrium> extreme_equilibria(const GeneralSumGame *g) {
    const long m = g->rows, n = g->cols;
    check_dims(m, n, "extreme_equilibria: rows and cols must be >= 1");
    const int *a = g->row_data, *b = g->col_data;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [b, n](long i, long j) { return (long)b[i * n + j]; };
    return collect_all(m, n, A, B, g->den);
  }

 private:
  // The solvers need a real game. (nash.hpp's lrs systems also rely on it:
  // with rows, cols >= 1 the payoff rows give them full column rank.)
  static void check_dims(long rows, long cols, const char *msg) {
    if (rows < 1 || cols < 1) throw std::invalid_argument(msg);
  }

  // LP in type I; before a pivot that might not fit, promote the dictionary
  // to type I+1 and continue there (no work is repeated). Types that cannot
  // be converted from (bignums) restart in the next type instead.
  template <size_t I, class F>
  double lp_chain(long m, long n, F A, FloatEquilibrium *out) {
    auto &w = std::get<I>(ws).lp;
    if constexpr (I + 1 == sizeof...(Ars)) {
      w.run(false);  // last type: overflow propagates
    } else {
      using Cur = std::tuple_element_t<I, std::tuple<Ars...>>;
      constexpr bool convertible = std::is_same_v<Cur, arith::I64> ||
                                   std::is_same_v<Cur, arith::I128>;
      bool done;
      try {
        done = w.run(convertible);
      } catch (overflow_error &) {
        done = false;
      }
      if (!done) {
        auto &nx = std::get<I + 1>(ws).lp;
        if constexpr (convertible)
          nx.promote_from(w);
        else
          nx.setup(m, n, A);
        return lp_chain<I + 1>(m, n, A, out);
      }
    }
    return w.solve_float(out->row_strategy, out->col_strategy);
  }

  // Mode::Last / Mode::First with the lrsnash enumeration
  template <class FA, class FB>
  void enumerate(long m, long n, FA A, FB B, Mode mode, FloatEquilibrium *out) {
    // scratch first: an overflow part-way must not leave a mix of
    // results from two precisions in the caller's arrays
    xs.assign(m, 0.f);
    ys.assign(n, 0.f);
    attempt([&](auto &w) {
      using Ar = typename std::decay_t<decltype(w.en)>::T_ar;
      detail::FloatEq<Ar> vis{xs.data(), ys.data(), m, n, mode == Mode::First};
      w.en.solve(m, n, A, B, vis);
      out_value = vis.value;
    });
    for (long i = 0; i < m; i++) out->row_strategy[i] = xs[i];
    for (long j = 0; j < n; j++) out->col_strategy[j] = ys[j];
  }

  // every extreme equilibrium, payoffs divided by den
  template <class FA, class FB>
  std::vector<Equilibrium> collect_all(long m, long n, FA A, FB B, long den) {
    std::vector<Equilibrium> all;
    attempt([&](auto &w) {
      using Ar = typename std::decay_t<decltype(w.en)>::T_ar;
      all.clear();  // drop a partial result from a type that overflowed
      detail::AllEq<Ar> vis{&all, m, n, (double)den};
      w.en.solve(m, n, A, B, vis);
    });
    return all;
  }

  // run f(workspace) with each arithmetic type in turn until one fits
  template <class F> void attempt(F &&f) { attempt_from<0>(f); }
  template <size_t I, class F> void attempt_from(F &f) {
    if constexpr (I + 1 == sizeof...(Ars)) {
      f(std::get<I>(ws));  // last type: overflow propagates
    } else {
      try {
        f(std::get<I>(ws));
        return;
      } catch (overflow_error &) {
      }
      attempt_from<I + 1>(f);
    }
  }

  std::tuple<detail::Workspace<Ars>...> ws;
  std::vector<float> xs, ys;
  double out_value = 0;
};

#ifdef LRSNASH_HAVE_GMP
using Solver = SolverT<arith::I64, arith::I128, arith::GMP>;
using SolverGMP = SolverT<arith::GMP>;
#else
using Solver = SolverT<arith::I64, arith::I128, arith::MP<>>;
#endif
using Solver64 = SolverT<arith::I64>;
using Solver128 = SolverT<arith::I128>;
using SolverMP = SolverT<arith::MP<>>;

}  // namespace LRSNash
