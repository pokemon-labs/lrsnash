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
//   ConstantSumGame    solve(&g, &out, mode)       extreme_equilibria(&g[, f])
//   GeneralSumGame     solve(&g, &out, mode, k)    extreme_equilibria(&g[, f])
//   RationalGame       solve(&g, &out, mode, k)    extreme_equilibria(&g[, f])
//
// RationalGame is a general-sum game with num/den payoffs. extreme_equilibria
// returns doubles, or with a callback f gives each equilibrium exactly
// (ExactEquilibrium).
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
// legacy. Every method throws std::invalid_argument unless rows >= 1,
// cols >= 1 and den >= 1 (every denominator >= 1 for a RationalGame).
#pragma once

#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "arith.hpp"
#include "lemke.hpp"
#include "nash.hpp"
#include "zerosum.hpp"

namespace LRSNash {

// The game and output types (the legacy FastInput, Input and
// FloatOneSumOutput, with the same fields). den must be >= 1. A zero-sum
// game is a constant-sum game with any den >= 1: the strategies are the
// same, and value is scaled by 1/den.

// Player 1's payoff is data (row-major rows x cols), player 2's is
// den - data.
struct ConstantSumGame {
  long rows;
  long cols;
  const int *data;
  int den;
};

// Player 1's and player 2's payoffs (row-major rows x cols). den only
// scales the reported payoffs.
struct GeneralSumGame {
  long rows;
  long cols;
  int den;
  const int *row_data;
  const int *col_data;
};

// A general-sum game with rational payoffs: player 1's payoff at (i, j) is
// row_num[k] / row_den[k] and player 2's is col_num[k] / col_den[k], with
// k = i * cols + j. Every denominator must be >= 1. The solvers scale each
// player's matrix by the lcm of its denominators, which changes no best
// response, and report payoffs and value in these true units. They throw
// overflow_error if an lcm does not fit in long, or if a scaled payoff a has
// |a| > 2^62 (one bit of headroom for the solvers' shifts and negations).
struct RationalGame {
  long rows;
  long cols;
  const long *row_num, *row_den;
  const long *col_num, *col_den;
};

// Caller-owned arrays of rows and cols floats, and value = player 1's
// payoff / den (true units for a RationalGame).
struct FloatEquilibrium {
  float *row_strategy;
  float *col_strategy;
  float value;
};

enum class Mode { Last, First, Any };

// One extreme equilibrium. Payoffs are divided by den, like
// FloatEquilibrium::value (true units for a RationalGame).
struct Equilibrium {
  std::vector<double> x, y;
  double payoff1 = 0, payoff2 = 0;
};

// One extreme equilibrium, exactly, in the arithmetic type Ar that computed
// it (Ar::T is int64_t, __int128, lrsmp's digit array or a GMP integer; use
// Ar::to_double, or the type itself). Integers over a common denominator:
//   x[1..rows] / x[0]   player 1's strategy,  x[rows+1] / x[0]  player 2's payoff
//   y[1..cols] / y[0]   player 2's strategy,  y[cols+1] / y[0]  player 1's payoff
// The denominators are positive. The payoffs are in the units of the
// matrices the solver used; divide player k's by scale_k for the value the
// double-valued API reports: den for a ConstantSumGame or GeneralSumGame,
// the lcm of player k's denominators for a RationalGame. The arrays are
// valid only during the callback.
template <class Ar> struct ExactEquilibrium {
  using arith = Ar;
  using T = typename Ar::T;
  long rows, cols;
  const T *x, *y;
  long scale1, scale2;
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
  double den1, den2;  // player 1's and player 2's payoff scale
  size_t pending = 0;  // first equilibrium of the current x
  bool on_eq2(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    Equilibrium e;
    e.y.resize(cols);
    for (long j = 0; j < cols; j++) e.y[j] = Ar::to_double(o[j + 1]) / d;
    e.payoff1 = Ar::to_double(o[cols + 1]) / d / den1;
    out->push_back(std::move(e));
    return true;
  }
  bool on_eq1(const typename Ar::T *o, long) {
    const double d = Ar::to_double(o[0]);
    std::vector<double> x(rows);
    for (long i = 0; i < rows; i++) x[i] = Ar::to_double(o[i + 1]) / d;
    const double v = Ar::to_double(o[rows + 1]) / d / den2;
    for (size_t k = pending; k < out->size(); k++) {
      (*out)[k].x = x;
      (*out)[k].payoff2 = v;
    }
    pending = out->size();
    return true;
  }
};

// Visitor that passes each (x, y) pair to f as computed. x is read from the
// enumeration's current player-1 vertex, which exists before its face.
template <class Ar, class F> struct ExactEq {
  const Nash<Ar> *en;
  long rows, cols, scale1, scale2;
  F *f;
  bool on_eq2(const typename Ar::T *y, long) {
    (*f)(ExactEquilibrium<Ar>{rows, cols, en->vertex1(), y, scale1, scale2});
    return true;
  }
  bool on_eq1(const typename Ar::T *, long) { return true; }
};

// Visitor that copies each (x, y) pair into buf (rows+2 then cols+2 values).
template <class Ar> struct BufferEq {
  const Nash<Ar> *en;
  long rows, cols;
  std::vector<typename Ar::T> *buf;
  bool on_eq2(const typename Ar::T *y, long) {
    const typename Ar::T *x = en->vertex1();
    buf->insert(buf->end(), x, x + rows + 2);
    buf->insert(buf->end(), y, y + cols + 2);
    return true;
  }
  bool on_eq1(const typename Ar::T *, long) { return true; }
};

template <class Ar> struct Workspace {
  ZeroSumLP<Ar> lp;
  LemkeHowson<Ar> lh;
  Nash<Ar> en;
  std::vector<typename Ar::T> buf;  // exact pairs of a type that may overflow
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
    require(m >= 1 && n >= 1, "solve: rows and cols must be >= 1");
    require(g->den >= 1, "solve: den must be >= 1");
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
    require(m >= 1 && n >= 1, "solve: rows and cols must be >= 1");
    require(g->den >= 1, "solve: den must be >= 1");
    const int *a = g->row_data, *b = g->col_data;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [b, n](long i, long j) { return (long)b[i * n + j]; };
    solve_general(m, n, A, B, out, mode, missing_label, g->den);
  }

  // One equilibrium of a game with rational payoffs, as for GeneralSumGame.
  void solve(const RationalGame *g, FloatEquilibrium *out, Mode mode = Mode::Any,
             long missing_label = 0) {
    const long m = g->rows, n = g->cols;
    require(m >= 1 && n >= 1, "solve: rows and cols must be >= 1");
    const long s1 = scale_rational(g->row_num, g->row_den, m * n, ra);
    scale_rational(g->col_num, g->col_den, m * n, rb);
    const long *a = ra.data(), *b = rb.data();
    auto A = [a, n](long i, long j) { return a[i * n + j]; };
    auto B = [b, n](long i, long j) { return b[i * n + j]; };
    solve_general(m, n, A, B, out, mode, missing_label, s1);
  }

  // Every extreme equilibrium of a constant-sum game: the enumeration on
  // (data, den - data).
  std::vector<Equilibrium> extreme_equilibria(const ConstantSumGame *g) {
    std::vector<Equilibrium> all;
    with_game(g, [&](long m, long n, auto A, auto B, long s1, long s2) {
      collect_all(m, n, A, B, s1, s2, all);
    });
    return all;
  }
  // Every extreme equilibrium of a general-sum game.
  std::vector<Equilibrium> extreme_equilibria(const GeneralSumGame *g) {
    std::vector<Equilibrium> all;
    with_game(g, [&](long m, long n, auto A, auto B, long s1, long s2) {
      collect_all(m, n, A, B, s1, s2, all);
    });
    return all;
  }
  // Every extreme equilibrium of a game with rational payoffs.
  std::vector<Equilibrium> extreme_equilibria(const RationalGame *g) {
    std::vector<Equilibrium> all;
    with_game(g, [&](long m, long n, auto A, auto B, long s1, long s2) {
      collect_all(m, n, A, B, s1, s2, all);
    });
    return all;
  }

  // Every extreme equilibrium, exactly: f(const ExactEquilibrium<Ar> &) is
  // called once per equilibrium, in the order of extreme_equilibria. Ar is
  // the arithmetic type that finished, so f is usually a generic lambda.
  // Types before the last one keep their results until they finish without
  // overflow, so f never sees a partial result. The last (or only) type
  // passes them on as computed: if it throws overflow_error, f has already
  // seen some.
  template <class G, class F> void extreme_equilibria(const G *g, F &&f) {
    with_game(g, [&](long m, long n, auto A, auto B, long s1, long s2) {
      visit_exact<0>(m, n, A, B, s1, s2, f);
    });
  }

 private:
  static void require(bool ok, const char *msg) {
    if (!ok) throw std::invalid_argument(msg);
  }

  // Check a game and call body(m, n, A, B, scale1, scale2) with its payoff
  // functions. The solvers need a real game (nash.hpp's lrs systems also
  // rely on it: with rows, cols >= 1 the payoff rows give them full column
  // rank), and the reported payoffs are divided by the scales.
  template <class Body> void with_game(const ConstantSumGame *g, Body &&body) {
    const long m = g->rows, n = g->cols;
    require(m >= 1 && n >= 1, "extreme_equilibria: rows and cols must be >= 1");
    require(g->den >= 1, "extreme_equilibria: den must be >= 1");
    const int *a = g->data;
    const long den = g->den;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [a, n, den](long i, long j) { return den - (long)a[i * n + j]; };
    body(m, n, A, B, den, den);
  }
  template <class Body> void with_game(const GeneralSumGame *g, Body &&body) {
    const long m = g->rows, n = g->cols;
    require(m >= 1 && n >= 1, "extreme_equilibria: rows and cols must be >= 1");
    require(g->den >= 1, "extreme_equilibria: den must be >= 1");
    const int *a = g->row_data, *b = g->col_data;
    auto A = [a, n](long i, long j) { return (long)a[i * n + j]; };
    auto B = [b, n](long i, long j) { return (long)b[i * n + j]; };
    body(m, n, A, B, (long)g->den, (long)g->den);
  }
  template <class Body> void with_game(const RationalGame *g, Body &&body) {
    const long m = g->rows, n = g->cols;
    require(m >= 1 && n >= 1, "extreme_equilibria: rows and cols must be >= 1");
    const long s1 = scale_rational(g->row_num, g->row_den, m * n, ra);
    const long s2 = scale_rational(g->col_num, g->col_den, m * n, rb);
    const long *a = ra.data(), *b = rb.data();
    auto A = [a, n](long i, long j) { return a[i * n + j]; };
    auto B = [b, n](long i, long j) { return b[i * n + j]; };
    body(m, n, A, B, s1, s2);
  }

  // out[k] = num[k] * (L / den[k]) with L the lcm of den[0..count-1];
  // returns L. See RationalGame for the limits.
  static long scale_rational(const long *num, const long *den, long count,
                             std::vector<long> &out) {
    long L = 1;
    for (long k = 0; k < count; k++) {
      require(den[k] >= 1, "RationalGame: denominators must be >= 1");
      const long g = std::gcd(L, den[k]);
      if (__builtin_mul_overflow(L / g, den[k], &L))
        throw overflow_error("RationalGame: lcm of the denominators");
    }
    out.resize(count);
    constexpr long lim = 1L << 62;
    for (long k = 0; k < count; k++) {
      long a;
      if (__builtin_mul_overflow(num[k], L / den[k], &a) || a > lim || a < -lim)
        throw overflow_error("RationalGame: scaled payoff exceeds 2^62");
      out[k] = a;
    }
    return L;
  }

  // solve for a general-sum game; value is player 1's payoff / scale1
  template <class FA, class FB>
  void solve_general(long m, long n, FA A, FB B, FloatEquilibrium *out, Mode mode,
                     long missing_label, long scale1) {
    double v = 0;
    if (mode == Mode::Any) {
      require(missing_label >= 0 && missing_label < m + n,
              "solve: missing_label must be in [0, rows + cols)");
      attempt([&](auto &w) {
        w.lh.setup(m, n, A, B);
        w.lh.run(missing_label + 1);  // lemke.hpp labels are 1-based
        v = w.lh.solve_float(out->row_strategy, out->col_strategy);
      });
    } else {
      enumerate(m, n, A, B, mode, out);
      v = out_value;
    }
    out->value = (float)(v / scale1);
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

  // every extreme equilibrium, player k's payoffs divided by scale_k
  template <class FA, class FB>
  void collect_all(long m, long n, FA A, FB B, long scale1, long scale2,
                   std::vector<Equilibrium> &all) {
    attempt([&](auto &w) {
      using Ar = typename std::decay_t<decltype(w.en)>::T_ar;
      all.clear();  // drop a partial result from a type that overflowed
      detail::AllEq<Ar> vis{&all, m, n, (double)scale1, (double)scale2};
      w.en.solve(m, n, A, B, vis);
    });
  }

  // extreme_equilibria with an exact visitor, from type I on
  template <size_t I, class FA, class FB, class F>
  void visit_exact(long m, long n, FA A, FB B, long scale1, long scale2, F &f) {
    auto &w = std::get<I>(ws);
    using Ar = typename std::decay_t<decltype(w.en)>::T_ar;
    if constexpr (I + 1 == sizeof...(Ars)) {
      detail::ExactEq<Ar, F> vis{&w.en, m, n, scale1, scale2, &f};
      w.en.solve(m, n, A, B, vis);  // last type: overflow propagates
    } else {
      w.buf.clear();
      detail::BufferEq<Ar> vis{&w.en, m, n, &w.buf};
      try {
        w.en.solve(m, n, A, B, vis);
      } catch (overflow_error &) {
        visit_exact<I + 1>(m, n, A, B, scale1, scale2, f);
        return;
      }
      const std::size_t len = (std::size_t)(m + n + 4);
      for (std::size_t k = 0; k < w.buf.size(); k += len)
        f(ExactEquilibrium<Ar>{m, n, &w.buf[k], &w.buf[k + m + 2], scale1, scale2});
    }
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
  std::vector<long> ra, rb;  // a RationalGame's scaled payoffs
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
