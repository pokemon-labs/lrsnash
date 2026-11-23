// Exact solution of a two-player constant-sum (matrix) game by one LP.
//
// For a constant-sum game every Nash equilibrium is a pair of optimal
// (maximin / minimax) strategies, the value is unique, and the equilibria are
// interchangeable. So a single equilibrium can be read off one optimal simplex
// dictionary, with no vertex enumeration. The legacy lrsnash path instead
// enumerates every vertex of player 1's best-response polytope (about 480 for
// a 9x9 game and 48,000 for 15x15) and runs a player-2 subproblem for each
// one.
//
// Formulation (row player maximizes x^T A y). Shift A' = A - min(A) + 1, so
// every A'_ij >= 1, and solve
//     max 1^T y'   s.t.   A' y' <= 1,  y' >= 0.
// y' = 0 is feasible, so no phase 1 is needed, and A' >= 1 makes the LP
// bounded. At an optimum with value z* = 1^T y' > 0:
//     y = y' / z*,   x = x' / z*  (x' = dual solution = reduced costs of the
//     slacks),   value(A') = 1 / z*,   value(A) = 1/z* + min(A) - 1.
// Strong duality gives 1^T x' = z*, so x and y are probability vectors, and
// (x, y) is an equilibrium of A.
//
// Arithmetic is exact. The dictionary is lrs's integer (fraction-free)
// dictionary, with row i meaning basic_i = (D[i][0] + sum_j D[i][j] cobasic_j)
// / det, and pivots use the lrs integer pivot (Ar::pivot). Entering rule:
// Dantzig (largest objective coefficient) until the first degenerate pivot,
// then Bland's smallest-index rule, which cannot cycle. Leaving rule: minimum
// ratio, ties broken by smallest basic variable index.
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "arith.hpp"

namespace LRSNash {
namespace detail {

template <class Ar> struct ZeroSumLP {
  using T = typename Ar::T;

  Ar ar;
  long m = 0, n = 0, stride = 0;
  long shift = 0;     // A = A' + shift
  std::vector<T> D;   // (m+1) x (n+1); row 0 = objective, column 0 = rhs
  std::vector<long> B;  // B[i], i = 1..m: variable basic in row i
  std::vector<long> N;  // N[j], j = 1..n: variable cobasic in column j
  T det;
  std::uint64_t bound = ~(std::uint64_t)0;  // see Ar::pivot
  long pivots = 0;
  bool bland = false;  // switched to Bland's rule
  // Variables: 1..n are y'_1..y'_n; n+1..n+m are the slacks of rows 1..m.

  T &at(long i, long j) { return D[i * stride + j]; }

  // payoff(i, j) returns the row player's payoff A[i][j] (0-based) as long.
  template <class F> void setup(long m_, long n_, F payoff) {
    m = m_;
    n = n_;
    stride = n + 1;
    D.resize((m + 1) * stride);
    B.resize(m + 1);
    N.resize(n + 1);
    long amin = std::numeric_limits<long>::max();
    for (long i = 0; i < m; i++)
      for (long j = 0; j < n; j++) amin = std::min(amin, (long)payoff(i, j));
    shift = amin - 1;
    ar.set(at(0, 0), 0);
    for (long j = 1; j <= n; j++) ar.set(at(0, j), 1);
    for (long i = 1; i <= m; i++) {
      ar.set(at(i, 0), 1);
      for (long j = 1; j <= n; j++)
        ar.set(at(i, j), -((long)payoff(i - 1, j - 1) - shift));
    }
    for (long i = 1; i <= m; i++) B[i] = n + i;
    for (long j = 1; j <= n; j++) N[j] = j;
    ar.set(det, 1);
    bound = Ar::magnitude_bound(D.data(), stride, m + 1, n + 1);
    pivots = 0;
    bland = false;
  }

  // Pivot to optimality. With stop_on_risk, stop (returning false) before a
  // pivot whose result is not provably representable, leaving the dictionary
  // intact so that it can be promoted to a wider type (promote_from) and
  // continued. Otherwise overflow throws overflow_error.
  bool run(bool stop_on_risk = false) {
    for (;;) {
      long s = 0;
      for (long j = 1; j <= n; j++) {
        if (!ar.positive(at(0, j))) continue;
        if (s == 0 || (bland ? N[j] < N[s] : ar.greater(at(0, j), at(0, s))))
          s = j;
      }
      if (s == 0) return true;  // optimal
      long r = 0;
      for (long i = 1; i <= m; i++) {
        if (!ar.negative(at(i, s))) continue;
        if (r == 0) {
          r = i;
          continue;
        }
        // ratio_i < ratio_r  <=>  D[i][0]*D[r][s] - D[r][0]*D[i][s] > 0
        long c = ar.comprod(at(i, 0), at(r, s), at(r, 0), at(i, s));
        if (c > 0 || (c == 0 && B[i] < B[r])) r = i;
      }
      if (r == 0) throw std::logic_error("zero-sum LP unbounded");
      if (stop_on_risk &&
          !Ar::pivot_safe(D.data(), stride, n + 1, r, s, det, bound))
        return false;
      if (ar.zero(at(r, 0))) bland = true;
      ar.pivot(D.data(), stride, m + 1, n + 1, r, s, det, bound);
      std::swap(B[r], N[s]);
      ++pivots;
    }
  }

  // Continue the solve of a narrower LP (I64 -> I128 -> bignum) after its
  // run(true) returned false. Values are converted exactly, so the rest of
  // the pivot sequence, and the result, is the same as for a solve done
  // entirely in Ar.
  template <class From> void promote_from(const ZeroSumLP<From> &o) {
    m = o.m;
    n = o.n;
    stride = o.stride;
    shift = o.shift;
    B = o.B;
    N = o.N;
    pivots = o.pivots;
    bland = o.bland;
    D.resize(o.D.size());
    for (size_t i = 0; i < D.size(); i++) convert(D[i], o.D[i]);
    convert(det, o.det);
    bound = Ar::magnitude_bound(D.data(), stride, m + 1, n + 1);
  }
  static void convert(T &dst, std::int64_t v) { Ar::from_i64(dst, v); }
  static void convert(T &dst, __int128 v) { Ar::from_i128(dst, v); }

  // Exact solution: x_i = X[i]/S, y_j = Y[j]/S, value(A) = det/S + shift.
  // X, Y must have room for m and n entries.
  void solution(T *X, T *Y, T &S, T &Det) {
    for (long i = 0; i < m; i++) ar.set(X[i], 0);
    for (long j = 0; j < n; j++) ar.set(Y[j], 0);
    for (long i = 1; i <= m; i++)
      if (B[i] <= n) ar.copy(Y[B[i] - 1], at(i, 0));
    for (long j = 1; j <= n; j++)
      if (N[j] > n) {
        ar.copy(X[N[j] - n - 1], at(0, j));
        ar.changesign(X[N[j] - n - 1]);
      }
    ar.copy(S, at(0, 0));
    ar.copy(Det, det);
  }

  // Float output: x (m entries), y (n entries), value of A (row player).
  double solve_float(float *x, float *y) {
    const double S = Ar::to_double(at(0, 0));
    for (long i = 0; i < m; i++) x[i] = 0.f;
    for (long j = 0; j < n; j++) y[j] = 0.f;
    for (long i = 1; i <= m; i++)
      if (B[i] <= n) y[B[i] - 1] = (float)(Ar::to_double(at(i, 0)) / S);
    for (long j = 1; j <= n; j++)
      if (N[j] > n) x[N[j] - n - 1] = (float)(-Ar::to_double(at(0, j)) / S);
    return Ar::to_double(det) / S + (double)shift;
  }
};

}  // namespace detail
}  // namespace LRSNash
