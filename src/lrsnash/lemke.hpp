// Exact Lemke-Howson for one Nash equilibrium of a bimatrix game.
//
// Reference: B. von Stengel, "Computing equilibria for two-person games",
// Handbook of Game Theory vol. 3 (2002), sections 3 and 6.
//
// Shift the payoffs so that A' = A - min(A) + 1 >= 1 and B' = B - min(B) + 1 >= 1.
// Then
//   P = {x >= 0 : B'^T x <= 1}   labels: x_i = 0 -> i, (B'^T x)_j = 1 -> m+j
//   Q = {y >= 0 : A' y <= 1}     labels: (A' y)_i = 1 -> i, y_j = 0 -> m+j
// are polytopes, and the completely labeled pairs other than (0,0) are the
// Nash equilibria (after normalizing x and y). Lemke-Howson starts at (0,0),
// drops label k0, and alternates pivots between the two tableaus: each
// leaving variable's label is duplicated, so the variable with that label in
// the other tableau enters. It stops when label k0 is picked up again.
//
// Each tableau stores every variable as a column, and column c holds the
// variable with label c, so "enter the complement" means "enter the same
// column in the other tableau". Pivots are fraction-free (Bareiss) with
// exact division (Ar::tableau_pivot). Degenerate games are handled by the
// lexicographic minimum-ratio test over the slack columns (= det * B^-1),
// which is the standard perturbation argument. It never ties, so the path is
// unique and terminates.
//
// The equilibrium found depends on k0 and is in general NOT the one legacy
// lrsnash reports last. For constant-sum games every equilibrium has the same
// value.
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "arith.hpp"

namespace LRSNash {
namespace detail {

template <class Ar> struct LemkeHowson {
  using T = typename Ar::T;

  Ar ar;
  long m = 0, n = 0, W = 0;  // W = 1 + m + n columns: rhs, labels 1..m+n
  long shiftA = 0, shiftB = 0;
  std::vector<T> TP, TQ;           // n x W and m x W
  std::vector<long> basP, basQ;    // label of the basic variable per row
  T detP, detQ;
  std::uint64_t boundP = ~(std::uint64_t)0, boundQ = ~(std::uint64_t)0;
  long pivots = 0;

  T &P(long r, long c) { return TP[r * W + c]; }
  T &Q(long r, long c) { return TQ[r * W + c]; }

  // A(i, j), B(i, j): payoffs of players 1 and 2 (0-based)
  template <class FA, class FB> void setup(long m_, long n_, FA A, FB B) {
    m = m_;
    n = n_;
    W = 1 + m + n;
    long amin = std::numeric_limits<long>::max(), bmin = amin;
    for (long i = 0; i < m; i++)
      for (long j = 0; j < n; j++) {
        amin = std::min(amin, (long)A(i, j));
        bmin = std::min(bmin, (long)B(i, j));
      }
    shiftA = amin - 1;
    shiftB = bmin - 1;
    TP.resize(n * W);
    TQ.resize(m * W);
    basP.resize(n);
    basQ.resize(m);
    for (auto &t : TP) ar.set(t, 0);
    for (auto &t : TQ) ar.set(t, 0);
    // P row j: sum_i B'_ij x_i + s_j = 1     (x_i: column i, s_j: column m+j)
    for (long j = 0; j < n; j++) {
      ar.set(P(j, 0), 1);
      for (long i = 0; i < m; i++) ar.set(P(j, 1 + i), (long)B(i, j) - shiftB);
      ar.set(P(j, 1 + m + j), 1);
      basP[j] = 1 + m + j;
    }
    // Q row i: r_i + sum_j A'_ij y_j = 1     (r_i: column i, y_j: column m+j)
    for (long i = 0; i < m; i++) {
      ar.set(Q(i, 0), 1);
      ar.set(Q(i, 1 + i), 1);
      for (long j = 0; j < n; j++) ar.set(Q(i, 1 + m + j), (long)A(i, j) - shiftA);
      basQ[i] = 1 + i;
    }
    ar.set(detP, 1);
    ar.set(detQ, 1);
    boundP = boundQ = ~(std::uint64_t)0;
    pivots = 0;
  }

  // Lexicographic minimum ratio row for entering column c. The ratio vector of
  // row r is (T[r][0], T[r][s0], ..., T[r][s0+R-1]) / T[r][c] over rows with
  // T[r][c] > 0, where columns s0.. are the initial slack columns.
  long lexmin_ratio(T *Tb, long R, long c, long s0) {
    long best = -1;
    for (long r = 0; r < R; r++) {
      const T &arc = Tb[r * W + c];
      if (!ar.positive(arc)) continue;
      if (best < 0) {
        best = r;
        continue;
      }
      const T &abc = Tb[best * W + c];
      for (long t = -1; t < R; t++) {
        const long col = t < 0 ? 0 : s0 + t;
        // sign(ratio_r - ratio_best) for this column
        long cmp = ar.comprod(Tb[r * W + col], abc, Tb[best * W + col], arc);
        if (cmp < 0) {
          best = r;
          break;
        }
        if (cmp > 0) break;
      }
    }
    return best;
  }

  // Follow the path that drops label k0 (1..m+n).
  void run(long k0) {
    bool inP = k0 <= m;  // x_k0 enters P, or y_(k0-m) enters Q
    long enter = k0;
    for (;;) {
      T *Tb = inP ? TP.data() : TQ.data();
      const long R = inP ? n : m;
      const long s0 = inP ? 1 + m : 1;
      long row = lexmin_ratio(Tb, R, enter, s0);
      if (row < 0) throw std::logic_error("Lemke-Howson: unbounded ray");
      long *bas = inP ? basP.data() : basQ.data();
      const long leave = bas[row];
      if (inP)
        ar.tableau_pivot(Tb, W, R, W, row, enter, detP, boundP);
      else
        ar.tableau_pivot(Tb, W, R, W, row, enter, detQ, boundQ);
      bas[row] = enter;
      ++pivots;
      if (leave == k0) return;
      enter = leave;
      inP = !inP;
    }
  }

  // Exact result: x_i = X[i]/SX, y_j = Y[j]/SY with SX = sum X, SY = sum Y;
  // payoffs in the shifted games are u1' = detQ/SY, u2' = detP/SX, and the
  // original payoffs are u' + shift.
  void solution(T *X, T *Y) {
    for (long i = 0; i < m; i++) ar.set(X[i], 0);
    for (long j = 0; j < n; j++) ar.set(Y[j], 0);
    for (long r = 0; r < n; r++)
      if (basP[r] <= m) ar.copy(X[basP[r] - 1], P(r, 0));
    for (long r = 0; r < m; r++)
      if (basQ[r] > m) ar.copy(Y[basQ[r] - m - 1], Q(r, 0));
  }

  // Float output; returns player 1's payoff in the original game.
  double solve_float(float *x, float *y, double *u2 = nullptr) {
    double sx = 0, sy = 0;
    for (long i = 0; i < m; i++) x[i] = 0.f;
    for (long j = 0; j < n; j++) y[j] = 0.f;
    for (long r = 0; r < n; r++)
      if (basP[r] <= m) sx += Ar::to_double(P(r, 0));
    for (long r = 0; r < m; r++)
      if (basQ[r] > m) sy += Ar::to_double(Q(r, 0));
    for (long r = 0; r < n; r++)
      if (basP[r] <= m) x[basP[r] - 1] = (float)(Ar::to_double(P(r, 0)) / sx);
    for (long r = 0; r < m; r++)
      if (basQ[r] > m) y[basQ[r] - m - 1] = (float)(Ar::to_double(Q(r, 0)) / sy);
    if (u2) *u2 = Ar::to_double(detP) / sx + (double)shiftB;
    return Ar::to_double(detQ) / sy + (double)shiftA;
  }
};

}  // namespace detail
}  // namespace LRSNash
