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
// drops label k0, and alternates pivots between the two polytopes: each
// leaving variable's label is duplicated, so the variable with that label in
// the other polytope enters. It stops when label k0 is picked up again.
//
// Variables are named by their labels. In P, label i is x_i and label m+j is
// the slack s_j of constraint j (sum_i B'_ij x_i + s_j = 1); in Q, label i is
// the slack r_i of constraint i (r_i + sum_j A'_ij y_j = 1) and label m+j is
// y_j. A variable enters only where it is cobasic.
//
// Each polytope is an lrs integer dictionary that stores only the cobasic
// columns. Row k belongs to constraint k and reads
//   basic_k = (D[k][0] + sum_c D[k][c] * cobasic_c) / det,
// and a pivot is the lrs fraction-free pivot with exact division (Ar::pivot;
// the pivot element is always negative here, so the kernel negates column s,
// stores -det at (r, s) and keeps row r). In terms of the full tableau
// (rows sum_c T[k][c] var_c = T[k][0], basic coefficient det), the
// dictionary holds T[k][0] = D[k][0] and T[k][c] = -D[k][column of c] for a
// cobasic c; a basic variable's tableau column is det times a unit vector
// and is not stored. This does the same arithmetic as pivoting the full
// tableau, without the m (or n) basic columns that stay unit vectors.
//
// Degenerate games are handled by the lexicographic minimum-ratio test over
// the tableau columns of the initial slacks (= det * B^-1), the standard
// perturbation argument. It never ties, so the path is unique and
// terminates. The tableau entries it compares are read from the dictionary
// as above (det or 0 for a basic slack).
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

  // One polytope's dictionary: `rows` constraints, columns 0 (right-hand
  // side) and 1..stride-1 (the cobasic variables).
  struct Side {
    std::vector<T> D;
    long rows = 0, stride = 0;
    std::vector<long> bas;  // bas[k]: label of the basic variable of row k
    std::vector<long> col;  // col[label]: its column if cobasic, else 0
    std::vector<long> row;  // row[label]: its row if basic, else -1
    long lex0 = 0, lexn = 0;  // the initial slacks: labels lex0..lex0+lexn-1
    T det;
    std::uint64_t bound = ~(std::uint64_t)0;  // see Ar::pivot
    T &d(long k, long c) { return D[k * stride + c]; }
  };

  Ar ar;
  long m = 0, n = 0;
  long shiftA = 0, shiftB = 0;
  Side P, Q;
  T detP, detQ;  // the determinants after run()
  long pivots = 0;

  // A(i, j), B(i, j): payoffs of players 1 and 2 (0-based)
  template <class FA, class FB> void setup(long m_, long n_, FA A, FB B) {
    m = m_;
    n = n_;
    long amin = std::numeric_limits<long>::max(), bmin = amin;
    for (long i = 0; i < m; i++)
      for (long j = 0; j < n; j++) {
        amin = std::min(amin, (long)A(i, j));
        bmin = std::min(bmin, (long)B(i, j));
      }
    shiftA = amin - 1;
    shiftB = bmin - 1;
    // P: row j is s_j = 1 - sum_i B'_ij x_i; x_i (label i) in column i
    init_side(P, n, m, m + 1, n);
    for (long j = 0; j < n; j++) {
      ar.set(P.d(j, 0), 1);
      for (long i = 0; i < m; i++) ar.set(P.d(j, 1 + i), -((long)B(i, j) - shiftB));
      basic(P, j, m + 1 + j);
    }
    for (long i = 0; i < m; i++) cobasic(P, 1 + i, 1 + i);
    // Q: row i is r_i = 1 - sum_j A'_ij y_j; y_j (label m+1+j) in column j+1
    init_side(Q, m, n, 1, m);
    for (long i = 0; i < m; i++) {
      ar.set(Q.d(i, 0), 1);
      for (long j = 0; j < n; j++) ar.set(Q.d(i, 1 + j), -((long)A(i, j) - shiftA));
      basic(Q, i, 1 + i);
    }
    for (long j = 0; j < n; j++) cobasic(Q, m + 1 + j, 1 + j);
    // unknown: each side's first pivot checks its rows one by one, which
    // costs about what a magnitude_bound pass would
    P.bound = Q.bound = ~(std::uint64_t)0;
    pivots = 0;
  }

  // Follow the path that drops label k0 (1..m+n).
  void run(long k0) {
    bool inP = k0 <= m;  // x_k0 enters P, or y_(k0-m) enters Q
    long enter = k0;
    for (;;) {
      Side &S = inP ? P : Q;
      const long r = lexmin_ratio(S, enter);
      if (r < 0) throw std::logic_error("Lemke-Howson: unbounded ray");
      const long leave = S.bas[r], c = S.col[enter];
      ar.pivot(S.D.data(), S.stride, S.rows, S.stride, r, c, S.det, S.bound);
      basic(S, r, enter);
      cobasic(S, leave, c);
      ++pivots;
      if (leave == k0) break;
      enter = leave;
      inP = !inP;
    }
    ar.copy(detP, P.det);
    ar.copy(detQ, Q.det);
  }

  // Exact result: x_i = X[i]/SX, y_j = Y[j]/SY with SX = sum X, SY = sum Y;
  // payoffs in the shifted games are u1' = detQ/SY, u2' = detP/SX, and the
  // original payoffs are u' + shift.
  void solution(T *X, T *Y) {
    for (long i = 0; i < m; i++) ar.set(X[i], 0);
    for (long j = 0; j < n; j++) ar.set(Y[j], 0);
    for (long r = 0; r < n; r++)
      if (P.bas[r] <= m) ar.copy(X[P.bas[r] - 1], P.d(r, 0));
    for (long r = 0; r < m; r++)
      if (Q.bas[r] > m) ar.copy(Y[Q.bas[r] - m - 1], Q.d(r, 0));
  }

  // Float output; returns player 1's payoff in the original game.
  double solve_float(float *x, float *y, double *u2 = nullptr) {
    double sx = 0, sy = 0;
    for (long i = 0; i < m; i++) x[i] = 0.f;
    for (long j = 0; j < n; j++) y[j] = 0.f;
    for (long r = 0; r < n; r++)
      if (P.bas[r] <= m) sx += Ar::to_double(P.d(r, 0));
    for (long r = 0; r < m; r++)
      if (Q.bas[r] > m) sy += Ar::to_double(Q.d(r, 0));
    for (long r = 0; r < n; r++)
      if (P.bas[r] <= m) x[P.bas[r] - 1] = (float)(Ar::to_double(P.d(r, 0)) / sx);
    for (long r = 0; r < m; r++)
      if (Q.bas[r] > m) y[Q.bas[r] - m - 1] = (float)(Ar::to_double(Q.d(r, 0)) / sy);
    if (u2) *u2 = Ar::to_double(detP) / sx + (double)shiftB;
    return Ar::to_double(detQ) / sy + (double)shiftA;
  }

 private:
  void init_side(Side &S, long rows, long ncob, long lex0, long lexn) {
    S.rows = rows;
    S.stride = 1 + ncob;
    S.D.resize(rows * S.stride);
    S.bas.resize(rows);
    S.col.resize(m + n + 1);  // every label 1..m+n is set by basic/cobasic
    S.row.resize(m + n + 1);
    S.lex0 = lex0;
    S.lexn = lexn;
    ar.set(S.det, 1);
  }
  static void basic(Side &S, long r, long label) {
    S.bas[r] = label;
    S.row[label] = r;
    S.col[label] = 0;
  }
  static void cobasic(Side &S, long label, long c) {
    S.col[label] = c;
    S.row[label] = -1;
  }

  // Leaving row for the entering (cobasic) label e: the lexicographic minimum
  // of (T[k][0], T[k][s_1], ..., T[k][s_R]) / T[k][e] over the rows with
  // T[k][e] > 0, where s_1..s_R are the initial slacks. Rows k, b compare by
  // the sign of ratio_k - ratio_b column by column, which is
  //   sign(T[b][t] D[k][ce] - T[k][t] D[b][ce])      (ce = column of e)
  // since T[.][e] = -D[.][ce] > 0. With T[.][0] = D[.][0], T[.][t] =
  // -D[.][ct] for a cobasic slack t, and T[.][t] = det [. == row of t] for a
  // basic one (det > 0, D[.][ce] < 0) this gives the cases below.
  long lexmin_ratio(Side &S, long e) {
    const long ce = S.col[e];
    long best = -1;
    for (long k = 0; k < S.rows; k++) {
      const T &dke = S.d(k, ce);
      if (!ar.negative(dke)) continue;
      if (best < 0) {
        best = k;
        continue;
      }
      const T &dbe = S.d(best, ce);
      long cmp = ar.comprod(S.d(best, 0), dke, S.d(k, 0), dbe);
      for (long t = S.lex0; cmp == 0 && t < S.lex0 + S.lexn; t++) {
        if (const long ct = S.col[t])
          cmp = ar.comprod(S.d(k, ct), dbe, S.d(best, ct), dke);
        else
          cmp = S.row[t] == k ? 1 : S.row[t] == best ? -1 : 0;
      }
      if (cmp < 0) best = k;
    }
    return best;
  }
};

}  // namespace detail
}  // namespace LRSNash
