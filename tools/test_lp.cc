// Exactness tests for the zero-sum LP (src/lrsnash/zerosum.hpp):
//  * fast pivot kernels (I64, I128) give dictionaries bit-identical to the
//    reference lrs pivot (dense_pivot) and to GMP and lrsmp;
//  * every result is an exact equilibrium (GMP rational check);
//  * optional: the value matches the legacy MP library (built with
//    -DLRSNASH_TEST_LEGACY, pass any 4th arg).
// Build (repo root):
//   g++ -O2 -std=c++17 -Isrc tools/test_lp.cc -lgmp -o test_lp
// With the legacy check (L from tools/instrument_legacy.py L, MP objects):
//   g++ -O2 -std=c++17 -DLRSNASH_TEST_LEGACY -Isrc -IL tools/test_lp.cc L/legacy/mp_*.o -lgmp -o test_lp
//   ./test_lp trials maxdim seed [legacy]
// Tests for ZeroSumLP: (1) fast pivot kernels give bit-identical dictionaries
// to the reference lrs pivot, (2) the result is an exact equilibrium,
// (3) the value matches legacy solve_fast.
#include "lrsnash/zerosum.hpp"
#ifdef LRSNASH_TEST_LEGACY
#include "legacy/lib.h"
#endif
#include <gmp.h>
#include <cstdio>
#include <random>
#include <string>
using namespace LRSNash;
using namespace LRSNash::arith;
template <class Base> struct RefPivot : Base {
  void pivot(typename Base::T *A, long st, long rows, long cols, long r, long s, typename Base::T &det, std::uint64_t &) {
    dense_pivot(static_cast<Base &>(*this), A, st, rows, cols, r, s, det);
  }
};
static void to_mpz(mpz_t z, const int64_t &a) { mpz_set_si(z, a); }
static void to_mpz(mpz_t z, const __int128 &a) {
  unsigned __int128 u = a < 0 ? -(unsigned __int128)a : a;
  mpz_set_ui(z, (uint64_t)(u >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (uint64_t)u);
  if (a < 0) mpz_neg(z, z);
}
static void to_mpz(mpz_t z, const GMP::T &a) { mpz_set(z, a.z); }
static void to_mpz(mpz_t z, const MP<>::T &a) {  // base 10^9 digits, a[0] = +-length
  long la = a[0] > 0 ? a[0] : -a[0]; mpz_set_ui(z, 0);
  for (long i = la - 1; i >= 1; i--) { mpz_mul_ui(z, z, 1000000000UL); mpz_add_ui(z, z, (unsigned long)a[i]); }
  if (a[0] < 0) mpz_neg(z, z);
}
static std::string str(const mpz_t z) { char *c = mpz_get_str(nullptr, 10, z); std::string s(c); free(c); return s; }

struct Exact { std::vector<std::string> X, Y; std::string S, det; long pivots; bool ok = false; };

template <class Ar> Exact run(long m, long n, const std::vector<int> &A) {
  Exact ex;
  try {
    detail::ZeroSumLP<Ar> lp;
    lp.setup(m, n, [&](long i, long j) { return (long)A[i * n + j]; });
    lp.run();
    std::vector<typename Ar::T> X(m), Y(n); typename Ar::T S, det;
    lp.solution(X.data(), Y.data(), S, det);
    mpz_t z; mpz_init(z);
    for (auto &x : X) { to_mpz(z, x); ex.X.push_back(str(z)); }
    for (auto &y : Y) { to_mpz(z, y); ex.Y.push_back(str(z)); }
    to_mpz(z, S); ex.S = str(z); to_mpz(z, det); ex.det = str(z);
    mpz_clear(z);
    ex.pivots = lp.pivots; ex.ok = true;
  } catch (overflow_error &) {}
  return ex;
}
static bool same(const Exact &a, const Exact &b) { return a.X == b.X && a.Y == b.Y && a.S == b.S && a.det == b.det && a.pivots == b.pivots; }

// exact equilibrium check in the shifted game A' = A - shift
static bool verify(long m, long n, const std::vector<int> &A, const Exact &e, double *value) {
  long amin = 1L << 40; for (int a : A) amin = std::min(amin, (long)a);
  long shift = amin - 1;
  mpz_t S, det, t, acc; mpz_inits(S, det, t, acc, nullptr);
  mpz_set_str(S, e.S.c_str(), 10); mpz_set_str(det, e.det.c_str(), 10);
  bool ok = mpz_sgn(S) > 0 && mpz_sgn(det) > 0;
  std::vector<mpz_t> X(m), Y(n);
  mpz_set_ui(acc, 0);
  for (long i = 0; i < m; i++) { mpz_init_set_str(X[i], e.X[i].c_str(), 10); ok &= mpz_sgn(X[i]) >= 0; mpz_add(acc, acc, X[i]); }
  ok &= mpz_cmp(acc, S) == 0;
  mpz_set_ui(acc, 0);
  for (long j = 0; j < n; j++) { mpz_init_set_str(Y[j], e.Y[j].c_str(), 10); ok &= mpz_sgn(Y[j]) >= 0; mpz_add(acc, acc, Y[j]); }
  ok &= mpz_cmp(acc, S) == 0;
  // value(A') = det/S.  x^T A' >= det/S * 1 :  sum_i X_i A'_ij >= det  ;  A' y <= det/S: sum_j A'_ij Y_j <= det
  for (long j = 0; j < n; j++) { mpz_set_ui(acc, 0); for (long i = 0; i < m; i++) { mpz_mul_si(t, X[i], A[i * n + j] - shift); mpz_add(acc, acc, t); } ok &= mpz_cmp(acc, det) >= 0; }
  for (long i = 0; i < m; i++) { mpz_set_ui(acc, 0); for (long j = 0; j < n; j++) { mpz_mul_si(t, Y[j], A[i * n + j] - shift); mpz_add(acc, acc, t); } ok &= mpz_cmp(acc, det) <= 0; }
  *value = mpz_get_d(det) / mpz_get_d(S) + shift;
  for (auto &x : X) mpz_clear(x); for (auto &y : Y) mpz_clear(y);
  mpz_clears(S, det, t, acc, nullptr);
  return ok;
}

int main(int argc, char **argv) {
  int trials = atoi(argv[1]), maxdim = atoi(argv[2]); unsigned seed = atoi(argv[3]);
#ifdef LRSNASH_TEST_LEGACY
  bool legacy = argc > 4;
#else
  if (argc > 4) { fprintf(stderr, "built without -DLRSNASH_TEST_LEGACY\n"); return 1; }
#endif
  std::mt19937 gen(seed);
  long bad_kernel = 0, bad_eq = 0, bad_val = 0, ovf64 = 0, ovf128 = 0, tot_piv = 0;
  for (int t = 0; t < trials; ++t) {
    long m = std::uniform_int_distribution<int>(1, maxdim)(gen), n = std::uniform_int_distribution<int>(1, maxdim)(gen);
    int kind = t % 4;
    int vmax = kind == 0 ? 2 : kind == 1 ? 100 : kind == 2 ? 1000000 : 81;
    std::uniform_int_distribution<int> dis(kind == 2 ? -vmax : 0, vmax);
    std::vector<int> A(m * n); for (auto &a : A) a = dis(gen);
    if (t % 5 == 1 && m > 1) for (long j = 0; j < n; j++) A[n + j] = A[j];
    if (t % 7 == 3 && n > 1) for (long i = 0; i < m; i++) A[i * n + 1] = A[i * n];
    Exact g = run<GMP>(m, n, A), gr = run<RefPivot<GMP>>(m, n, A);
    Exact a = run<I64>(m, n, A), ar = run<RefPivot<I64>>(m, n, A);
    Exact b = run<I128>(m, n, A), br = run<RefPivot<I128>>(m, n, A);
    Exact c = run<MP<>>(m, n, A);
    if (!same(c, g)) { ++bad_kernel; printf("MP mismatch t=%d\n", t); }
    if (!same(g, gr)) { ++bad_kernel; printf("GMP ref mismatch t=%d\n", t); }
    if (a.ok) { if (!same(a, g)) { ++bad_kernel; printf("I64 mismatch t=%d\n", t); } } else ++ovf64;
    if (ar.ok && !same(ar, g)) { ++bad_kernel; printf("I64 ref mismatch t=%d\n", t); }
    if (a.ok != ar.ok) printf("note: I64 fast/ref overflow differ t=%d (fast %d ref %d)\n", t, a.ok, ar.ok);
    if (b.ok) { if (!same(b, g)) { ++bad_kernel; printf("I128 mismatch t=%d\n", t); } } else ++ovf128;
    if (br.ok && !same(br, g)) { ++bad_kernel; printf("I128 ref mismatch t=%d\n", t); }
    if (b.ok != br.ok) printf("note: I128 fast/ref overflow differ t=%d (fast %d ref %d)\n", t, b.ok, br.ok);
    double v;
    if (!verify(m, n, A, g, &v)) { ++bad_eq; printf("NOT AN EQUILIBRIUM t=%d m=%ld n=%ld\n", t, m, n); }
    tot_piv += g.pivots;
#ifdef LRSNASH_TEST_LEGACY
    if (legacy && kind != 2) {  // legacy MP lib; den = vmax
      std::vector<float> rs(m + 2), cs(n + 2);
      FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
      FastInput in{m, n, A.data(), vmax};
      solve_fast(&in, &out);
      double lv = out.value * vmax;
      if (std::fabs(lv - v) > 1e-4 * std::max(1.0, std::fabs(v))) { ++bad_val; printf("VALUE MISMATCH t=%d legacy %g lp %g\n", t, lv, v); }
    }
#else
    (void)v;
#endif
  }
  printf("trials=%d kernel_mismatch=%ld not_eq=%ld value_mismatch=%ld I64_overflow=%ld I128_overflow=%ld avg_pivots=%.2f\n", trials, bad_kernel, bad_eq, bad_val, ovf64, ovf128, (double)tot_piv / trials);
}
