// Exactness tests for Lemke-Howson (src/lrsnash/lemke.hpp), for every
// missing label k0: fast tableau kernels bit-identical to the reference and
// across I64/I128/lrsmp/GMP; every result an exact Nash equilibrium with
// exactly the reported payoffs (GMP rational check). Also checks the public
// API: solve_full(Mode::Any, label) runs Lemke-Howson with missing label
// label + 1 and rejects labels outside [0, rows + cols); every solve_*
// rejects rows < 1 or cols < 1.
//   g++ -O2 -std=c++17 -Isrc tools/test_lh.cc -lgmp -o test_lh
//   ./test_lh trials maxdim seed
#include "lrsnash/lemke.hpp"
#include "lrsnash/solve.hpp"
#include <gmp.h>
#include <cstdio>
#include <random>
#include <string>
using namespace LRSNash; using namespace LRSNash::arith;
template <class Base> struct RefTab : Base {
  void tableau_pivot(typename Base::T *A, long st, long rows, long cols, long r, long s, typename Base::T &det, std::uint64_t &) {
    dense_tableau_pivot(static_cast<Base &>(*this), A, st, rows, cols, r, s, det);
  }
};
static void to_mpz(mpz_t z, const int64_t &a) { mpz_set_si(z, a); }
static void to_mpz(mpz_t z, const __int128 &a) { unsigned __int128 u = a < 0 ? -(unsigned __int128)a : a; mpz_set_ui(z, (uint64_t)(u >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (uint64_t)u); if (a < 0) mpz_neg(z, z); }
static void to_mpz(mpz_t z, const GMP::T &a) { mpz_set(z, a.z); }
static void to_mpz(mpz_t z, const MP<>::T &a) {  // base 10^9 digits, a[0] = +-length
  long la = a[0] > 0 ? a[0] : -a[0]; mpz_set_ui(z, 0);
  for (long i = la - 1; i >= 1; i--) { mpz_mul_ui(z, z, 1000000000UL); mpz_add_ui(z, z, (unsigned long)a[i]); }
  if (a[0] < 0) mpz_neg(z, z);
}
static std::string str(const mpz_t z) { char *c = mpz_get_str(nullptr, 10, z); std::string s(c); free(c); return s; }
struct Ex { std::vector<std::string> X, Y; std::string dP, dQ; long piv = 0; bool ok = false; };
template <class Ar> Ex run(long m, long n, const std::vector<int> &A, const std::vector<int> &B, long k0) {
  Ex e;
  try {
    detail::LemkeHowson<Ar> lh;
    lh.setup(m, n, [&](long i, long j) { return A[i * n + j]; }, [&](long i, long j) { return B[i * n + j]; });
    lh.run(k0);
    std::vector<typename Ar::T> X(m), Y(n); lh.solution(X.data(), Y.data());
    mpz_t z; mpz_init(z);
    for (auto &x : X) { to_mpz(z, x); e.X.push_back(str(z)); }
    for (auto &y : Y) { to_mpz(z, y); e.Y.push_back(str(z)); }
    to_mpz(z, lh.detP); e.dP = str(z); to_mpz(z, lh.detQ); e.dQ = str(z); mpz_clear(z);
    e.piv = lh.pivots; e.ok = true;
  } catch (overflow_error &) {}
  return e;
}
static bool same(const Ex &a, const Ex &b) { return a.X == b.X && a.Y == b.Y && a.dP == b.dP && a.dQ == b.dQ && a.piv == b.piv; }
// exact NE check + payoff check:  u1 = detQ/SY + shiftA == x^T A y ;  u2 = detP/SX + shiftB == x^T B y
static bool verify(long m, long n, const std::vector<int> &A, const std::vector<int> &B, const Ex &e) {
  long amin = 1L << 40, bmin = 1L << 40; for (int a : A) amin = std::min(amin, (long)a); for (int b : B) bmin = std::min(bmin, (long)b);
  mpz_t SX, SY, t, acc, best, tot, lhs, rhs; mpz_inits(SX, SY, t, acc, best, tot, lhs, rhs, nullptr);
  std::vector<mpz_t> X(m), Y(n); bool ok = true;
  for (long i = 0; i < m; i++) { mpz_init_set_str(X[i], e.X[i].c_str(), 10); ok &= mpz_sgn(X[i]) >= 0; mpz_add(SX, SX, X[i]); }
  for (long j = 0; j < n; j++) { mpz_init_set_str(Y[j], e.Y[j].c_str(), 10); ok &= mpz_sgn(Y[j]) >= 0; mpz_add(SY, SY, Y[j]); }
  ok &= mpz_sgn(SX) > 0 && mpz_sgn(SY) > 0;
  // player 1: (A Y)_i * SX <= sum_i' X_i' (A Y)_i'   for all i
  std::vector<mpz_t> AY(m); mpz_set_ui(tot, 0);
  for (long i = 0; i < m; i++) { mpz_init(AY[i]); for (long j = 0; j < n; j++) { mpz_mul_si(t, Y[j], A[i * n + j]); mpz_add(AY[i], AY[i], t); } mpz_mul(t, X[i], AY[i]); mpz_add(tot, tot, t); }
  for (long i = 0; i < m; i++) { mpz_mul(t, AY[i], SX); ok &= mpz_cmp(t, tot) <= 0; }
  // u1 = tot/(SX*SY) must equal detQ/SY + shiftA  <=>  tot == SX*(detQ + shiftA*SY)
  mpz_set_str(t, e.dQ.c_str(), 10); mpz_set_si(acc, amin - 1); mpz_mul(acc, acc, SY); mpz_add(t, t, acc); mpz_mul(t, t, SX); ok &= mpz_cmp(t, tot) == 0;
  // player 2
  std::vector<mpz_t> XB(n); mpz_set_ui(tot, 0);
  for (long j = 0; j < n; j++) { mpz_init(XB[j]); for (long i = 0; i < m; i++) { mpz_mul_si(t, X[i], B[i * n + j]); mpz_add(XB[j], XB[j], t); } mpz_mul(t, Y[j], XB[j]); mpz_add(tot, tot, t); }
  for (long j = 0; j < n; j++) { mpz_mul(t, XB[j], SY); ok &= mpz_cmp(t, tot) <= 0; }
  mpz_set_str(t, e.dP.c_str(), 10); mpz_set_si(acc, bmin - 1); mpz_mul(acc, acc, SX); mpz_add(t, t, acc); mpz_mul(t, t, SY); ok &= mpz_cmp(t, tot) == 0;
  for (auto &v : X) mpz_clear(v); for (auto &v : Y) mpz_clear(v); for (auto &v : AY) mpz_clear(v); for (auto &v : XB) mpz_clear(v);
  mpz_clears(SX, SY, t, acc, best, tot, lhs, rhs, nullptr);
  return ok;
}
// public API: solve_full(Mode::Any, k0 - 1) must equal LemkeHowson.run(k0) exactly
static SolverGMP api;
static bool api_same(long m, long n, std::vector<int> &A, std::vector<int> &B, long k0) {
  std::vector<float> rs(m), cs(n), rr(m), cr(n);
  Input in{m, n, 7, A.data(), B.data()};
  FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
  api.solve_full(&in, &out, Mode::Any, k0 - 1);
  detail::LemkeHowson<GMP> lh;
  lh.setup(m, n, [&](long i, long j) { return A[i * n + j]; }, [&](long i, long j) { return B[i * n + j]; });
  lh.run(k0);
  const double v = lh.solve_float(rr.data(), cr.data());
  return rs == rr && cs == cr && out.value == (float)(v / 7);
}
static bool api_rejects(long m, long n, std::vector<int> &A, std::vector<int> &B) {
  std::vector<float> rs(m), cs(n);
  Input in{m, n, 1, A.data(), B.data()};
  FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
  for (long bad : {-1L, m + n}) {
    try { api.solve_full(&in, &out, Mode::Any, bad); return false; } catch (std::invalid_argument &) {}
  }
  return true;
}
// every solve_* rejects rows < 1 or cols < 1 with invalid_argument
static bool dims_rejected() {
  std::vector<int> A(4, 1), B(4, 1); std::vector<float> rs(4), cs(4);
  FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
  for (long m : {0L, 2L}) for (long n : {0L, 2L}) {
    if (m > 0 && n > 0) continue;
    int caught = 0;
    Input in{m, n, 1, A.data(), B.data()}; FastInput fin{m, n, A.data(), 1};
    for (Mode md : {Mode::Any, Mode::Last}) {
      try { api.solve_full(&in, &out, md); } catch (std::invalid_argument &) { ++caught; }
      try { api.solve_fast(&fin, &out, md); } catch (std::invalid_argument &) { ++caught; }
    }
    try { api.solve_all(m, n, A.data(), B.data()); } catch (std::invalid_argument &) { ++caught; }
    if (caught != 5) return false;
  }
  return true;
}
int main(int argc, char **argv) {
  int trials = atoi(argv[1]), maxdim = atoi(argv[2]); std::mt19937 gen(atoi(argv[3]));
  long bad_k = 0, bad_eq = 0, bad_api = 0, o64 = 0, o128 = 0, runs = 0, piv = 0, maxpiv = 0;
  if (!dims_rejected()) ++bad_api;
  for (int t = 0; t < trials; ++t) {
    long m = std::uniform_int_distribution<int>(1, maxdim)(gen), n = std::uniform_int_distribution<int>(1, maxdim)(gen);
    int kind = t % 4; int vmax = kind == 0 ? 2 : kind == 1 ? 100 : kind == 2 ? 1000000 : 81;
    std::uniform_int_distribution<int> dis(kind == 2 ? -vmax : 0, vmax);
    std::vector<int> A(m * n), B(m * n); for (auto &a : A) a = dis(gen); for (auto &b : B) b = dis(gen);
    if (t % 5 == 1 && m > 1) for (long j = 0; j < n; j++) { A[n + j] = A[j]; B[n + j] = B[j]; }
    if (t % 7 == 3 && n > 1) for (long i = 0; i < m; i++) { A[i * n + 1] = A[i * n]; B[i * n + 1] = B[i * n]; }
    if (t % 11 == 4) for (long k = 0; k < m * n; k++) B[k] = vmax - A[k];  // constant-sum
    if (!api_rejects(m, n, A, B)) ++bad_api;
    for (long k0 = 1; k0 <= m + n; k0++) {
      if (!api_same(m, n, A, B, k0)) ++bad_api;
      Ex g = run<GMP>(m, n, A, B, k0), gr = run<RefTab<GMP>>(m, n, A, B, k0);
      Ex a = run<I64>(m, n, A, B, k0), ar = run<RefTab<I64>>(m, n, A, B, k0);
      Ex b = run<I128>(m, n, A, B, k0), br = run<RefTab<I128>>(m, n, A, B, k0);
      Ex c = run<MP<>>(m, n, A, B, k0);
      if (!same(c, g)) ++bad_k;
      if (!same(g, gr)) ++bad_k;
      if (a.ok) { if (!same(a, g)) ++bad_k; } else ++o64;
      if (ar.ok && !same(ar, g)) ++bad_k;
      if (a.ok != ar.ok) printf("note: I64 overflow differs t=%d\n", t);
      if (b.ok) { if (!same(b, g)) ++bad_k; } else ++o128;
      if (br.ok && !same(br, g)) ++bad_k;
      if (b.ok != br.ok) printf("note: I128 overflow differs t=%d\n", t);
      if (!verify(m, n, A, B, g)) { if (++bad_eq < 5) printf("NOT AN EQUILIBRIUM t=%d m=%ld n=%ld k0=%ld\n", t, m, n, k0); }
      ++runs; piv += g.piv; maxpiv = std::max(maxpiv, g.piv);
    }
  }
  printf("games=%d runs(all k0)=%ld kernel_mismatch=%ld not_eq=%ld api_mismatch=%ld I64_overflow=%ld I128_overflow=%ld avg_pivots=%.2f max_pivots=%ld\n", trials, runs, bad_k, bad_eq, bad_api, o64, o128, (double)piv / runs, maxpiv);
}
