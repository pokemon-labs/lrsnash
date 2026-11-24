// Public API (src/lrsnash/solve.hpp) vs the legacy library, for all four
// arithmetic types: extreme_equilibria(GeneralSumGame*) must reproduce the
// full legacy output sequence, solve(GeneralSumGame*, Last) the last legacy
// pair, solve(GeneralSumGame*, First) the first legacy pair. On constant-sum
// games the ConstantSumGame overloads (extreme_equilibria, solve Last/First)
// must give exactly the general-sum results on (a, den - a). "legacy float
// output differs" counts games where legacy's own float output is wrong
// (mptoi truncates numerators that exceed long); it is informational, not a
// mismatch.
// Needs a legacy GMP build that records every lrs_nashoutput call into
// g_seq; create it with tools/instrument_legacy.py:
//   python3 tools/instrument_legacy.py /tmp/legacy && cd /tmp/legacy &&
//   for f in lib lrslib lrsgmp; do g++ -O2 -c -I. $f.cc -DGMP -w -fpermissive; done
//   g++ -O2 -std=c++17 -Isrc tools/test_api.cc /tmp/legacy/*.o -lgmp -o test_api
//   ./test_api trials maxdim seed
#include "lrsnash/solve.hpp"
#include <cmath>
#include <cstdio>
#include <random>
extern std::vector<std::vector<double>> g_seq;  // legacy: every lrs_nashoutput, in order
using namespace LRSNash;
static bool close(double a, double b) { return std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(a)); }
static bool closef(double a, double b) { return std::fabs(a - b) <= 2e-6 * std::max(1.0, std::fabs(a)); }
// rebuild upstream's output sequence from extreme_equilibria (whose payoffs
// are divided by den; upstream's are not)
static std::vector<std::vector<double>> as_seq(const std::vector<Equilibrium> &all, int den) {
  std::vector<std::vector<double>> s;
  for (size_t k = 0; k < all.size(); k++) {
    std::vector<double> r{2.0}; r.insert(r.end(), all[k].y.begin(), all[k].y.end()); r.push_back(all[k].payoff1 * den); s.push_back(r);
    if (k + 1 == all.size() || all[k + 1].x != all[k].x) { std::vector<double> q{1.0}; q.insert(q.end(), all[k].x.begin(), all[k].x.end()); q.push_back(all[k].payoff2 * den); s.push_back(q); }
  }
  return s;
}
static bool same(const std::vector<std::vector<double>> &x, const std::vector<std::vector<double>> &y) {
  if (x.size() != y.size()) return false;
  for (size_t i = 0; i < x.size(); i++) { if (x[i].size() != y[i].size()) return false; for (size_t j = 0; j < x[i].size(); j++) if (!close(x[i][j], y[i][j])) return false; }
  return true;
}
static bool same_eqs(const std::vector<Equilibrium> &p, const std::vector<Equilibrium> &q) {
  if (p.size() != q.size()) return false;
  for (size_t k = 0; k < p.size(); k++)
    if (p[k].x != q[k].x || p[k].y != q[k].y || p[k].payoff1 != q[k].payoff1 || p[k].payoff2 != q[k].payoff2) return false;
  return true;
}
template <class S> struct Check {
  const char *name; S s; long bad_all = 0, bad_last = 0, bad_first = 0, bad_cs = 0, legacy_float_diff = 0, ovf = 0, n = 0;
  void run(long m, long n_, std::vector<int> &a, std::vector<int> &b, int den, bool constant_sum,
           const std::vector<std::vector<double>> &L, const std::vector<float> &lrs, const std::vector<float> &lcs, float lval) {
    try {
      GeneralSumGame in{m, n_, den, a.data(), b.data()};
      auto all = s.extreme_equilibria(&in);
      if (!same(as_seq(all, den), L)) ++bad_all;
      std::vector<float> rs(m + 2, -7.f), cs(n_ + 2, -7.f); FloatEquilibrium out{rs.data(), cs.data(), 0.f};
      s.solve(&in, &out, Mode::Last);
      // vs the last pair of legacy's exact output sequence (the real gate)
      const std::vector<double> *yl = nullptr, *xl = nullptr;
      for (auto &r : L) { if (r[0] == 2.0) yl = &r; if (r[0] == 1.0) xl = &r; }
      bool ok = yl && xl && closef(out.value, (*yl)[n_ + 1] / den) && rs[m] == -7.f && cs[n_] == -7.f;
      for (long i = 0; ok && i < m; i++) ok &= closef(rs[i], (*xl)[i + 1]);
      for (long j = 0; ok && j < n_; j++) ok &= closef(cs[j], (*yl)[j + 1]);
      if (!ok) ++bad_last;
      // vs legacy's float output, which mptoi truncation corrupts when a
      // numerator exceeds long (informational, not a regression)
      ok = closef(out.value, lval);
      for (long i = 0; i < m; i++) ok &= closef(rs[i], lrs[i]);
      for (long j = 0; j < n_; j++) ok &= closef(cs[j], lcs[j]);
      if (!ok) ++legacy_float_diff;
      s.solve(&in, &out, Mode::First);   // first pair: first "2" record and first "1" record
      const std::vector<double> *y0 = nullptr, *x0 = nullptr;
      for (auto &r : L) { if (r[0] == 2.0 && !y0) y0 = &r; if (r[0] == 1.0 && !x0) x0 = &r; }
      ok = y0 && x0 && closef(out.value, (*y0)[n_ + 1] / den);
      for (long i = 0; ok && i < m; i++) ok &= closef(rs[i], (*x0)[i + 1]);
      for (long j = 0; ok && j < n_; j++) ok &= closef(cs[j], (*y0)[j + 1]);
      if (!ok) ++bad_first;
      if (constant_sum) {  // b == den - a: the ConstantSumGame overloads must agree exactly
        ConstantSumGame cg{m, n_, a.data(), den};
        ok = same_eqs(s.extreme_equilibria(&cg), all);
        for (Mode md : {Mode::Last, Mode::First}) {
          std::vector<float> r1(m), c1(n_), r2(m), c2(n_);
          FloatEquilibrium o1{r1.data(), c1.data(), 0.f}, o2{r2.data(), c2.data(), 0.f};
          s.solve(&in, &o1, md);
          s.solve(&cg, &o2, md);
          ok &= r1 == r2 && c1 == c2 && o1.value == o2.value;
        }
        if (!ok) ++bad_cs;
      }
      ++n;
    } catch (overflow_error &) { ++ovf; }
  }
  void report() { printf("  %-10s games=%ld  all-seq mismatches=%ld  last mismatches=%ld  first mismatches=%ld  constant-sum overload mismatches=%ld  overflow=%ld  (legacy float output differs: %ld)\n", name, n, bad_all, bad_last, bad_first, bad_cs, ovf, legacy_float_diff); }
};
int main(int argc, char **argv) {
  int trials = atoi(argv[1]), maxdim = atoi(argv[2]); std::mt19937 gen(atoi(argv[3]));
  Check<Solver> c0{"auto"}; Check<Solver64> c1{"int64"}; Check<Solver128> c2{"int128"}; Check<SolverMP> c3{"lrsmp"}; Check<SolverGMP> c4{"gmp"};
  for (int t = 0; t < trials; ++t) {
    long m = std::uniform_int_distribution<int>(1, maxdim)(gen), n = std::uniform_int_distribution<int>(1, maxdim)(gen);
    int vmax = (t % 3 == 0) ? 2 : (t % 3 == 1 ? 80 : 1000);
    std::uniform_int_distribution<int> dis(0, vmax);
    std::vector<int> a(m * n), b(m * n); for (auto &x : a) x = dis(gen); for (auto &x : b) x = dis(gen);
    if (t % 5 == 1 && m > 1) for (long j = 0; j < n; ++j) { a[n + j] = a[j]; b[n + j] = b[j]; }
    if (t % 7 == 2 && n > 1) for (long i = 0; i < m; ++i) { a[i * n + 1] = a[i * n]; b[i * n + 1] = b[i * n]; }
    if (t % 2) for (long k = 0; k < m * n; k++) b[k] = vmax - a[k];
    g_seq.clear();
    std::vector<float> rs(m + 2), cs(n + 2); FloatOneSumOutput out{rs.data(), cs.data(), 0.f};
    Input in{m, n, vmax, a.data(), b.data()}; solve_full(&in, &out);   // legacy lib.h API (GMP build)
    auto L = g_seq;
    const bool cs_game = t % 2;
    c0.run(m, n, a, b, vmax, cs_game, L, rs, cs, out.value); c1.run(m, n, a, b, vmax, cs_game, L, rs, cs, out.value);
    c2.run(m, n, a, b, vmax, cs_game, L, rs, cs, out.value); c3.run(m, n, a, b, vmax, cs_game, L, rs, cs, out.value);
    c4.run(m, n, a, b, vmax, cs_game, L, rs, cs, out.value);
  }
  printf("trials=%d maxdim=%d\n", trials, maxdim); c0.report(); c1.report(); c2.report(); c3.report(); c4.report();
}
