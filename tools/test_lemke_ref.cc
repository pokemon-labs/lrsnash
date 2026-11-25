// Gate for internal changes to Lemke-Howson (src/lrsnash/lemke.hpp): for
// every missing label and all four arithmetic types, the current code must
// follow exactly the same path as a frozen reference copy: the same number of
// pivots, the same exact solution (X, Y, detP, detQ), the same float output,
// and the same overflow behaviour. One solver object per type is reused for
// games of changing shape.
//
// The reference is the pre-change lemke.hpp with namespace LRSNash renamed
// (its arith.hpp is the fff2737 copy that test_refactor uses):
//   REF=/path/to/ref; mkdir -p $REF/lrsnash_ref
//   for f in arith lrs nash; do git show fff2737:src/lrsnash/$f.hpp |
//     sed 's/\bLRSNash\b/LRSNashRef/g' > $REF/lrsnash_ref/$f.hpp; done
//   git show 18834e9:src/lrsnash/lemke.hpp |
//     sed 's/\bLRSNash\b/LRSNashRef/g' > $REF/lrsnash_ref/lemke.hpp
//   g++ -O2 -std=c++17 -Isrc -I$REF tools/test_lemke_ref.cc -lgmp -o test_lemke_ref
//   ./test_lemke_ref trials maxdim seed
#include "lrsnash/lemke.hpp"
#include "lrsnash_ref/lemke.hpp"
#include <gmp.h>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

static void to_mpz(mpz_t z, const int64_t &a) { mpz_set_si(z, a); }
static void to_mpz(mpz_t z, const __int128 &a) {
  unsigned __int128 u = a < 0 ? -(unsigned __int128)a : a;
  mpz_set_ui(z, (uint64_t)(u >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (uint64_t)u);
  if (a < 0) mpz_neg(z, z);
}
template <class G> static auto to_mpz(mpz_t z, const G &a) -> decltype(mpz_set(z, a.z)) { mpz_set(z, a.z); }
template <class M> static auto to_mpz(mpz_t z, const M &a) -> decltype((void)a.v[0]) {  // lrsmp
  long la = a[0] > 0 ? a[0] : -a[0]; mpz_set_ui(z, 0);
  for (long i = la - 1; i >= 1; i--) { mpz_mul_ui(z, z, 1000000000UL); mpz_add_ui(z, z, (unsigned long)a[i]); }
  if (a[0] < 0) mpz_neg(z, z);
}
template <class T> static std::string str(const T &a) {
  mpz_t z; mpz_init(z); to_mpz(z, a);
  char *c = mpz_get_str(nullptr, 10, z); std::string s(c); free(c); mpz_clear(z);
  return s;
}

struct Res {
  std::string err, x, y, dp, dq;
  long piv = -1;
  std::vector<float> fx, fy;
  double v1 = 0, v2 = 0;
  bool operator==(const Res &o) const {
    return err == o.err && x == o.x && y == o.y && dp == o.dp && dq == o.dq && piv == o.piv && fx == o.fx &&
           fy == o.fy && v1 == o.v1 && v2 == o.v2;
  }
};

template <class LH> static Res run(LH &lh, long m, long n, const std::vector<int> &A, const std::vector<int> &B, long k0) {
  Res r;
  try {
    lh.setup(m, n, [&](long i, long j) { return (long)A[i * n + j]; }, [&](long i, long j) { return (long)B[i * n + j]; });
    lh.run(k0);
    using T = typename std::decay_t<decltype(lh.detP)>;
    std::vector<T> X(m), Y(n);
    lh.solution(X.data(), Y.data());
    for (auto &v : X) r.x += str(v) + ",";
    for (auto &v : Y) r.y += str(v) + ",";
    r.dp = str(lh.detP); r.dq = str(lh.detQ); r.piv = lh.pivots;
    r.fx.resize(m); r.fy.resize(n);
    r.v1 = lh.solve_float(r.fx.data(), r.fy.data(), &r.v2);
  } catch (std::overflow_error &) { r.err = "overflow"; }
  catch (std::exception &e) { r.err = e.what(); }
  return r;
}

template <class New, class Ref> struct Pair {
  const char *name;
  LRSNash::detail::LemkeHowson<New> a;
  LRSNashRef::detail::LemkeHowson<Ref> b;
  long runs = 0, bad = 0, ovf = 0, piv = 0;
  explicit Pair(const char *nm) : name(nm) {}
  void check(long m, long n, const std::vector<int> &A, const std::vector<int> &B, long k0) {
    Res x = run(a, m, n, A, B, k0), y = run(b, m, n, A, B, k0);
    ++runs;
    if (!(x == y) && ++bad <= 3)
      printf("  MISMATCH %s %ldx%ld k0=%ld: pivots %ld vs %ld, err '%s' vs '%s'\n", name, m, n, k0, x.piv, y.piv, x.err.c_str(), y.err.c_str());
    if (!x.err.empty()) ++ovf;
    piv += x.piv > 0 ? x.piv : 0;
  }
  void report() { printf("  %-7s runs=%ld mismatches=%ld overflow=%ld pivots=%ld\n", name, runs, bad, ovf, piv); }
};

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: %s trials maxdim seed\n", argv[0]); return 2; }
  int trials = atoi(argv[1]), maxdim = atoi(argv[2]); std::mt19937 gen(atoi(argv[3]));
  namespace N = LRSNash::arith; namespace R = LRSNashRef::arith;
  Pair<N::I64, R::I64> p64{"int64"};
  Pair<N::I128, R::I128> p128{"int128"};
  Pair<N::MP<>, R::MP<>> pmp{"lrsmp"};
  Pair<N::GMP, R::GMP> pgmp{"gmp"};
  for (int t = 0; t < trials; ++t) {
    long m = std::uniform_int_distribution<int>(1, maxdim)(gen), n = std::uniform_int_distribution<int>(1, maxdim)(gen);
    int kind = t % 5; int vmax = kind == 0 ? 2 : kind == 1 ? 9 : kind == 2 ? 100 : kind == 3 ? 1000000 : 1000;
    std::uniform_int_distribution<int> u(kind == 4 ? -vmax : 0, vmax);
    std::vector<int> A(m * n), B(m * n);
    for (auto &x : A) x = u(gen);
    for (auto &x : B) x = u(gen);
    if (t % 3 == 1 && m > 1) for (long j = 0; j < n; j++) { A[n + j] = A[j]; B[n + j] = B[j]; }
    if (t % 7 == 3 && n > 1) for (long i = 0; i < m; i++) { A[i * n + 1] = A[i * n]; B[i * n + 1] = B[i * n]; }
    if (t % 4 == 2) for (long k = 0; k < m * n; k++) B[k] = vmax - A[k];
    if (t % 23 == 5) for (long k = 0; k < m * n; k++) A[k] = B[k] = 0;  // every profile an equilibrium
    for (long k0 = 1; k0 <= m + n; k0++) {
      p64.check(m, n, A, B, k0); p128.check(m, n, A, B, k0);
      pmp.check(m, n, A, B, k0); pgmp.check(m, n, A, B, k0);
    }
  }
  printf("trials=%d maxdim=%d\n", trials, maxdim);
  p64.report(); p128.report(); pmp.report(); pgmp.report();
  return p64.bad + p128.bad + pmp.bad + pgmp.bad ? 1 : 0;
}
