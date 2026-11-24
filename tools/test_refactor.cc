// Refactor gate for the lrsnash enumeration (src/lrsnash/nash.hpp, lrs.hpp):
// the current headers must produce exactly the same visitor call sequence
// (exact integer outputs, in order), the same early-stop behaviour and the
// same pivot count as a frozen reference copy, for all four arithmetic types.
// Solver objects are reused across games of changing shape, so dictionary
// pool reuse is exercised too. Each game is then solved again on the warm
// new solver: the result must repeat and no dictionary may be allocated.
//
// The reference is the pre-refactor commit with namespace LRSNash renamed:
//   REF=/path/to/ref; mkdir -p $REF/lrsnash_ref
//   for f in arith lrs nash; do git show fff2737:src/lrsnash/$f.hpp |
//     sed 's/\bLRSNash\b/LRSNashRef/g' > $REF/lrsnash_ref/$f.hpp; done
// (fff2737 has no Nash::pivot_count; the reference side reads count[3].)
//   g++ -O2 -std=c++17 -Isrc -I$REF tools/test_refactor.cc -lgmp -o test_refactor
//   ./test_refactor trials maxdim seed
#include "lrsnash/nash.hpp"
#include "lrsnash_ref/nash.hpp"
#include <gmp.h>
#include <algorithm>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static void to_mpz(mpz_t z, const int64_t &a) { mpz_set_si(z, a); }
static void to_mpz(mpz_t z, const __int128 &a) {
  unsigned __int128 u = a < 0 ? -(unsigned __int128)a : a;
  mpz_set_ui(z, (uint64_t)(u >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (uint64_t)u);
  if (a < 0) mpz_neg(z, z);
}
template <class G> static auto to_mpz(mpz_t z, const G &a) -> decltype(mpz_set(z, a.z)) { mpz_set(z, a.z); }
template <class M> static auto to_mpz(mpz_t z, const M &a) -> decltype((void)a.v[0]) {  // lrsmp: base 1e9, v[0] = +-length
  long la = a[0] > 0 ? a[0] : -a[0]; mpz_set_ui(z, 0);
  for (long i = la - 1; i >= 1; i--) { mpz_mul_ui(z, z, 1000000000UL); mpz_add_ui(z, z, (unsigned long)a[i]); }
  if (a[0] < 0) mpz_neg(z, z);
}

// Records every visitor call as "<player>:<n>:v0,v1,..." and stops after
// `limit` calls (limit < 0: never).
struct Rec {
  std::vector<std::string> ev;
  long limit = -1;
  template <class T> bool put(char who, const T *o, long n) {
    std::string s(1, who); s += ':' + std::to_string(n) + ':';
    mpz_t z; mpz_init(z);
    for (long k = 0; k < n; k++) { to_mpz(z, o[k]); char *c = mpz_get_str(nullptr, 10, z); s += c; s += ','; free(c); }
    mpz_clear(z);
    ev.push_back(std::move(s));
    return limit < 0 || (long)ev.size() < limit;
  }
  template <class T> bool on_eq2(const T *o, long n) { return put('2', o, n); }
  template <class T> bool on_eq1(const T *o, long n) { return put('1', o, n); }
};

struct Result { std::vector<std::string> ev; long piv = -1; std::string err; };
template <class En> auto pivots(En &e, int) -> decltype(e.pivot_count()) { return e.pivot_count(); }
template <class En> long pivots(En &e, long) { return e.Q1.count[3] + e.Q2.count[3]; }

template <class En> Result run(En &en, long m, long n, const std::vector<int> &A, const std::vector<int> &B, long limit) {
  Result r; Rec v; v.limit = limit;
  try {
    en.solve(m, n, [&](long i, long j) { return (long)A[i * n + j]; }, [&](long i, long j) { return (long)B[i * n + j]; }, v);
    r.piv = pivots(en, 0);
  } catch (std::overflow_error &) { r.err = "overflow"; }
  catch (std::exception &e) { r.err = e.what(); }
  r.ev = std::move(v.ev);
  return r;
}

template <class New, class Ref> struct Pair {
  const char *name;
  LRSNash::detail::Nash<New> a;
  LRSNashRef::detail::Nash<Ref> b;
  long games = 0, bad = 0, ovf = 0, events = 0, piv = 0, bad_repeat = 0, bad_alloc = 0;
  explicit Pair(const char *nm) : name(nm) {}
  void check(long m, long n, const std::vector<int> &A, const std::vector<int> &B, long limit) {
    Result x = run(a, m, n, A, B, limit), y = run(b, m, n, A, B, limit);
    ++games;
    if (x.ev != y.ev || x.piv != y.piv || x.err != y.err) {
      if (++bad <= 3) printf("  MISMATCH %s %ldx%ld limit=%ld: events %zu vs %zu, pivots %ld vs %ld, err '%s' vs '%s'\n",
                             name, m, n, limit, x.ev.size(), y.ev.size(), x.piv, y.piv, x.err.c_str(), y.err.c_str());
      if (bad == 1)  // the first differing event, new vs reference
        for (size_t k = 0; k < std::max(x.ev.size(), y.ev.size()); k++) {
          const char *ex = k < x.ev.size() ? x.ev[k].c_str() : "-", *ey = k < y.ev.size() ? y.ev[k].c_str() : "-";
          if (std::string(ex) != ey) { printf("    event %zu: new %s\n             ref %s\n", k, ex, ey); break; }
        }
    }
    // repeat on the warm solver: same result, and the pool is reused (the
    // pre-refactor code reallocated every dictionary on every solve)
    const auto alloc = a.dict_allocations();
    Result x2 = run(a, m, n, A, B, limit);
    if (x2.ev != x.ev || x2.piv != x.piv || x2.err != x.err) ++bad_repeat;
    if (a.dict_allocations() != alloc) ++bad_alloc;
    if (!x.err.empty()) ++ovf;
    events += x.ev.size(); piv += x.piv > 0 ? x.piv : 0;
  }
  void report() {
    printf("  %-7s games=%ld mismatches=%ld repeat_mismatches=%ld repeat_allocations=%ld overflow=%ld events=%ld pivots=%ld\n",
           name, games, bad, bad_repeat, bad_alloc, ovf, events, piv);
  }
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
    // degenerate structure: duplicated rows/columns, constant-sum
    if (t % 3 == 1 && m > 1) for (long j = 0; j < n; j++) { A[n + j] = A[j]; B[n + j] = B[j]; }
    if (t % 7 == 3 && n > 1) for (long i = 0; i < m; i++) { A[i * n + 1] = A[i * n]; B[i * n + 1] = B[i * n]; }
    if (t % 4 == 2) for (long k = 0; k < m * n; k++) B[k] = vmax - A[k];
    // full enumeration, plus an early stop after 1..3 visitor calls
    long limit = t % 2 ? -1 : std::uniform_int_distribution<int>(1, 3)(gen);
    p64.check(m, n, A, B, limit); p128.check(m, n, A, B, limit);
    pmp.check(m, n, A, B, limit); pgmp.check(m, n, A, B, limit);
  }
  printf("trials=%d maxdim=%d\n", trials, maxdim);
  p64.report(); p128.report(); pmp.report(); pgmp.report();
  long fail = 0;
  for (long v : {p64.bad, p128.bad, pmp.bad, pgmp.bad, p64.bad_repeat, p128.bad_repeat, pmp.bad_repeat, pgmp.bad_repeat,
                 p64.bad_alloc, p128.bad_alloc, pmp.bad_alloc, pgmp.bad_alloc})
    fail += v;
  return fail ? 1 : 0;
}
