// Benchmark for the header-only solvers (src/lrsnash/solve.hpp).
// Same game generator as benchmark.cc / the legacy bench: mt19937(seed),
// payoffs uniform in [0, den+1]. For a general-sum game ("gs") the column
// player's matrix is an independent draw; a constant-sum game ("cs") gives
// the column player den - A.
//
//   g++ -O2 -std=c++17 -Isrc tools/bench.cc -lgmp -o bench
//   ./bench rows cols den trials seed MODE auto|i64|i128|mp|gmp
// MODE is <game>-<answer>: game cs|gs, answer any|first|last|all (all =
// extreme_equilibria, else solve with that Mode). The old names still work:
// fast = cs-any, full = gs-last, first = gs-first, any = gs-any, all = gs-all.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>

#include "lrsnash/solve.hpp"

using namespace LRSNash;

// MODE -> (constant-sum?, answer); false if unknown
static bool parse_mode(std::string mode, bool &cs, std::string &ans) {
  if (mode == "fast") mode = "cs-any";
  else if (mode == "full") mode = "gs-last";
  else if (mode == "first" || mode == "any" || mode == "all") mode = "gs-" + mode;
  if (mode.size() < 4 || mode[2] != '-') return false;
  const std::string game = mode.substr(0, 2);
  ans = mode.substr(3);
  cs = game == "cs";
  return (cs || game == "gs") && (ans == "any" || ans == "first" || ans == "last" || ans == "all");
}

template <class S>
static void run(S &s, int rows, int cols, int den, int trials, unsigned seed,
                const std::string &mode) {
  bool cs;
  std::string ans;
  if (!parse_mode(mode, cs, ans)) {
    fprintf(stderr, "unknown MODE %s\n", mode.c_str());
    exit(1);
  }
  const Mode md = ans == "first" ? Mode::First : ans == "last" ? Mode::Last : Mode::Any;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dis(0, den + 1);
  std::vector<int> a(rows * cols), b(rows * cols);
  std::vector<float> rs(rows), cols_out(cols);
  double chk = 0;
  long ovf = 0, neq = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int t = 0; t < trials; ++t) {
    for (auto &x : a) x = dis(gen);
    if (!cs)
      for (auto &x : b) x = dis(gen);
    FloatEquilibrium out{rs.data(), cols_out.data(), 0.f};
    try {
      if (cs) {
        ConstantSumGame in{rows, cols, a.data(), den};
        if (ans == "all") neq += (long)s.extreme_equilibria(&in).size();
        else s.solve(&in, &out, md);
      } else {
        GeneralSumGame in{rows, cols, den, a.data(), b.data()};
        if (ans == "all") neq += (long)s.extreme_equilibria(&in).size();
        else s.solve(&in, &out, md);
      }
    } catch (overflow_error &) {
      ++ovf;
      continue;
    }
    chk += out.value;
  }
  double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  printf("%dx%d den=%d %-5s trials=%d: %10.2f us/game  overflow=%ld%s chk=%.6f\n",
         rows, cols, den, mode.c_str(), trials, 1e6 * sec / trials, ovf,
         ans == "all" ? (" eq/game=" + std::to_string((double)neq / trials)).c_str() : "", chk);
}

int main(int argc, char **argv) {
  if (argc < 8) {
    fprintf(stderr, "usage: %s rows cols den trials seed {cs,gs}-{any,first,last,all} auto|i64|i128|mp|gmp\n", argv[0]);
    return 1;
  }
  int rows = atoi(argv[1]), cols = atoi(argv[2]), den = atoi(argv[3]), trials = atoi(argv[4]);
  unsigned seed = atoi(argv[5]);
  std::string mode = argv[6], ar = argv[7];
  bool cs_unused;
  std::string ans_unused;
  if (!parse_mode(mode, cs_unused, ans_unused)) {
    fprintf(stderr, "unknown MODE %s\n", mode.c_str());
    return 1;
  }
  printf("[%-4s] ", ar.c_str());
  if (ar == "i64") { Solver64 s; run(s, rows, cols, den, trials, seed, mode); }
  else if (ar == "i128") { Solver128 s; run(s, rows, cols, den, trials, seed, mode); }
  else if (ar == "mp") { SolverMP s; run(s, rows, cols, den, trials, seed, mode); }
#ifdef LRSNASH_HAVE_GMP
  else if (ar == "gmp") { SolverGMP s; run(s, rows, cols, den, trials, seed, mode); }
#endif
  else { Solver s; run(s, rows, cols, den, trials, seed, mode); }
}
