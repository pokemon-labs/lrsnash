// Differential test: stripped LRSNash vs upstream lrsnash binary.
#include "lib.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <sstream>
#include <exception>

static double parse_rat(const std::string &t) {
  auto p = t.find('/');
  if (p == std::string::npos) return std::stod(t);
  return std::stod(t.substr(0, p)) / std::stod(t.substr(p + 1));
}

struct Eq { std::vector<double> x, y; double v1 = 0; bool ok = false; int neq = -1; };

static Eq run_orig(const char *bin, int m, int n, const std::vector<long> &A, const std::vector<long> &B, const char *path) {
  FILE *f = fopen(path, "w");
  fprintf(f, "%d %d\n\n", m, n);
  for (int i = 0; i < m; ++i) { for (int j = 0; j < n; ++j) fprintf(f, "%ld ", A[i*n+j]); fprintf(f, "\n"); }
  fprintf(f, "\n");
  for (int i = 0; i < m; ++i) { for (int j = 0; j < n; ++j) fprintf(f, "%ld ", B[i*n+j]); fprintf(f, "\n"); }
  fclose(f);
  std::string cmd = std::string(bin) + " " + path + " 2>&1";
  FILE *p = popen(cmd.c_str(), "r");
  char buf[1 << 16];
  Eq e; std::vector<double> last2; double last2v = 0;
  while (fgets(buf, sizeof buf, p)) {
    std::istringstream ss(buf); std::string tok; std::vector<std::string> toks;
    while (ss >> tok) toks.push_back(tok);
    if (toks.empty()) continue;
    if (toks[0] == "*Number") { e.neq = std::atoi(toks[4].c_str()); continue; }
    if (toks[0] == "2" && (int)toks.size() == n + 2) {
      last2.clear(); for (int j = 1; j <= n; ++j) last2.push_back(parse_rat(toks[j])); last2v = parse_rat(toks[n+1]);
    } else if (toks[0] == "1" && (int)toks.size() == m + 2) {
      e.x.clear(); for (int i = 1; i <= m; ++i) e.x.push_back(parse_rat(toks[i]));
      e.y = last2; e.v1 = last2v; e.ok = true;
    }
  }
  pclose(p);
  return e;
}

int main(int argc, char **argv) {
  if (argc < 7) { fprintf(stderr, "usage: %s orig_bin mode(fast|full) trials maxdim maxval seed\n", argv[0]); return 1; }
  const char *bin = argv[1]; bool fast = !strcmp(argv[2], "fast");
  int trials = atoi(argv[3]), maxdim = atoi(argv[4]), maxval = atoi(argv[5]);
  std::mt19937 gen(atoi(argv[6]));
  std::string path = std::string("g") + argv[6] + (fast ? "f" : "u");
  int fails = 0, exc = 0, noeq = 0, overflow_writes = 0;
  for (int t = 0; t < trials; ++t) {
    int m = std::uniform_int_distribution<int>(1, maxdim)(gen), n = std::uniform_int_distribution<int>(1, maxdim)(gen);
    // vary range: sometimes tiny range (degenerate), sometimes duplicate rows/cols
    int vmax = (t % 3 == 0) ? 2 : maxval;
    std::uniform_int_distribution<int> dis(0, vmax);
    std::vector<int> a(m*n), b(m*n);
    for (auto &x : a) x = dis(gen);
    for (auto &x : b) x = dis(gen);
    if (t % 5 == 1 && m > 1) for (int j = 0; j < n; ++j) { a[n + j] = a[j]; b[n + j] = b[j]; }  // duplicate row
    int den = vmax;
    std::vector<long> A(m*n), B(m*n);
    for (int k = 0; k < m*n; ++k) { A[k] = a[k]; B[k] = fast ? den - a[k] : b[k]; }
    // canary-padded outputs: exactly m / n floats then sentinels
    std::vector<float> rs(m + 4, -777.f), cs(n + 4, -777.f);
    LRSNash::FloatOneSumOutput out{rs.data(), cs.data(), -999.f};
    try {
      if (fast) { LRSNash::FastInput in{m, n, a.data(), den}; LRSNash::solve_fast(&in, &out); }
      else { LRSNash::Input in{m, n, den, a.data(), b.data()}; LRSNash::solve_full(&in, &out); }
    } catch (std::exception &) { ++exc; continue; }
    if (rs[m] != -777.f || cs[n] != -777.f) ++overflow_writes;
    Eq e = run_orig(bin, m, n, A, B, path.c_str());
    if (!e.ok) { ++noeq; fprintf(stderr, "no eq parsed t=%d m=%d n=%d\n", t, m, n); continue; }
    double err = 0;
    for (int i = 0; i < m; ++i) err = std::max(err, std::fabs(e.x[i] - rs[i]));
    for (int j = 0; j < n; ++j) err = std::max(err, std::fabs(e.y[j] - cs[j]));
    double v = e.v1 / den;
    double verr = std::fabs(v - out.value);
    if (err > 1e-5 || verr > 1e-4 * std::max(1.0, std::fabs(v))) {
      if (++fails <= 5) fprintf(stderr, "MISMATCH t=%d m=%d n=%d err=%g v_orig=%g v_strip=%g\n", t, m, n, err, v, out.value);
    }
  }
  printf("%s %s trials=%d fails=%d exceptions=%d noeq=%d writes_past_end=%d\n", bin, argv[2], trials, fails, exc, noeq, overflow_writes);
  return fails != 0;
}
