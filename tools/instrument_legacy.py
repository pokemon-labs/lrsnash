#!/usr/bin/env python3
"""Copy the legacy sources (src/*.cc, src/*.h) into OUTDIR and patch lib.cc so
that every lrs_nashoutput call appends (player, values...) to the global
std::vector<std::vector<double>> g_seq (GMP build only). Used by test_api.cc.
usage: python3 tools/instrument_legacy.py OUTDIR"""
import glob, os, shutil, sys
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src')
for f in glob.glob(os.path.join(src, '*.cc')) + glob.glob(os.path.join(src, '*.h')):
    shutil.copy(f, out)
p = os.path.join(out, 'lib.cc')
t = open(p).read()
t = '#include <vector>\nstd::vector<std::vector<double>> g_seq;\n' + t
old = '  if (origin)\n    return FALSE;\n'
assert old in t
t = t.replace(old, old + '''  { std::vector<double> rec; rec.push_back((double)player);
    for (i = 1; i < Q->n; i++) rec.push_back(mpz_get_d(output[i]) / mpz_get_d(output[0]));
    g_seq.push_back(rec); }
''', 1)
open(p, 'w').write(t)
print('instrumented legacy sources in', out)
