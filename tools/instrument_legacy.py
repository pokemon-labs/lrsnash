#!/usr/bin/env python3
"""Extract the legacy library into OUTDIR/legacy/ and instrument it.

The legacy sources (src/lib.cc, lrslib, lrsmp/lrslong/lrsgmp) were deleted
from the tree when lib.h was switched to the header-only solvers. They are
kept at the git tag `legacy`; this script takes them from there (or from
REF). Tests that compare against legacy include "legacy/lib.h" with
-IOUTDIR and link the objects built below.

lib.cc is patched so that, in the GMP build, every lrs_nashoutput call
appends (player, values...) to the global
std::vector<std::vector<double>> g_seq. test_api.cc reads it. Other builds
do not record.

usage: python3 tools/instrument_legacy.py OUTDIR [REF]
then, in OUTDIR/legacy:
  for f in lib lrslib lrsgmp; do g++ -O2 -c -I. $f.cc -DGMP -w -fpermissive -o gmp_$f.o; done
  for f in lib lrslib lrsmp; do g++ -O2 -c -I. $f.cc -DMP -DLRS_QUIET -w -fpermissive -o mp_$f.o; done
test_api.cc links gmp_*.o; test_lp.cc's legacy mode links mp_*.o."""
import os
import subprocess
import sys

out = os.path.join(sys.argv[1], 'legacy')
ref = sys.argv[2] if len(sys.argv) > 2 else 'legacy'
repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def git(*args):
    return subprocess.run(['git', '-C', repo] + list(args), check=True,
                          capture_output=True).stdout


os.makedirs(out, exist_ok=True)
names = git('ls-tree', '--name-only', ref + ':src').decode().split()
for f in names:
    if f.endswith('.cc') or f.endswith('.h'):
        with open(os.path.join(out, f), 'wb') as fh:
            fh.write(git('show', '%s:src/%s' % (ref, f)))
p = os.path.join(out, 'lib.cc')
t = open(p).read()
t = '#include <vector>\nstd::vector<std::vector<double>> g_seq;\n' + t
old = '  if (origin)\n    return FALSE;\n'
assert old in t
t = t.replace(old, old + '''#ifdef GMP
  { std::vector<double> rec; rec.push_back((double)player);
    for (i = 1; i < Q->n; i++) rec.push_back(mpz_get_d(output[i]) / mpz_get_d(output[0]));
    g_seq.push_back(rec); }
#endif
''', 1)
open(p, 'w').write(t)
print('legacy sources from %s, instrumented, in %s' % (ref, out))
