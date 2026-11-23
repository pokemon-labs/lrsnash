import re, sys, subprocess
def strip(path):
    src=open(path).read()
    src=re.sub(r'/\*.*?\*/', ' ', src, flags=re.S)
    src=re.sub(r'//[^\n]*', ' ', src)
    return src
def funcs(src):
    src=fmt(src)
    out={}
    for m in re.finditer(r'^(?:[A-Za-z_][\w \t\*]*?[\s\*])?(\w+)\s*\(([^;{)]*(?:\([^)]*\)[^;{)]*)*)\)\s*\{', src, re.M):
        name=m.group(1)
        if name in ('if','while','for','switch'): continue
        i=m.end()-1; depth=0
        for j in range(i,len(src)):
            if src[j]=='{': depth+=1
            elif src[j]=='}':
                depth-=1
                if depth==0: break
        out.setdefault(name,[]).append(src[m.start():j+1])
    return out
def fmt(code):
    return subprocess.run(['clang-format','--style={BasedOnStyle: LLVM, ColumnLimit: 200}'],input=code,capture_output=True,text=True).stdout
a=funcs(strip(sys.argv[1])); b=funcs(strip(sys.argv[2])); outdir=sys.argv[3]
import os, difflib
os.makedirs(outdir,exist_ok=True)
for name in sorted(b):
    if name not in a: print("ONLY-IN-STRIPPED", name); continue
    fa=fmt(a[name][-1]).splitlines(); fb=fmt(b[name][-1]).splitlines()
    d=list(difflib.unified_diff(fa,fb,'orig/'+name,'strip/'+name,lineterm='',n=2))
    open(f'{outdir}/{name}.diff','w').write('\n'.join(d))
    print(f"{name:25s} orig={len(fa):4d} strip={len(fb):4d} difflines={sum(1 for l in d if l[:1] in '+-' and not l.startswith(('+++','---')))}  multi_orig={len(a[name])}")
