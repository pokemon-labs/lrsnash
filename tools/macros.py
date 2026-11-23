import re,sys
def defs(p):
    s=open(p).read().replace('\\\n',' ')
    s=re.sub(r'/\*.*?\*/',' ',s,flags=re.S); s=re.sub(r'//[^\n]*','',s)
    d={}
    # track #if context crudely
    ctx=[]
    for line in s.splitlines():
        t=line.strip()
        m=re.match(r'#\s*(ifdef|ifndef|if)\s+(.*)',t)
        if m: ctx.append(m.group(1)+' '+m.group(2).strip()); continue
        if re.match(r'#\s*else',t): ctx[-1]='!'+ctx[-1]; continue
        if re.match(r'#\s*endif',t): ctx and ctx.pop(); continue
        m=re.match(r'#\s*define\s+(\w+(?:\([^)]*\))?)\s*(.*)',t)
        if m:
            key=m.group(1).split('(')[0]+' ['+' & '.join(ctx)+']'
            d.setdefault(key,[]).append(re.sub(r'\s+','',m.group(1)+m.group(2)))
    return d
a=defs(sys.argv[1]); b=defs(sys.argv[2])
for k in sorted(set(a)|set(b)):
    if a.get(k)!=b.get(k):
        print('==',k); print('  O:',a.get(k)); print('  S:',b.get(k))
