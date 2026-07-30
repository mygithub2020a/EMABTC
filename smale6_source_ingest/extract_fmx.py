#!/usr/bin/env python3
from __future__ import annotations
import hashlib,json,math,re,sys
from pathlib import Path

def match(s,start):
    depth=0; i=start; ins=False; esc=False
    while i<len(s):
        c=s[i]
        if ins:
            if esc: esc=False
            elif c=='\\': esc=True
            elif c=='"': ins=False
        else:
            if c=='"': ins=True
            elif c=='[': depth+=1
            elif c==']':
                depth-=1
                if depth==0: return i
        i+=1
    raise ValueError('unmatched')

def main():
    if len(sys.argv)!=4: raise SystemExit('usage: extract_fmx.py INPUT.nb FMX.wl DIAGRAMS.json')
    src=Path(sys.argv[1]); text=src.read_text(errors='replace')
    blocks=[]; pos=0
    while True:
        j=text.find('GraphicsBox[',pos)
        if j<0: break
        a=j+len('GraphicsBox'); b=match(text,a); block=text[j:b+1]
        if block.count('DiskBox[')==6 and '{1.1, -1}' in block: blocks.append(block)
        pos=b+1
    coords={1:(1.,0.),2:(.5,math.sqrt(3)/2),3:(-.5,math.sqrt(3)/2),4:(-1.,0.),5:(-.5,-math.sqrt(3)/2),6:(.5,-math.sqrt(3)/2)}
    def vertex(x,y): return min(coords,key=lambda k:(coords[k][0]-x)**2+(coords[k][1]-y)**2)
    def kind(c):
        c=''.join(c.split())
        if c.startswith('GrayLevel[1]') or c=='White': return 'none'
        if 'RGBColor[1,0.3,0.3]' in c:return 'z'
        if 'RGBColor[0.2,0.7,1]' in c:return 'w'
        if 'RGBColor[0.5,0.3,0.9]' in c:return 'zw'
        return 'unknown'
    epat=re.compile(r'\{(RGBColor\[[^\]]+\])\s*,\s*Thickness\[0\.025\]\s*,\s*LineBox\[\{\{([^,{}]+),\s*([^{}]+?)\},\s*\{([^,{}]+),\s*([^{}]+?)\}\}\]\}',re.S)
    vpat=re.compile(r'\{(GrayLevel\[1\]|RGBColor\[[^\]]+\])\s*,\s*EdgeForm\[GrayLevel\[0\]\]\s*,\s*DiskBox\[\{[^\]]+\},\s*0\.24\].*?InsetBox\["([1-6])",',re.S)
    by={}
    for block in blocks:
        nums=re.findall(r'StyleBox\[InsetBox\["(\d+)", \{1\.1, -1\}\]',block)
        if not nums: continue
        num=int(max(set(nums),key=nums.count)); sep=block.find('{Hue[0.6, 0.2, 0.8]')
        zedges=set(); wedges=set(); zcirc=set(); wcirc=set(); fills={}
        for m in epat.finditer(block[:sep]):
            vals=[float(x.replace('`','')) for x in m.groups()[1:]]
            e=tuple(sorted((vertex(vals[0],vals[1]),vertex(vals[2],vals[3])))); k=kind(m.group(1))
            if k in ('z','zw'): zedges.add(e)
            if k in ('w','zw'): wedges.add(e)
        for m in vpat.finditer(block[sep:]):
            k=kind(m.group(1)); i=int(m.group(2)); fills[i]=k
            if k in ('z','zw'): zcirc.add(i)
            if k in ('w','zw'): wcirc.add(i)
        if len(fills)!=6: raise SystemExit(f'vertex parse failure diagram {num}: {fills}')
        def matrix(edges,circles):
            M=[[0]*6 for _ in range(6)]
            for i in circles:M[i-1][i-1]=1
            for i,j in edges:M[i-1][j-1]=M[j-1][i-1]=1
            return M
        d={'number':num,'z':matrix(zedges,zcirc),'w':matrix(wedges,wcirc),'z_edges':sorted(map(list,zedges)),'w_edges':sorted(map(list,wedges)),'z_circles':sorted(zcirc),'w_circles':sorted(wcirc)}
        if num in by and d!=by[num]: raise SystemExit(f'duplicate mismatch {num}')
        by[num]=d
    if set(by)!=set(range(1,118)): raise SystemExit(f'diagram labels mismatch: {sorted(set(range(1,118))-set(by))}')
    def wl(obj):
        if isinstance(obj,list): return '{'+','.join(wl(x) for x in obj)+'}'
        return str(obj)
    fmx=[[by[i]['z'],by[i]['w']] for i in range(1,118)]
    Path(sys.argv[2]).write_text('FMX='+wl(fmx)+';\n')
    out={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'count':117,'diagrams':[by[i] for i in range(1,118)]}
    Path(sys.argv[3]).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'count':117,'source_sha256':out['source_sha256'],'fmx_bytes':Path(sys.argv[2]).stat().st_size}))
if __name__=='__main__': main()
