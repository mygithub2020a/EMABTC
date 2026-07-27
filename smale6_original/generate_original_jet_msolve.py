#!/usr/bin/env python3
"""Generate a fixed-mass second-jet system for the regular symmetric 5-body gate.

The system uses the exact inverse-distance cover and a Rabinowitsch equation
zz*Delta-1, so msolve is invoked without its native -S mode.  A genuine
nonconstant fixed-mass analytic curve supplies a one-dimensional component in
at least one of the s, t, w projective tangent charts.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sympy as sp

ap=argparse.ArgumentParser()
ap.add_argument('--prime',type=int,required=True)
ap.add_argument('--chart',choices=['s','t','w'],required=True)
ap.add_argument('--output',required=True)
args=ap.parse_args()

s,t,w,a,b,c,d,A,C,xs,xt,xw,ys,yt,yw,zz=sp.symbols(
    's t w a b c d A C xs xt xw ys yt yw zz'
)
vars_=[s,t,w,a,b,c,d,A,C,xs,xt,xw,ys,yt,yw,zz]
T=t+w; aa=a**3; u=b**3; v=c**3; h=d**3
F=4*T*(2*h-u-v)+w*(8*A*aa-A-4*s*(v-u))
G=(4*A*s**3*u*w+4*A*s**3*v*w-A*s**3*w
   -4*A*s**2*t*u+4*A*s**2*t*v-4*A*s**2*u*w+4*A*s**2*v*w
   +4*s**4*u*w-4*s**4*v*w-4*s**3*t*u-4*s**3*t*v
   -4*s**3*u*w-4*s**3*v*w+t+w)
H=4*C*T*(aa-h)-t*(A*(4*u+4*v-1)-4*s*(v-u))
g1=a**2*(1+w**2)-1
g2=b**2*((1-s)**2+t**2)-1
g3=c**2*((1+s)**2+t**2)-1
g4=d**2*(s**2+T**2)-1

def DS(f): return sp.diff(f,s)+(1-s)*b**3*sp.diff(f,b)-(1+s)*c**3*sp.diff(f,c)-s*d**3*sp.diff(f,d)
def DT(f): return sp.diff(f,t)-t*b**3*sp.diff(f,b)-t*c**3*sp.diff(f,c)-T*d**3*sp.diff(f,d)
def DW(f): return sp.diff(f,w)-w*a**3*sp.diff(f,a)-T*d**3*sp.diff(f,d)

Fs,Ft,Fw=DS(F),DT(F),DW(F)
Gs,Gt,Gw=DS(G),DT(G),DW(G)
Hs,Ht,Hw=DS(H),DT(H),DW(H)
kF=sp.expand(Fs*xs+Ft*xt+Fw*xw)
kG=sp.expand(Gs*xs+Gt*xt+Gw*xw)
kH=sp.expand(Hs*xs+Ht*xt+Hw*xw)
qF=sp.expand(DS(kF)*xs+DT(kF)*xt+DW(kF)*xw+Fs*ys+Ft*yt+Fw*yw)
qG=sp.expand(DS(kG)*xs+DT(kG)*xt+DW(kG)*xw+Gs*ys+Gt*yt+Gw*yw)
qH=sp.expand(DS(kH)*xs+DT(kH)*xt+DW(kH)*xw+Hs*ys+Ht*yt+Hw*yw)
charts={'s':(xs-1,ys),'t':(xt-1,yt),'w':(xw-1,yw)}
normx,normy=charts[args.chart]
P=w*(s-1)-t; Q=t+w*(s+1)
Delta=A*C*s*w*T*(aa-h)*P*Q*(8*aa-1)
eqs=[F,G,H,g1,g2,g3,g4,kF,kG,kH,qF,qG,qH,normx,normy,
     sp.expand(zz*Delta-1)]

def pstr(expr,pmod):
    poly=sp.Poly(sp.expand(expr),*vars_,domain=sp.ZZ)
    out=[]
    for mon,coef in poly.terms():
        ci=int(coef)
        if pmod: ci%=pmod
        if ci==0: continue
        fac=[]
        for var,e in zip(vars_,mon):
            if e==1: fac.append(str(var))
            elif e>1: fac.append(f'{var}^{e}')
        m='*'.join(fac)
        if not m: term=str(ci)
        elif ci==1: term=m
        elif not pmod and ci==-1: term='-'+m
        else: term=f'{ci}*{m}'
        out.append(term)
    return '+'.join(out).replace('+-','-') if out else '0'

ss=[pstr(e,args.prime) for e in eqs]
lines=[','.join(map(str,vars_)),str(args.prime)]
lines += [q+(',' if i<len(ss)-1 else '') for i,q in enumerate(ss)]
Path(args.output).write_text('\n'.join(lines)+'\n')
print(f'output={args.output}')
print(f'variables={len(vars_)} equations={len(eqs)} bytes={Path(args.output).stat().st_size}')
print('base_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [F,G,H,g1,g2,g3,g4]]))
print('first_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [kF,kG,kH]]))
print('second_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [qF,qG,qH]]))
