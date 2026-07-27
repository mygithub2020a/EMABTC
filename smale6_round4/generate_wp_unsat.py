#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sympy as sp

ap=argparse.ArgumentParser()
ap.add_argument('--prime',type=int,required=True)
ap.add_argument('--chart',choices=['W','s','p'],required=True)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args()
W,s,p,aa,bb,cc,dd,A,C,xW,xs,xp,yW,ys,yp=sp.symbols('W s p a b c d A C xW xs xp yW ys yp')
vars0=[W,s,p,aa,bb,cc,dd,A,C,xW,xs,xp,yW,ys,yp]
a3=aa**3;u=bb**3;v=cc**3;h=dd**3
F=4*(p+1)*(2*h-u-v)+(8*A*a3-A-4*s*(v-u))
G=4*A*s**3*u+4*A*s**3*v-A*s**3-4*A*s**2*p*u+4*A*s**2*p*v-4*A*s**2*u+4*A*s**2*v+4*s**4*u-4*s**4*v-4*s**3*p*u-4*s**3*p*v-4*s**3*u-4*s**3*v+p+1
H=4*C*(p+1)*(a3-h)-p*(A*(4*u+4*v-1)-4*s*(v-u))
g=[aa**2*(1+W)-1,bb**2*((1-s)**2+p**2*W)-1,cc**2*((1+s)**2+p**2*W)-1,dd**2*(s**2+(p+1)**2*W)-1]
def DW(f):return sp.expand(2*sp.diff(f,W)-aa**3*sp.diff(f,aa)-p**2*bb**3*sp.diff(f,bb)-p**2*cc**3*sp.diff(f,cc)-(p+1)**2*dd**3*sp.diff(f,dd))
def DS(f):return sp.expand(sp.diff(f,s)+(1-s)*bb**3*sp.diff(f,bb)-(1+s)*cc**3*sp.diff(f,cc)-s*dd**3*sp.diff(f,dd))
def DP(f):return sp.expand(sp.diff(f,p)-p*W*bb**3*sp.diff(f,bb)-p*W*cc**3*sp.diff(f,cc)-(p+1)*W*dd**3*sp.diff(f,dd))
D=[DW,DS,DP]; xi=[xW,xs,xp]; eta=[yW,ys,yp]
J=[[op(f) for op in D] for f in (F,G,H)]
k=[sp.expand(sum(J[i][j]*xi[j] for j in range(3))) for i in range(3)]
q=[sp.expand(sum(D[j](k[i])*xi[j] for j in range(3))+sum(J[i][j]*eta[j] for j in range(3))) for i in range(3)]
idx={'W':0,'s':1,'p':2}[a.chart]
eq=[F,G,H,*g,*k,*q,xi[idx]-1,eta[idx]]
# Linear variables first was the fastest tested DRL order.
order=[yW,ys,yp,A,C,xW,xs,xp,W,s,p,aa,bb,cc,dd]
def text(f):
 poly=sp.Poly(sp.expand(f),*order,domain=sp.ZZ); out=[]
 for mon,coef in poly.terms():
  c=int(coef); c=c%a.prime if a.prime else c
  if not c: continue
  fac=[]
  for z,e in zip(order,mon):
   if e==1:fac.append(str(z))
   elif e>1:fac.append(f'{z}^{e}')
  m='*'.join(fac)
  term=str(c) if not m else (m if c==1 else f'{c}*{m}')
  out.append(term)
 return '+'.join(out) if out else '0'
ss=[text(f) for f in eq]
a.output.write_text('\n'.join([','.join(map(str,order)),str(a.prime)]+[s+(',' if i+1<len(ss) else '') for i,s in enumerate(ss)])+'\n')
print({'vars':len(order),'eqs':len(eq),'terms':[len(sp.Poly(e,*order).terms()) for e in eq],'bytes':a.output.stat().st_size})
