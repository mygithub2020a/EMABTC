#!/usr/bin/env python3
"""Generate the mass-eliminated, natively saturated five-body second-jet system."""
from __future__ import annotations
import argparse
from pathlib import Path
import sympy as sp

ap=argparse.ArgumentParser()
ap.add_argument('--prime',type=int,required=True)
ap.add_argument('--chart',choices=['s','t','w'],required=True)
ap.add_argument('--output',required=True)
args=ap.parse_args()

s,t,w,a,b,c,d,xs,xt,xw,ys,yt,yw=sp.symbols('s t w a b c d xs xt xw ys yt yw')
vars_=[s,t,w,a,b,c,d,xs,xt,xw,ys,yt,yw]
T=t+w; aa=a**3; u=b**3; v=c**3; h=d**3

# With masses (A,C,A,1,1), write F=FA*A+F0, G=GA*A+G0,
# H=HA*A+HC*C+H0.  On FA*HC != 0 the incidence is the shape
# hypersurface R=0 with rational mass map A=-F0/FA, C=Cnum/Cden.
FA=w*(8*aa-1)
F0=4*T*(2*h-u-v)-4*s*w*(v-u)
GA=s**2*(4*u*s*w-4*u*t-4*u*w+4*v*s*w+4*v*t+4*v*w-s*w)
G0=4*u*s**4*w-4*u*s**3*t-4*u*s**3*w-4*v*s**4*w-4*v*s**3*t-4*v*s**3*w+t+w
HA=-t*(4*u+4*v-1)
HC=4*T*(aa-h)
H0=4*s*t*(v-u)
R=sp.expand(G0*FA-GA*F0)
CN=sp.expand(HA*F0-H0*FA)
CD=sp.expand(FA*HC)

g1=a**2*(1+w**2)-1
g2=b**2*((1-s)**2+t**2)-1
g3=c**2*((1+s)**2+t**2)-1
g4=d**2*(s**2+T**2)-1

def DS(P): return sp.expand(sp.diff(P,s)+sp.diff(P,b)*(1-s)*b**3-sp.diff(P,c)*(1+s)*c**3-sp.diff(P,d)*s*d**3)
def DT(P): return sp.expand(sp.diff(P,t)-sp.diff(P,b)*t*b**3-sp.diff(P,c)*t*c**3-sp.diff(P,d)*T*d**3)
def DW(P): return sp.expand(sp.diff(P,w)-sp.diff(P,a)*w*a**3-sp.diff(P,d)*T*d**3)

Rs,Rt,Rw=DS(R),DT(R),DW(R)
F0s,F0t,F0w=DS(F0),DT(F0),DW(F0)
FAs,FAt,FAw=DS(FA),DT(FA),DW(FA)
CNs,CNt,CNw=DS(CN),DT(CN),DW(CN)
CDs,CDt,CDw=DS(CD),DT(CD),DW(CD)
As=sp.expand(FA*F0s-F0*FAs); At=sp.expand(FA*F0t-F0*FAt); Aw=sp.expand(FA*F0w-F0*FAw)
Cs=sp.expand(CD*CNs-CN*CDs); Ct=sp.expand(CD*CNt-CN*CDt); Cw=sp.expand(CD*CNw-CN*CDw)
kR=sp.expand(Rs*xs+Rt*xt+Rw*xw)
kA=sp.expand(As*xs+At*xt+Aw*xw)
kC=sp.expand(Cs*xs+Ct*xt+Cw*xw)
qR=sp.expand(DS(kR)*xs+DT(kR)*xt+DW(kR)*xw+Rs*ys+Rt*yt+Rw*yw)
qA=sp.expand(DS(kA)*xs+DT(kA)*xt+DW(kA)*xw+As*ys+At*yt+Aw*yw)
qC=sp.expand(DS(kC)*xs+DT(kC)*xt+DW(kC)*xw+Cs*ys+Ct*yt+Cw*yw)

if args.chart=='s': normx,normy=xs-1,ys
elif args.chart=='t': normx,normy=xt-1,yt
else: normx,normy=xw-1,yw

P=w*(s-1)-t; Q=t+w*(s+1)
Delta=sp.expand(s*w*T*(aa-h)*P*Q*FA*F0*CN)
# msolve -S interprets the final polynomial as the saturation element.
eqs=[R,g1,g2,g3,g4,kR,kA,kC,qR,qA,qC,normx,normy,Delta]

def pstr(expr):
    poly=sp.Poly(sp.expand(expr),*vars_,domain=sp.ZZ)
    out=[]
    for mon,coef in poly.terms():
        c=int(coef)
        fac=[]
        for var,e in zip(vars_,mon):
            if e==1: fac.append(str(var))
            elif e>1: fac.append(f'{var}^{e}')
        m='*'.join(fac)
        term=str(abs(c)) if not m else (m if abs(c)==1 else f'{abs(c)}*{m}')
        out.append((('-' if c<0 else '+') if out else ('-' if c<0 else ''))+term)
    return ''.join(out) if out else '0'

ss=[pstr(e) for e in eqs]
lines=[','.join(map(str,vars_)),str(args.prime)]+[q+(',' if i<len(ss)-1 else '') for i,q in enumerate(ss)]
Path(args.output).write_text('\n'.join(lines)+'\n')
print({'file':args.output,'chart':args.chart,'prime':args.prime,'variables':len(vars_),'equations':len(eqs),'bytes':Path(args.output).stat().st_size})
