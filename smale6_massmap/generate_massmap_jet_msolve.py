#!/usr/bin/env python3
"""Generate the exact mass-eliminated five-body fixed-mass second-jet scheme.

On the regular mass-solving chart, A=-F0/FA and C=CN/CD.  A fixed-mass
analytic curve forces R=0 and the first and second directional derivatives of
R,A,C to vanish. Projective tangent gauges are substituted before output. An
exact Rabinowitsch equation performs all stated saturations.
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

s,t,w,a,b,c,d=sp.symbols('s t w a b c d')
xs,xt,xw,ys,yt,yw,zz=sp.symbols('xs xt xw ys yt yw zz')
T=t+w; aa=a**3; u=b**3; v=c**3; h=d**3
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
kR=Rs*xs+Rt*xt+Rw*xw
kA=As*xs+At*xt+Aw*xw
kC=Cs*xs+Ct*xt+Cw*xw
qR=DS(kR)*xs+DT(kR)*xt+DW(kR)*xw+Rs*ys+Rt*yt+Rw*yw
qA=DS(kA)*xs+DT(kA)*xt+DW(kA)*xw+As*ys+At*yt+Aw*yw
qC=DS(kC)*xs+DT(kC)*xt+DW(kC)*xw+Cs*ys+Ct*yt+Cw*yw
if args.chart=='s': sub={xs:1,ys:0}; tangent=[xt,xw,yt,yw]
elif args.chart=='t': sub={xt:1,yt:0}; tangent=[xs,xw,ys,yw]
else: sub={xw:1,yw:0}; tangent=[xs,xt,ys,yt]
kR,kA,kC=[sp.expand(e.subs(sub)) for e in (kR,kA,kC)]
qR,qA,qC=[sp.expand(e.subs(sub)) for e in (qR,qA,qC)]
P=w*(s-1)-t; Q=t+w*(s+1)
Delta=s*w*T*(aa-h)*P*Q*FA*F0*CN
vars_=[s,t,w,a,b,c,d]+tangent+[zz]
eqs=[R,g1,g2,g3,g4,kR,kA,kC,qR,qA,qC,sp.expand(zz*Delta-1)]

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
Path(args.output).write_text(','.join(map(str,vars_))+'\n'+str(args.prime)+'\n'+'\n'.join(q+(',' if i<len(ss)-1 else '') for i,q in enumerate(ss))+'\n')
print(f'output={args.output}')
print(f'variables={len(vars_)} equations={len(eqs)} bytes={Path(args.output).stat().st_size}')
print('base_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [R,g1,g2,g3,g4]]))
print('first_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [kR,kA,kC]]))
print('second_terms='+str([len(sp.Poly(e,*vars_).terms()) for e in [qR,qA,qC]]))
