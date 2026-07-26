#!/usr/bin/env python3
from pathlib import Path
import argparse
import sympy as sp

ap=argparse.ArgumentParser()
ap.add_argument('--prime',type=int,required=True)
ap.add_argument('--chart',choices=['s','t','w'],required=True)
ap.add_argument('--output',required=True)
args=ap.parse_args()

s,t,w,a,b,c,d,rs,rt,rw,ys,yt,yw,z,A,C=sp.symbols('s t w a b c d rs rt rw ys yt yw z A C')
vars_=[s,t,w,a,b,c,d,rs,rt,rw,ys,yt,yw,z,A,C]
T=t+w; aa=a**3; u=b**3; v=c**3; h=d**3
F4=4*T*(2*h-u-v)+w*(8*A*aa-A-4*s*(v-u))
G4=(4*A*s**3*u*w+4*A*s**3*v*w-A*s**3*w-4*A*s**2*t*u+4*A*s**2*t*v-4*A*s**2*u*w+4*A*s**2*v*w+4*s**4*u*w-4*s**4*v*w-4*s**3*t*u-4*s**3*t*v-4*s**3*u*w-4*s**3*v*w+t+w)
H4=4*C*T*(aa-h)-t*(A*(4*u+4*v-1)-4*s*(v-u))
g1=a**2*(1+w**2)-1
g2=b**2*((1-s)**2+t**2)-1
g3=c**2*((1+s)**2+t**2)-1
g4=d**2*(s**2+T**2)-1

def DS(f): return sp.diff(f,s)+(1-s)*b**3*sp.diff(f,b)-(1+s)*c**3*sp.diff(f,c)-s*d**3*sp.diff(f,d)
def DT(f): return sp.diff(f,t)-t*b**3*sp.diff(f,b)-t*c**3*sp.diff(f,c)-T*d**3*sp.diff(f,d)
def DW(f): return sp.diff(f,w)-w*a**3*sp.diff(f,a)-T*d**3*sp.diff(f,d)
Fs,Ft,Fw=DS(F4),DT(F4),DW(F4)
Gs,Gt,Gw=DS(G4),DT(G4),DW(G4)
Hs,Ht,Hw=DS(H4),DT(H4),DW(H4)
Fr=Fs*rs+Ft*rt+Fw*rw; Gr=Gs*rs+Gt*rt+Gw*rw; Hr=Hs*rs+Ht*rt+Hw*rw
Frr=rs*DS(Fr)+rt*DT(Fr)+rw*DW(Fr)
Grr=rs*DS(Gr)+rt*DT(Gr)+rw*DW(Gr)
Hrr=rs*DS(Hr)+rt*DT(Hr)+rw*DW(Hr)
Fy=Fs*ys+Ft*yt+Fw*yw+Frr
Gy=Gs*ys+Gt*yt+Gw*yw+Grr
Hy=Hs*ys+Ht*yt+Hw*yw+Hrr
P=w*(s-1)-t; Q=t+w*(s+1); Delta=A*C*s*w*T*(aa-h)*P*Q
charts={'s':(rs,ys),'t':(rt,yt),'w':(rw,yw)}
rv,yv=charts[args.chart]
eqs=[F4,G4,H4,g1,g2,g3,g4,Fr,Gr,Hr,Fy,Gy,Hy,rv-1,yv,z*Delta-1]

def pstr(expr,p):
    poly=sp.Poly(sp.expand(expr),*vars_,domain=sp.ZZ)
    out=[]
    for mon,coef in poly.terms():
        ci=int(coef)
        if p: ci%=p
        if ci==0: continue
        fac=[]
        for var,e in zip(vars_,mon):
            if e==1: fac.append(str(var))
            elif e>1: fac.append(f'{var}^{e}')
        m='*'.join(fac)
        if not m: term=str(ci)
        elif ci==1: term=m
        elif not p and ci==-1: term='-'+m
        else: term=f'{ci}*{m}'
        out.append(term)
    return '+'.join(out).replace('+-','-') if out else '0'

ss=[pstr(e,args.prime) for e in eqs]
lines=[','.join(map(str,vars_)),str(args.prime)]+[q+(',' if i<len(ss)-1 else '') for i,q in enumerate(ss)]
Path(args.output).write_text('\n'.join(lines)+'\n')
print(args.output,len(eqs),'equations',Path(args.output).stat().st_size,'bytes')
