#!/usr/bin/env python3
"""Generate the exact four-distance second-jet system for msolve.

The system is necessary for a fixed-mass analytic curve in the regular
reflection-symmetric five-body family. It uses one projective tangent chart
and the induced parameterization gauge on the acceleration. The final input
polynomial is the native msolve saturation factor used with `-S`.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument("--prime", type=int, required=True)
ap.add_argument("--chart", choices=["a", "b", "c", "d"], required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

a,b,c,d,A,C,xa,xb,xc,xd,ya,yb,yc,yd = sp.symbols(
    "a b c d A C xa xb xc xd ya yb yc yd"
)
shape = [a,b,c,d]
xvec = [xa,xb,xc,xd]
yvec = [ya,yb,yc,yd]
vars_ = [a,b,c,d,A,C,xa,xb,xc,xd,ya,yb,yc,yd]

Y = -(2*b*c-b-c)*(2*b*c-b+c)*(2*b*c+b-c)*(2*b*c+b+c)
NE = (4*a**2*b**2*c**2*d**2 + 2*a**2*b**2*c**2
      - a**2*b**2*d**2 - a**2*c**2*d**2 - 2*b**2*c**2*d**2)
G0 = sp.expand(NE**2-a**2*(1-a**2)*Y*d**4)

s = (b**2-c**2)/(4*b**2*c**2)
W = (1-a**2)/a**2
E = 1/d**2-sp.Rational(1,2)/b**2-sp.Rational(1,2)/c**2+2-1/a**2
p = sp.cancel(E/(2*W))
u=b**3; v=c**3; h=d**3
F = 4*(p+1)*(2*h-u-v)+A*(8*a**3-1)-4*s*(v-u)
G = (4*A*s**3*(u+v)-A*s**3-4*A*s**2*p*u+4*A*s**2*p*v
     -4*A*s**2*u+4*A*s**2*v+4*s**4*u-4*s**4*v
     -4*s**3*(p+1)*(u+v)+p+1)
H = 4*C*(p+1)*(a**3-h)-p*(A*(4*u+4*v-1)-4*s*(v-u))

def numerator(expr: sp.Expr) -> sp.Expr:
    return sp.expand(sp.fraction(sp.cancel(expr))[0])

base = [G0, numerator(F), numerator(G), numerator(H)]
J = sp.Matrix([[sp.diff(f,q) for q in shape] for f in base])
first = [sp.expand(sum(J[i,j]*xvec[j] for j in range(4))) for i in range(4)]
second=[]
for i,f in enumerate(base):
    hxx=0
    for j in range(4):
        for k in range(4):
            hxx += sp.diff(f,shape[j],shape[k])*xvec[j]*xvec[k]
    jy=sum(J[i,j]*yvec[j] for j in range(4))
    second.append(sp.expand(jy+hxx))

idx={"a":0,"b":1,"c":2,"d":3}[args.chart]
norm=[xvec[idx]-1,yvec[idx]]

# Y=0 is collinearity; b^2-c^2=0 is s=0; NE=0 is t*w=0;
# 8*a^3-1 is a separately handled singular mass-solving chart.
Delta=A*C*a*b*c*d*(1-a**2)*Y*(b**2-c**2)*NE*(8*a**3-1)
# With msolve -S, the final polynomial is the saturation factor itself.
eqs=base+first+second+norm+[sp.expand(Delta)]

def pstr(expr: sp.Expr,pmod: int) -> str:
    poly=sp.Poly(sp.expand(expr),*vars_,domain=sp.ZZ)
    out=[]
    for mon,coef in poly.terms():
        ci=int(coef)
        if pmod:
            ci%=pmod
        if ci==0:
            continue
        factors=[]
        for var,ex in zip(vars_,mon):
            if ex==1:
                factors.append(str(var))
            elif ex>1:
                factors.append(f"{var}^{ex}")
        m="*".join(factors)
        if not m:
            term=str(ci)
        elif ci==1:
            term=m
        elif not pmod and ci==-1:
            term="-"+m
        else:
            term=f"{ci}*{m}"
        out.append(term)
    return "+".join(out).replace("+-","-") if out else "0"

ss=[pstr(e,args.prime) for e in eqs]
lines=[",".join(map(str,vars_)),str(args.prime)]
lines += [q+("," if i<len(ss)-1 else "") for i,q in enumerate(ss)]
Path(args.output).write_text("\n".join(lines)+"\n")
print(f"output={args.output}")
print(f"variables={len(vars_)} generators={len(eqs)-1} saturation_factor=1 total_input_polynomials={len(eqs)} bytes={Path(args.output).stat().st_size}")
print("base_terms="+str([len(sp.Poly(e,*vars_).terms()) for e in base]))
print("first_terms="+str([len(sp.Poly(e,*vars_).terms()) for e in first]))
print("second_terms="+str([len(sp.Poly(e,*vars_).terms()) for e in second]))
