#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument("--prime", type=int, required=True)
ap.add_argument("--chart", choices=["W", "s", "p"], required=True)
ap.add_argument("--output", required=True)
args = ap.parse_args()

W, s, p, a, b, c, d = sp.symbols("W s p a b c d")
if args.chart == "W":
    xs, xp, ys, yp = sp.symbols("xs xp ys yp")
    x = (sp.Integer(1), xs, xp)
    y = (sp.Integer(0), ys, yp)
    jetvars = [xs, xp, ys, yp]
elif args.chart == "s":
    xW, xp, yW, yp = sp.symbols("xW xp yW yp")
    x = (xW, sp.Integer(1), xp)
    y = (yW, sp.Integer(0), yp)
    jetvars = [xW, xp, yW, yp]
else:
    xW, xs, yW, ys = sp.symbols("xW xs yW ys")
    x = (xW, xs, sp.Integer(1))
    y = (yW, ys, sp.Integer(0))
    jetvars = [xW, xs, yW, ys]

vars_ = [W, s, p, a, b, c, d] + jetvars
A, C = sp.symbols("A C")
aa, u, v, h = a**3, b**3, c**3, d**3

F = 4*(p+1)*(2*h-u-v) + (8*A*aa-A-4*s*(v-u))
G = (
    4*A*s**3*u + 4*A*s**3*v - A*s**3
    - 4*A*s**2*p*u + 4*A*s**2*p*v
    - 4*A*s**2*u + 4*A*s**2*v
    + 4*s**4*u - 4*s**4*v
    - 4*s**3*p*u - 4*s**3*p*v
    - 4*s**3*u - 4*s**3*v + p + 1
)
H = 4*C*(p+1)*(aa-h) - p*(A*(4*u+4*v-1)-4*s*(v-u))

FA = sp.diff(F, A)
F0 = F.subs(A, 0)
GA = sp.diff(G, A)
G0 = G.subs(A, 0)
HA = sp.diff(H, A)
HC = sp.diff(H, C)
H0 = H.subs({A: 0, C: 0})

# Shape incidence and rational mass map:
# A=-F0/FA, C=Cnum/Cden, and R=0 is compatibility with G=0.
R = sp.expand(G0*FA-GA*F0)
Cnum = sp.expand(HA*F0-H0*FA)
Cden = sp.expand(FA*HC)

g1 = a**2*(1+W)-1
g2 = b**2*((1-s)**2+p**2*W)-1
g3 = c**2*((1+s)**2+p**2*W)-1
g4 = d**2*(s**2+(p+1)**2*W)-1

def DW(f):
    return sp.expand(
        2*sp.diff(f, W)-a**3*sp.diff(f, a)
        -p**2*b**3*sp.diff(f, b)-p**2*c**3*sp.diff(f, c)
        -(p+1)**2*d**3*sp.diff(f, d)
    )

def DS(f):
    return sp.expand(
        sp.diff(f, s)+(1-s)*b**3*sp.diff(f, b)
        -(1+s)*c**3*sp.diff(f, c)-s*d**3*sp.diff(f, d)
    )

def DP(f):
    return sp.expand(
        sp.diff(f, p)-p*W*b**3*sp.diff(f, b)
        -p*W*c**3*sp.diff(f, c)-(p+1)*W*d**3*sp.diff(f, d)
    )

Ds = [DW, DS, DP]
Rgrad = [D(R) for D in Ds]
Agrad = [sp.expand(FA*D(F0)-F0*D(FA)) for D in Ds]
Cgrad = [sp.expand(Cden*D(Cnum)-Cnum*D(Cden)) for D in Ds]
grads = [Rgrad, Agrad, Cgrad]

# First fixed-mass/tangency equations.
k = [sp.expand(sum(gr[j]*x[j] for j in range(3))) for gr in grads]

def Dx(e):
    return sp.expand(x[0]*DW(e)+x[1]*DS(e)+x[2]*DP(e))

# Second-jet equations. On k_A=k_C=0 these are precisely the numerators
# of A''=C''=0; q_R is the second derivative of the incidence equation.
q = [
    sp.expand(Dx(k[i])+sum(grads[i][j]*y[j] for j in range(3)))
    for i in range(3)
]

# Open physical regular chart. The last polynomial is supplied to msolve -S.
factors = [
    W, s, p+1, aa-h, s-1-p, p+s+1,
    FA, F0, Cnum, HC,
]
Delta = sp.expand(sp.prod(factors))
eqs = [R, g1, g2, g3, g4] + k + q + [Delta]

def pstr(expr):
    P = sp.Poly(sp.expand(expr), *vars_, domain=sp.ZZ)
    out = []
    for mon, coef in P.terms():
        ci = int(coef) % args.prime
        if ci == 0:
            continue
        fac = []
        for var, exponent in zip(vars_, mon):
            if exponent == 1:
                fac.append(str(var))
            elif exponent > 1:
                fac.append(f"{var}^{exponent}")
        monomial = "*".join(fac)
        term = str(ci) if not monomial else (monomial if ci == 1 else f"{ci}*{monomial}")
        out.append(("+" if out else "") + term)
    return "".join(out) if out else "0"

strings = [pstr(e) for e in eqs]
Path(args.output).write_text(
    "\n".join(
        [",".join(map(str, vars_)), str(args.prime)]
        + [poly + ("," if i < len(strings)-1 else "") for i, poly in enumerate(strings)]
    ) + "\n"
)
polys = [sp.Poly(sp.expand(e), *vars_) for e in eqs]
print(json.dumps({
    "chart": args.chart,
    "prime": args.prime,
    "variables": len(vars_),
    "equations_before_saturation": len(eqs)-1,
    "native_saturation_polynomial": True,
    "bytes": Path(args.output).stat().st_size,
    "terms": [len(P.terms()) for P in polys],
    "degrees": [P.total_degree() for P in polys],
    "saturation_factors": [str(f) for f in factors],
}, indent=2))
