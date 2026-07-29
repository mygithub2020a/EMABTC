#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--prime',type=int,required=True);ap.add_argument('--row',type=int,choices=[0,1,2],required=True);ap.add_argument('--col',type=int,choices=[0,1,2],required=True);ap.add_argument('--output',required=True);args=ap.parse_args();P=args.prime
names=['W','s','p','a','b','c','d','A','C'];n=9
def const(c):return {} if c%P==0 else {(0,)*n:c%P}
def var(i):m=[0]*n;m[i]=1;return {tuple(m):1}
V=[var(i) for i in range(n)];W,s,p,a,b,c,d,A,C=V
def add(x,y):
 z=dict(x)
 for m,k in y.items():
  z[m]=(z.get(m,0)+k)%P
  if not z[m]:del z[m]
 return z
def neg(x):return {m:(-k)%P for m,k in x.items()}
def sub(x,y):return add(x,neg(y))
def scale(x,k):k%=P;return {} if not k else {m:(v*k)%P for m,v in x.items() if (v*k)%P}
def mul(x,y):
 if not x or not y:return {}
 z={}
 for m1,c1 in x.items():
  for m2,c2 in y.items():
   m=tuple(a0+b0 for a0,b0 in zip(m1,m2));z[m]=(z.get(m,0)+c1*c2)%P
 return {m:c0 for m,c0 in z.items() if c0}
def powp(x,e):
 z=const(1);bb=x
 while e:
  if e&1:z=mul(z,bb)
  e//=2
  if e:bb=mul(bb,bb)
 return z
def diffp(x,i):
 z={}
 for m,c0 in x.items():
  if m[i]:
   mm=list(m);e=mm[i];mm[i]-=1;mm=tuple(mm);z[mm]=(z.get(mm,0)+c0*e)%P
 return {m:c0 for m,c0 in z.items() if c0}
def prod(xs):
 z=const(1)
 for x in xs:z=mul(z,x)
 return z
def det2(a,b,c,d):return sub(mul(a,d),mul(b,c))
def terms(x):return len(x)
def degree(x):return max((sum(m) for m in x),default=-1)
one=const(1);aa=powp(a,3);u=powp(b,3);v=powp(c,3);h=powp(d,3);p1=add(p,one)
F=add(scale(mul(p1,add(scale(h,2),neg(add(u,v)))),4),add(sub(scale(mul(A,aa),8),A),neg(scale(mul(s,sub(v,u)),4))))
G={}
def acc(co,*xs):
 global G
 G=add(G,scale(prod(xs),co))
acc(4,A,powp(s,3),u);acc(4,A,powp(s,3),v);acc(-1,A,powp(s,3));acc(-4,A,powp(s,2),p,u);acc(4,A,powp(s,2),p,v);acc(-4,A,powp(s,2),u);acc(4,A,powp(s,2),v);acc(4,powp(s,4),u);acc(-4,powp(s,4),v);acc(-4,powp(s,3),p,u);acc(-4,powp(s,3),p,v);acc(-4,powp(s,3),u);acc(-4,powp(s,3),v);G=add(G,p1)
H=add(scale(prod([C,p1,sub(aa,h)]),4),neg(mul(p,sub(mul(A,add(scale(add(u,v),4),neg(one))),scale(mul(s,sub(v,u)),4)))))
g1=sub(mul(powp(a,2),add(one,W)),one);g2=sub(mul(powp(b,2),add(powp(sub(one,s),2),mul(powp(p,2),W))),one);g3=sub(mul(powp(c,2),add(powp(add(one,s),2),mul(powp(p,2),W))),one);g4=sub(mul(powp(d,2),add(powp(s,2),mul(powp(p1,2),W))),one)
EWc=[const(2),{},{},neg(aa),neg(mul(powp(p,2),u)),neg(mul(powp(p,2),v)),neg(mul(powp(p1,2),h)),{},{}]
DSc=[{},one,{}, {},mul(sub(one,s),u),neg(mul(add(one,s),v)),neg(mul(s,h)),{},{}]
DPc=[{},{},one,{},neg(mul(mul(p,W),u)),neg(mul(mul(p,W),v)),neg(mul(mul(p1,W),h)),{},{}]
Dcoeff=[EWc,DSc,DPc]
def Dapply(f,coefs):
 z={}
 for k,ck in enumerate(coefs):
  if ck:z=add(z,mul(ck,diffp(f,k)))
 return z
fs=[F,G,H];J=[[Dapply(f,Dcoeff[j]) for j in range(3)] for f in fs]
adj=[[None]*3 for _ in range(3)]
for i in range(3):
 for j in range(3):
  rs=[r for r in range(3) if r!=j];cs=[cc for cc in range(3) if cc!=i]
  q=det2(J[rs[0]][cs[0]],J[rs[0]][cs[1]],J[rs[1]][cs[0]],J[rs[1]][cs[1]])
  if (i+j)%2:q=neg(q)
  adj[i][j]=q
detJ=add(add(mul(J[0][0],adj[0][0]),mul(J[0][1],adj[1][0])),mul(J[0][2],adj[2][0]))
i,j=args.row,args.col;chart=adj[i][j];vcol=[adj[k][j] for k in range(3)]
ddet={}
for k in range(3):ddet=add(ddet,mul(vcol[k],Dapply(detJ,Dcoeff[k])))
FA=sub(scale(aa,8),one);F0=add(scale(mul(p1,add(scale(h,2),neg(add(u,v)))),4),neg(scale(mul(s,sub(v,u)),4)));CN=sub(mul(neg(mul(p,add(scale(add(u,v),4),neg(one)))),F0),mul(scale(mul(mul(s,p),sub(v,u)),4),FA))
Delta=prod([A,C,W,s,p1,a,b,c,d,sub(aa,h),sub(sub(s,one),p),add(add(p,s),one),FA,F0,CN])
baseeq=[F,G,H,g1,g2,g3,g4,detJ,ddet];N=11
def ext(q):return {m+(0,0):co for m,co in q.items()}
def ztimes(q,idx):
 out={}
 for m,co in q.items():mm=list(m)+[0,0];mm[9+idx]=1;out[tuple(mm)]=co
 return out
E=[ext(q) for q in baseeq]+[add(ztimes(Delta,0),{(0,)*N:(-1)%P}),add(ztimes(chart,1),{(0,)*N:(-1)%P})]
def ser(q):
 vn=names+['z0','z1'];out=[]
 for m,co in sorted(q.items(),reverse=True):
  fac=[]
  for nm,e in zip(vn,m):
   if e==1:fac.append(nm)
   elif e>1:fac.append(f'{nm}^{e}')
  mm='*'.join(fac);out.append(str(co) if not mm else (mm if co==1 else f'{co}*{mm}'))
 return '+'.join(out) if out else '0'
path=Path(args.output);lines=[','.join(names+['z0','z1']),str(P)]+[ser(q)+(',' if k<len(E)-1 else '') for k,q in enumerate(E)];path.write_text('\n'.join(lines)+'\n')
meta={'prime':P,'row':i,'col':j,'variables':11,'equations':11,'bytes':path.stat().st_size,'terms':[terms(q) for q in E],'degrees':[degree(q) for q in E]};path.with_suffix(path.suffix+'.meta.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta))
