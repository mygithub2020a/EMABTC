#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,time,heapq
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ

EF='13d2dab0a1488981913fe9017ed554a0f5435f5a98d78c1bd95cc667ca14086a'
EG='8dd5c088d0aca3f8fad1bae5c0e32aea72f9259d7dee61e5b722c3dd49b4039a'
EH='0420d4163dc3d24e4adba9dba14a4218ee8dc50b30fec41cbb50a15c74ad09be'

def sha(dd):
 t='\n'.join(','.join(map(str,k))+':'+str(dd[k]) for k in sorted(dd));return hashlib.sha256(t.encode()).hexdigest()

def mul(A,B):
 C={}
 for ma,ca in A.items():
  for mb,cb in B.items():
   m=tuple(x+y for x,y in zip(ma,mb));C[m]=C.get(m,0)+ca*cb
 return {m:c for m,c in C.items() if c}

def divexact(P,D):
 R=dict(P);Q={};lm=max(D);lc=D[lm];heap=[tuple(-x for x in m) for m in R];heapq.heapify(heap)
 while heap:
  nm=heapq.heappop(heap);m=tuple(-x for x in nm);a=R.get(m,0)
  if not a:continue
  if any(m[i]<lm[i] for i in range(len(m))) or a%lc:raise AssertionError(('nonzero remainder',m,a,lm,lc))
  qm=tuple(m[i]-lm[i] for i in range(len(m)));qc=a//lc;Q[qm]=Q.get(qm,0)+qc
  for dm,dc in D.items():
   x=tuple(qm[i]+dm[i] for i in range(len(m)));old=R.get(x,0);new=old-qc*dc
   if new:
    R[x]=new
    if not old:heapq.heappush(heap,tuple(-z for z in x))
   elif x in R:del R[x]
 return {m:c for m,c in Q.items() if c}

def main():
 t0=time.time();R,s,u,v,b,c=ring('s,u,v,b,c',ZZ,order='grevlex');one=R(1)
 pp=u**3+v**3;d=v**3-u**3;D4=one-4*pp+4*b*s*d
 Nnum=4*s**2*d+4*s**3*pp-s**3+b-4*b*s**3*pp-4*b*s**4*d
 Xn=-Nnum;Xd=s**3*D4;Yn=Xd-Xn;Z=b*Xd+Yn
 A=c*Yn+4*b*c*pp*Xd+4*b*c*Yn*s*d-2*b*Xn*D4
 Pd=8*c*Z;Qn=A*Xd+2*Xn*D4*Z;Qd=Pd*Xd;Tn=one-u**2*(one-s)**2
 BA=u**2*Xn**2+Tn*Yn**2;BH=u**2*s**2*Xn**2+Tn*Xd**2
 F=A**2*BA**3-Xn**6*Pd**2*u**6;G=Qn**2*BH**3-Xn**6*Qd**2*u**6
 Sn=u**2-v**2;Sd=4*u**2*v**2
 def eliminate_s(poly):
  deg=max(mon[0] for mon,_ in poly.terms());coeff=[R.zero]*(deg+1)
  for (es,eu,ev,eb,ec),co in poly.terms():coeff[es]+=R({(0,eu,ev,eb,ec):co})
  out=R.zero
  for k,ck in enumerate(coeff):
   if ck:out+=ck*Sn**k*Sd**(deg-k)
  return out
 def strip(poly):
  ts=poly.terms();mu=min(mon[1] for mon,_ in ts);mv=min(mon[2] for mon,_ in ts)
  return {(mon[1]-mu,mon[2]-mv,mon[3],mon[4]):int(co) for mon,co in ts}
 fhat=strip(eliminate_s(F));ghat=strip(eliminate_s(G));assert sha(fhat)==EF and sha(ghat)==EG
 D1={(2,0,0,0):1,(0,2,0,0):-1};D1=mul(mul(mul(D1,D1),mul(D1,D1)),mul(D1,D1));G6=divexact(ghat,D1)
 QD={(5,0,1,0):1,(3,2,1,0):-1,(2,3,1,0):-1,(0,5,1,0):1,(5,2,0,0):4,(2,5,0,0):4,(2,2,0,0):-1}
 H=divexact({m:z//(2**29) for m,z in G6.items()},mul(QD,QD));assert sha(H)==EH
 def to_text(P):
  E={}
  for (iu,iv,ib,ic),z in P.items():
   key=(ic,iu,iv);E[key]=(E.get(key,0)+z)%127
  terms=[]
  for (ic,iu,iv),z in sorted(((k,z) for k,z in E.items() if z),reverse=True):
   fac=[]
   if ic:fac.append('c' if ic==1 else f'c^{ic}')
   if iu:fac.append('u' if iu==1 else f'u^{iu}')
   if iv:fac.append('v' if iv==1 else f'v^{iv}')
   mon='*'.join(fac) if fac else '1';terms.append(mon if z==1 else f'{z}*{mon}')
  return '+'.join(terms),len(terms)
 Ft,nF=to_text(fhat);Ht,nH=to_text(H)
 sing='ring R=127,(c,u,v),dp;\npoly F='+Ft+';\npoly H='+Ht+';\nprint("BEGIN_RAWQ_FACTOR");\nprint("F_terms="+string(size(F)));\nprint("H_terms="+string(size(H)));\npoly Q=resultant(F,H,c);\nprint("Q_terms="+string(size(Q)));\nprint("Q_degree="+string(deg(Q)));\nlist L=factorize(Q,1);\nL;\nprint("END_RAWQ_FACTOR");\nquit;\n'
 Path('factor_rawQ.sing').write_text(sing)
 meta={'F_terms_exact':len(fhat),'H_terms_exact':len(H),'F_terms_b1_mod127':nF,'H_terms_b1_mod127':nH,'generator_elapsed':time.time()-t0,'F_hash':sha(fhat),'H_hash':sha(H)}
 Path('rawQ_meta.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
if __name__=='__main__':main()
