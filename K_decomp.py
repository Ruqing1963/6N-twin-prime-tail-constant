#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
First-principles decomposition of the Part VII tail constant K (S10).
Claim: K(d) = rho * S(d) / <prod_{POOL} f_q>_merged, equivalently
       K(d) = rho * tail_S(d) * B,
where rho = twin-centre line density (Part VIII), tail_S(d) = S(d)/S_POOL(d) is the
>47 part of the Hardy-Littlewood singular series, and
       B = prod_{q in POOL} [ S_q(d) / <f_q>_q ]   (pure CRT, d-independent)
is the POOL basis-conversion constant (predicted ~6.54). We measure K(d) on real
data for several d, form K/(rho*tail_S), and check it is a d-independent constant
matching the CRT prediction B. If so, K is not a free constant: it is the twin
density times the singular-series tail times a closed-form CRT factor.
Default S10. Requires: numpy.
"""
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
POOL=[5,7,11,13,17,19,23,29,31,37,41,43,47]
def dead(q): q=int(q); inv=pow(6,q-2,q); return {inv%q,(-inv)%q}
def primes_list(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(5,n+1) if s[i]]
def Sfac(d,q):
    dq=dead(q); safe=sum(1 for r in range(q) if r not in dq and (r+d)%q not in dq)/q
    return safe/((q-2)/q)**2
def S_full(d,qmax=5000):
    p=1.0
    for q in primes_list(qmax): p*=Sfac(d,q)
    return p
def S_pool(d):
    p=1.0
    for q in POOL: p*=Sfac(d,q)
    return p
def f_qnot(d,q):
    dq=dead(q); A=[r for r in range(q) if r not in dq and r!=0]
    return sum(1 for r in A if (r+d)%q not in dq)/len(A)
def B_predict(d):
    # pure CRT POOL basis-conversion: prod S_q / <f_q>_merged
    pr=1.0
    for q in POOL:
        dq=dead(q); adm=[r for r in range(q) if r not in dq]  # twin-center admissible (incl 0)
        pdiv=1/len(adm); dsafe=1.0 if (d%q) not in dq else 0.0
        mean_fq=pdiv*dsafe+(1-pdiv)*f_qnot(d,q)
        pr*=Sfac(d,q)/mean_fq
    return pr
# scan twins
Ns=[]; n=LO; import time; t0=time.time()
while n<=HI:
    nh=min(n+SEG,HI+1)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    Narr=np.arange(n,nh,dtype=np.int64)
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    Ns.append(Narr[np.nonzero(tw)[0]]); n=nh
N_arr=np.concatenate(Ns); Nrange=HI-LO+1; rho=len(N_arr)/Nrange
print(f"S{MAXK} twins {len(N_arr):,}; scan {time.time()-t0:.0f}s; rho={rho:.6f}")
def is_twin(vals):
    idx=np.searchsorted(N_arr,vals); idx=np.clip(idx,0,len(N_arr)-1)
    return N_arr[idx]==vals
def measure_K(d):
    prod=np.ones(len(N_arr))
    for q in POOL:
        res=N_arr%q; dq=dead(q); dsafe=(d%q) not in dq
        prod*=np.where(res==0,(1.0 if dsafe else 0.0),f_qnot(d,q))
    return is_twin(N_arr+d).mean()/prod.mean()
print(f"\nK(d) = rho * tail_S(d) * B,  B_CRT predicted (d-indep, ~6.54):")
print(f"{'d':>3}{'K meas':>10}{'rho*tail':>11}{'ratio=K/(rho*tail)':>20}{'B_CRT pred':>12}{'ratio/B':>9}")
ratios=[]
for d in [1,5,7,11,13,17,35]:
    K=measure_K(d); tailS=S_full(d)/S_pool(d); base=rho*tailS
    ratio=K/base; Bp=B_predict(d)
    ratios.append(ratio)
    print(f"{d:>3}{K:>10.5f}{base:>11.5f}{ratio:>20.4f}{Bp:>12.4f}{ratio/Bp:>9.4f}")
ratios=np.array(ratios)
print(f"\nratio K/(rho*tail): mean={ratios.mean():.4f}, CV={100*ratios.std()/ratios.mean():.2f}%")
print(f"=> if CV small and ratio ~ B_CRT, K is fully decomposed: K = rho * tail_S(d) * B,")
print(f"   with B a closed-form CRT constant. K is then NOT a free parameter.")

# ---- emit CSV (K_decomp_S{K}.csv) ----
import csv as _csv
with open(f'K_decomp_S{MAXK}.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['shell','d','K_meas','rho_tail','ratio','B_CRT'])
    _sh=f'S{MAXK}'
    for d in [1,5,7,11,13,17,35]:
        K=measure_K(d); tailS=S_full(d)/S_pool(d); base=rho*tailS
        _w.writerow([_sh,d,f'{K:.5f}',f'{base:.5f}',f'{K/base:.4f}',f'{B_predict(d):.4f}'])
print(f"\n[ok] wrote K_decomp_S{MAXK}.csv")
