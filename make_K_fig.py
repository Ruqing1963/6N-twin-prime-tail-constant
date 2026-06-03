#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 3-panel K-decomposition figure from ../data/K_decomp_data.csv (S9+S10).
Left: K(d) vs rho*tail*B both shells. Centre: ratio K/(rho*tail) flat at B.
Right: K and 1/C0 both scale with rho between shells (one scale).
Run K_decomp.py at MAXK=9 and MAXK=10 and concatenate their CSVs into
../data/K_decomp_data.csv (shell column distinguishes them) to regenerate.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/K_decomp_data.csv')))
def sub(sh):
    d=[int(r['d']) for r in rows if r['shell']==sh]
    K=[float(r['K_meas']) for r in rows if r['shell']==sh]
    rt=[float(r['rho_tail']) for r in rows if r['shell']==sh]
    ra=[float(r['ratio']) for r in rows if r['shell']==sh]
    return np.array(d),np.array(K),np.array(rt),np.array(ra)
B=float(rows[0]['B_CRT'])
fig,axes=plt.subplots(1,3,figsize=(17,4.8))
for sh,col,mk in [('S9','#185FA5','o'),('S10','#c0392b','s')]:
    d,K,rt,ra=sub(sh)
    if len(d)==0: continue
    o=np.argsort(d)
    axes[0].plot(d[o],K[o],mk,color=col,ms=9,label=f'{sh} measured $K(d)$',zorder=4)
    axes[0].plot(d[o],(rt*B)[o],'--',color=col,lw=1.5,alpha=.7,label=f'{sh} $\\rho\\cdot\\mathrm{{tail}}\\cdot B$',zorder=3)
axes[0].set_xlabel('centre-step $d$',fontsize=11); axes[0].set_ylabel('$K(d)$',fontsize=11)
axes[0].set_title('$K(d)=\\rho\\cdot\\mathrm{tail}_\\mathfrak{S}(d)\\cdot B$  (two shells)',fontsize=11)
axes[0].legend(fontsize=8.5); axes[0].grid(alpha=.25)
for sh,col,mk in [('S9','#185FA5','o'),('S10','#c0392b','s')]:
    d,K,rt,ra=sub(sh)
    if len(d)==0: continue
    o=np.argsort(d)
    axes[1].plot(d[o],ra[o],mk+'-',color=col,lw=1.6,ms=8,label=f'{sh}: $K/(\\rho\\,\\mathrm{{tail}})$')
axes[1].axhline(B,color='gray',ls=':',lw=1.4,label=f'CRT $B={B:.3f}$')
axes[1].set_ylim(6.3,6.8)
axes[1].set_xlabel('centre-step $d$',fontsize=11); axes[1].set_ylabel(r'$K/(\rho\cdot\mathrm{tail})$',fontsize=11)
axes[1].set_title('ratio is $d$-independent $=B$ (CV $\\leq0.31\\%$)',fontsize=11)
axes[1].legend(fontsize=9); axes[1].grid(alpha=.25)
labels=['$\\rho$','$1/C_0$','$K$ (mean)']
rho={'S9':0.019895,'S10':0.015992}; C0={'S9':50.06,'S10':62.75}
Kmean={sh:np.mean([float(r['K_meas']) for r in rows if r['shell']==sh]) for sh in ['S9','S10']}
s9=[rho['S9'],1/C0['S9'],Kmean['S9']]; s10=[rho['S10'],1/C0['S10'],Kmean['S10']]
x=np.arange(3); w=0.35
s10n=[s10[i]/s9[i] for i in range(3)]
axes[2].bar(x-w/2,[1,1,1],w,color='#185FA5',label='S9 (norm to 1)')
axes[2].bar(x+w/2,s10n,w,color='#c0392b',label='S10 / S9')
axes[2].axhline(rho['S10']/rho['S9'],color='gray',ls=':',lw=1.3,label=f'$\\rho_{{S10}}/\\rho_{{S9}}={rho["S10"]/rho["S9"]:.3f}$')
axes[2].set_xticks(x); axes[2].set_xticklabels(labels,fontsize=11)
axes[2].set_ylabel('shell ratio (S10/S9)',fontsize=10.5)
axes[2].set_title(r'$K$, $1/C_0$ both track $\rho$: one scale',fontsize=11)
axes[2].legend(fontsize=8.5); axes[2].grid(alpha=.25,axis='y'); axes[2].set_ylim(0,1.2)
plt.suptitle('First-principles decomposition of $K$ in $S_9$ and $S_{10}$:  $K=\\rho\\cdot\\mathrm{tail}_\\mathfrak{S}(d)\\cdot B$, $B$ a pure-CRT constant',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper11_K.pdf',bbox_inches='tight')
plt.savefig('fig_paper11_K.png',dpi=160,bbox_inches='tight')
print("figure saved")
