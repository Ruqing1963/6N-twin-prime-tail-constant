# 6N Twin-Prime Tail Constant (Part XI)

The first-principles origin of the Part VII tail constant K ≈ 0.105 — and the
reduction of the whole conditional theory to a single scale.

**Background.** Part VII wrote the right-centre survival as
P(N+d twin|ω) = K·∏_{q∈POOL} f_q(d,N), POOL={5..47}, with one empirical constant
K. It absorbed the primes q>47 and the relative-to-absolute normalisation.

**The decomposition.**

```
    K(d) = rho · tail_S(d) · B
```

- **rho** = twin-centre line density (Part VIII).
- **tail_S(d)** = S(d) / S_POOL(d), the part of the Hardy–Littlewood singular
  series S(d) carried by primes q>47.
- **B** = ∏_{q∈POOL} S_q(d) / ⟨f_q⟩, a pure-CRT basis-conversion constant,
  predicted analytically as **B ≈ 6.535** and d-independent.

Derivation: Part VIII gives the merged identity P_merged(d) = rho·S(d); the
ω-merged Part VII form is P_merged(d) = K·⟨∏_POOL f_q⟩. Equating,
K = rho·S(d)/⟨∏_POOL f_q⟩ = rho·tail_S(d)·B.

**Result (S₉ and S₁₀).** The ratio K/(rho·tail_S) is d-independent:
CV 0.31% (S₉), **0.16% (S₁₀)**, and equals the CRT prediction B=6.535 to within
0.4%. So K is **not a free parameter**.

**Decisive evidence — rho is the only scale.** From S₉ to S₁₀ the density drops
0.0199 → 0.0160 (factor 0.80), and K drops by the same factor (0.128 → 0.103)
while K/(rho·tail) stays at 6.55. The shell dependence of K is entirely that of
rho; tail_S(d) and B are shell-independent arithmetic.

**One scale for the whole theory.** Part VIII gave C₀ = 1/rho; this gives
K = rho·tail_S(d)·B. Both constants reduce to rho. The assembled gap preference

```
    r(d|ω) = S(d) · rho · P(N+d twin|ω),   P = rho·tail_S(d)·B·∏_POOL f_q
```

contains **no fitted constant**: S and tail_S are the Hardy–Littlewood singular
series, f_q and B are closed-form CRT, and rho is the twin-centre density (given,
to leading order, by the Hardy–Littlewood twin constant and the shell's mean
1/ln²). The single empirical input is rho — itself a known density, not a
constant of this construction.

> **Scope.** Experimental / computational number theory; d tested on S₉, S₁₀.
> The 0.2–0.4% excess of the measured ratio over B is the uniform-residue
> approximation in the q∤N branch (Part VII). No claim about the infinitude of
> twin primes or any k-tuple conjecture.

Part I: doi:10.5281/zenodo.20470367 · VII: doi:10.5281/zenodo.20518470 ·
VIII: doi:10.5281/zenodo.20519998

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── data/
│   └── K_decomp_data.csv    shell, d, K_meas, rho_tail, ratio, B_CRT  (S9 + S10)
├── code/
│   ├── K_decomp.py          measures K(d), forms K/(rho·tail_S), compares to the
│   │                        CRT prediction B; emits K_decomp_S{K}.csv
│   └── make_K_fig.py        builds the 3-panel figure from ../data
├── figures/                fig_paper11_K.{pdf,png}
└── paper/                  Chen_6N_Paper11.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# 1. Decompose K. Default S10 (~4 min: only twins + per-d CRT products).
#    Emits K_decomp_S{K}.csv.
python code/K_decomp.py            # S10
MAXK=9 python code/K_decomp.py     # S9

#    Concatenate both CSVs (keep one header) into data/K_decomp_data.csv
#    to feed the figure across shells.

# 2. Figure (reads ../data/K_decomp_data.csv).
cd code && python make_K_fig.py
```

### Conventions (same as Parts I–X)

- Twin centre N: 6N−1, 6N+1 both prime. dead(q) = {±6⁻¹ mod q}.
- S_q(d) = [#{r: r, r+d both q-safe}/q] / ((q−2)/q)²  (singular-series factor).
- f_q: Part VII CRT survival factor (q|N deterministic; q∤N admissible fraction).
- ⟨f_q⟩ = (1/|Adm_q|)·[d q-safe] + (1−1/|Adm_q|)·f_q^(q∤N), Adm_q = twin-admissible residues.
- rho = (#twin centres)/(#N in shell). tail_S(d) = S(d)/S_POOL(d), POOL={5..47}.
- Engine: deterministic interval-sieve primality; S₁₀ twin count 23,988,173 matches Part I.

## License

MIT — see `LICENSE`.
