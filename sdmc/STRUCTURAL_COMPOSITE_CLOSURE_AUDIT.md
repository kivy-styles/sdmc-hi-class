# SDMC structural-composite closure audit

Date: 2026-09-24

## Status

This is a theory audit, not a new likelihood fit and not a claim of a unique microscopic completion.

The accepted manuscript-compatible covariant completion is

\[
G_4(\phi)=F(\phi)/2,\qquad
G_3(\phi,X)=g_3(\phi)X,\qquad
G_2(\phi,X)=k_1(\phi)X+k_2(\phi)X^2-V(\phi).
\]

For late266, the accepted linear-\(G_3\) free replay reproduces the parameterized target at
\(|\Delta H/H|_{\max,z<100}=2.14\times10^{-9}\),
\(|\Delta F/F|_{\max,z<100}=1.34\times10^{-10}\),
and No-Slip residual \(3.26\times10^{-10}\).
The free solution remains stable with \(D_{\min}=0.05035\) and
\(c_{s,\min}^2=0.05291\).

The same accepted linear-\(G_3\) action gives:

- Planck full Plik: 2784.6563603701597; fair residual vs optimized local021 \(=+3.2193316685643\).
- DESI DR1 full shape: 330.6942517693205; fair residual \(=-1.2430105567822\).
- Pantheon+: 1405.9357511694543; fair residual \(=-0.2789521529339\).
- Union3: 28.59274220750376; fair residual \(=-0.1903977946931\).
- DES-Y5: 1650.1821092143655; fair residual \(=-0.4532855004072\).

Thus the accepted linear-\(G_3\) covariant completion preserves the late266 phenomenology:
DESI and all three SN sets close, while Planck remains the open gate.

## Active versus bare structural density

The action-level positive structural source is

\[
\rho_X \propto H^2F(1+\alpha_M)-\rho_m-\rho_r.
\]

This is the active density. It should not be identified directly with the bare inverse-square SDMC density. Instead,

\[
\rho_X=\chi\,\rho_v,\qquad
\rho_v=\rho_P S^{-2}.
\]

Define

\[
p_X\equiv-\frac12\frac{d\ln\rho_X}{d\ln a},\qquad
p\equiv\frac{d\ln S}{d\ln a},\qquad
A_\chi\equiv\frac{d\ln\chi}{d\ln a}.
\]

Then exactly,

\[
p_X=p-\frac12A_\chi.
\]

The covariant action determines \(p_X\). A structural closure is required to separate the bare structural evolution \(p\) from the activation history \(A_\chi\).

## Chronological closure

The Part-III chronological hypothesis is

\[
\frac{S}{S_0}=\frac{t_c}{t_0}.
\]

It implies

\[
p=\frac{1}{Ht_c},\qquad
N\equiv t_P\dot S=\frac{t_PS_0}{t_0}=\mathrm{constant}
\]

through the classical branch.

For the accepted late266 covariant background and \(S_0=2.72\times10^{61}\),

\[
t_0=13.6129780262\ {\rm Gyr},\qquad
N=3.41350846465.
\]

At the present epoch,

\[
p_0=1.03421565856,\qquad
p_{X0}=0.03503108773.
\]

The present geometric coefficient is

\[
\Xi_0\equiv\frac{H_0R_0}{c}=t_PS_0H_0=3.30057704734,
\]

and therefore

\[
N_0=p_0\Xi_0=3.41350846465.
\]

The historical numerical identification \(N_0\simeq\Xi_0\simeq3.37\) was the synchronized \(p=1\) approximation. It need not be exact for the current late266 state. The earlier \(\Xi_0\simeq3.37\) also used the older \(H_0\simeq70.85\) working benchmark; the current late266 \(H_0\) lowers \(\Xi_0\) to 3.30058 before the \(p_0\) factor is applied.

The activation rate is then

\[
A_{\chi0}=2(p_0-p_{X0})=1.99836914167.
\]

This nearly reproduces the \(A_\chi\simeq2\) rate previously required for a cosmological-constant-like active density, now as a consequence of the covariant active source plus the chronological structural closure.

The action-derived present active fraction is

\[
\Omega_{X0}=0.72374015096.
\]

With \(S_0=2.72\times10^{61}\),

\[
\Omega_{v,{\rm raw},0}=0.76904491838,\qquad
\chi_0=0.94108956924.
\]

Across the computed classical history,

\[
0.00400015\le\chi\le0.94108957.
\]

The shallow minimum occurs near \(z\simeq12.69\); \(A_\chi\) reaches a minimum \(-0.5019\) near \(z\simeq18.31\), then rises and peaks near \(z\simeq3.01\). Thus the current closure is bounded, but not strictly monotonic over the whole classical history.

## Radiation-matter bridge

Using the actual equality redshift

\[
z_{\rm eq}=3466.5021,
\]

the chronological closure reproduces the analytic Part-III radiation-matter bridge over \(100<z<10^6\) with

\[
\max|\Delta p|=3.42\times10^{-4},\qquad
{\rm RMS}(\Delta p)=1.92\times10^{-4}.
\]

The expected limits are recovered:

\[
p\to2 \quad {\rm (radiation)},\qquad
p\to\frac32 \quad {\rm (matter)}.
\]

This is a strong internal consistency test of \(S\propto t_c\) on the classical tracking branches.

## Derived present mapping coefficients

Using the current late266 physical matter densities,

\[
K_{m0}=2.1992619194\times10^{20},
\]

\[
K_{p0}=1.1743136282\times10^{20},
\]

and

\[
\frac{K_{p0}}{K_{m0}}
=f_b^{1/3}
=0.53395806015.
\]

With the chronological lapse,

\[
\lambda_0=K_{m0}^{-1/6}N^{-1}
=1.1923932901\times10^{-4},
\]

and

\[
\lambda_0^3
=1.6953468762\times10^{-12}.
\]

These remain domain-projection diagnostics, not literal time-varying particle rest masses.

Compared with the historical Part-II anchors, the current values differ by about:
\(N_0:+1.29\%\), \(K_{m0}:+0.42\%\), \(K_{p0}:-0.48\%\),
\(\lambda_0:-1.37\%\), and \(\lambda_0^3:-4.00\%\).

If the legacy \(N_0=3.37\) benchmark is imposed exactly on the current covariant age, the required normalization is

\[
S_0=2.6853307310\times10^{61},
\]

only \(1.275\%\) below the old \(2.72\times10^{61}\) working value, giving

\[
R_0=4.3401792207\times10^{26}\ {\rm m},
\]

\[
K_{m0}=2.1712300065\times10^{20},\qquad
K_{p0}=1.1593457624\times10^{20},
\]

\[
\lambda_0=1.2103727208\times10^{-4},\qquad
\lambda_0^3=1.7731986061\times10^{-12}.
\]

Those last two are within about \(0.11\%\) and \(0.41\%\) of the historical Part-II values.

## Local-gravity proxy from the accepted action

For the same accepted linear-\(G_3\) action at \(z=0\),

\[
F_0=1.02346657,\qquad
\alpha_{M0}=1.09024\times10^{-3},
\]

\[
D_0=0.38666113,\qquad
c_s^2(0)=0.138335.
\]

The decoupling-limit canonical proxy gives

\[
Q_{\rm time}=8.67\times10^{-4},
\]

and an unscreened simple conformal estimate

\[
\gamma_{\rm PPN}-1\sim-3.0\times10^{-6}.
\]

The reconstructed cubic coefficient gives a conditional Vainshtein scale of roughly \(2.4\)--\(9.0\) pc for the Sun depending on time/radial normalization proxy, with a deep-screened scalar-force fraction around \(4.3\times10^{-15}\) at 1 AU. These are not yet a substitute for the full static spherical Horndeski solution.

The naive cosmological estimate is

\[
\dot G/G\simeq-7.74\times10^{-14}\ {\rm yr}^{-1}.
\]

A full local solution is still required because Vainshtein screening of spatial forces does not by itself determine the locally measured time variation of \(G_N\).

## Main conclusion

The accepted action now supports the hierarchy

\[
\{G_2,G_3,G_4\}\to \rho_X,
\]

\[
\rho_X=\chi\rho_v,\qquad
\rho_v=\rho_P S^{-2},
\]

\[
S/S_0=t_c/t_0\quad\text{(chronological closure)},
\]

\[
N=t_P\dot S,\qquad
K_m=S/S_m,\qquad
K_p=S/S_b.
\]

The remaining fundamental question is whether the chronological closure \(S\propto t_c\), or equivalently the split between \(\rho_X\) and \(\chi\rho_v\), can be derived from the same local covariant action rather than imposed as an additional structural condition.
