# SDMC structural-composite closure audit

Date: 2026-09-24

## Status

This note records a theory audit, not a new likelihood fit and not a claim of a unique microscopic completion.

The accepted manuscript-compatible covariant completion is

[
G_4(phi)=F(phi)/2,qquad
G_3(phi,X)=g_3(phi)X,qquad
G_2(phi,X)=k_1(phi)X+k_2(phi)X^2-V(phi).
]

For late266, the accepted linear-(G_3) free replay reproduces the parameterized target at
(|Delta H/H|_{max,z<100}=2.14	imes10^{-9}),
(|Delta F/F|_{max,z<100}=1.34	imes10^{-10}),
and No-Slip residual (3.26	imes10^{-10}).
The free solution remains stable with (D_{min}=0.05035) and
(c_{s,min}^2=0.05291).

The same accepted linear-(G_3) action gives the following exact likelihood values:

- Planck full Plik: 2784.6563603701597; fair residual vs optimized local021 = +3.2193316685643.
- DESI DR1 full shape: 330.6942517693205; fair residual = -1.2430105567822.
- Pantheon+: 1405.9357511694543; fair residual = -0.2789521529339.
- Union3: 28.59274220750376; fair residual = -0.1903977946931.
- DES-Y5: 1650.1821092143655; fair residual = -0.4532855004072.

Thus the accepted linear-(G_3) covariant completion preserves the late266 phenomenology:
DESI and all three SN sets close, while Planck remains the open gate.

## Active versus bare structural density

The action-level positive structural source is

[
ho_X propto H^2 F(1+alpha_M)-ho_m-ho_r.
]

This source is the *active* density.  It should not be identified directly with the
bare inverse-square SDMC density.  Instead,

[
ho_X=chi,ho_v,qquad
ho_v=ho_P S^{-2}.
]

Define

[
p_Xequiv-rac12rac{dlnho_X}{dln a},qquad
pequivrac{dln S}{dln a},qquad
A_chiequivrac{dlnchi}{dln a}.
]

Then exactly,

[
p_X=p-rac12A_chi.
]

The covariant action determines (p_X).  A second structural closure is required to
separate the bare structural evolution (p) from the activation history (A_chi).

## Chronological closure

The Part-III chronological hypothesis is

[
rac{S}{S_0}=rac{t_c}{t_0}.
]

It implies

[
p=rac{1}{Ht_c},qquad
Nequiv t_Pdot S=rac{t_P S_0}{t_0}=mathrm{constant}
]

through the classical branch.

For the accepted late266 covariant background and the working
(S_0=2.72	imes10^{61}),

[
t_0=13.6129780262 {m Gyr},qquad
N=3.41350846465.
]

At the present epoch,

[
p_0=1.03421565856,qquad
p_{X0}=0.03503108773,
]

hence

[
A_{chi0}=2(p_0-p_{X0})=1.99836914167.
]

This nearly reproduces the (A_chisimeq2) activation rate previously required for
a cosmological-constant-like active density, but here it follows from combining the
covariant active source with the chronological structural closure.

The action-derived active fraction is

[
Omega_{X0}=0.72374015096.
]

With the original inverse-square normalization (S_0=2.72	imes10^{61}),

[
Omega_{v,{m raw},0}=0.76904491838,
qquad
chi_0=0.94108956924.
]

## Radiation-matter bridge

The chronological closure reproduces the analytic Part-III radiation-matter bridge
extremely closely.  Using the actual equality redshift

[
z_{m eq}=3466.5021,
]

the difference between (p=1/(Ht)) from the covariant background and the analytic
radiation-matter formula over (100<z<10^6) is

[
max|Delta p|=3.42	imes10^{-4},qquad
{m RMS}(Delta p)=1.92	imes10^{-4}.
]

This is a strong internal consistency check.  The expected limiting behavior is recovered:
(p	o2) in radiation domination and (p	o3/2) in matter domination.

## Derived present mapping coefficients

Using the current late266 physical matter densities,

[
K_{m0}=2.1992619194	imes10^{20},
]

[
K_{p0}=1.1743136282	imes10^{20},
]

and

[
rac{K_{p0}}{K_{m0}}
=f_b^{1/3}
=0.53395806015.
]

With the chronological lapse,

[
lambda_0=K_{m0}^{-1/6}N^{-1}
=1.1923932901	imes10^{-4},
]

and the historical source-domain diagnostic is

[
lambda_0^3
=1.6953468762	imes10^{-12}.
]

These remain domain-projection diagnostics; they are not interpreted as literal
time-varying particle rest masses.

## Present N=3.37 calibration

If the legacy benchmark (N_0=3.37) is imposed exactly on the same covariant age,
the structural normalization required is

[
S_0=2.6853307310	imes10^{61},
]

only (1.275%) below the earlier (2.72	imes10^{61}) working value.  This gives

[
R_0=4.3401792207	imes10^{26} {m m},
]

[
K_{m0}=2.1712300065	imes10^{20},qquad
K_{p0}=1.1593457624	imes10^{20}.
]

Under this exact-(N_0) calibration,

[
lambda_0=1.2103727208	imes10^{-4},
qquad
lambda_0^3=1.7731986061	imes10^{-12},
]

which are within about (0.11%) and (0.41%), respectively, of the historical
Part-II values (1.209	imes10^{-4}) and (1.766	imes10^{-12}).

## Local-gravity proxy from the accepted action

For the same accepted linear-(G_3) action at (z=0),

[
F_0=1.02346657,qquad
alpha_{M0}=1.09024	imes10^{-3},
qquad
D_0=0.38666113,qquad
c_s^2(0)=0.138335.
]

The decoupling-limit canonical proxy gives

[
Q_{m time}=8.67	imes10^{-4},
]

and an unscreened simple conformal estimate

[
gamma_{m PPN}-1sim-3.0	imes10^{-6}.
]

The actual cubic coefficient gives a conditional Vainshtein scale of roughly
(2.4)--(9.0) pc for the Sun depending on time/radial normalization proxy, with
deep-screened scalar-force fraction around (4.3	imes10^{-15}) at 1 AU.
These are not yet a substitute for the full static spherical Horndeski solution.

The naive cosmological estimate is

[
dot G/Gsimeq-7.74	imes10^{-14} {m yr}^{-1}.
]

A full local solution is still required because Vainshtein screening of spatial
forces does not by itself determine the locally measured time variation of (G_N).

## Main conclusion

The accepted action now supports the following hierarchy without introducing
(N), (K_m), or (K_p) as independent fundamental couplings:

[
{G_2,G_3,G_4}	o ho_X,
]

[
ho_X=chiho_v,qquad
ho_v=ho_P S^{-2},
]

[
S/S_0=t_c/t_0quad	ext{(chronological closure)},
]

[
N=t_Pdot S,qquad
K_m=S/S_m,qquad
K_p=S/S_b.
]

The remaining fundamental question is whether the chronological closure
(Spropto t_c), or equivalently the required (chi) equation, can be derived
from the same local covariant action rather than imposed as an additional
structural condition.
