# SDMC structural-clock Lagrangian continuation

Date: 2026-09-24

## Purpose

This continuation starts from the accepted manuscript-compatible covariant action

\[
\mathcal S=\int d^4x\sqrt{-g}
\left[G_2(\phi,X)-G_3(\phi,X)\Box\phi+G_4(\phi)R\right]
+S_m+S_r,
\]

with

\[
G_4=F/2,\qquad
G_3=g(\phi)X,\qquad
G_2=k_1(\phi)X+k_2(\phi)X^2-V(\phi),
\]

and the on-trajectory No-Slip identity

\[
XG_{3X}=Xg=-\frac12F_{,\phi}.
\]

The goal is to place the SDMC structural clock directly inside the same
covariant scalar degree of freedom, without introducing a second metric or an
additional propagating activation field.

## 1. Logarithmic structural-clock field

On the accepted monotonic cosmological branch define

\[
\psi\equiv\ln\frac{S}{S_0}.
\]

For the chronological branch,

\[
\psi=\ln\frac{t}{t_0}.
\]

Let \(\phi=f(\psi)\) and

\[
A(\psi)\equiv\frac{d\phi}{d\psi}.
\]

Then, with

\[
Y\equiv-\frac12\nabla_\mu\psi\nabla^\mu\psi,
\qquad X=A^2Y,
\]

the accepted action is closed under the field redefinition:

\[
\widetilde G_4=\frac{\widetilde F}{2},
\qquad
\widetilde G_3=\widetilde g\,Y,
\qquad
\widetilde G_2=\widetilde k_1Y+\widetilde k_2Y^2-\widetilde V,
\]

with

\[
\widetilde F=F,
\qquad
\widetilde g=A^3g,
\qquad
\widetilde k_1=A^2k_1,
\]

\[
\widetilde k_2=A^4k_2+2A^2A_{,\psi}g,
\qquad
\widetilde V=V.
\]

The No-Slip identity is invariant:

\[
Y\widetilde G_{3Y}
=Y\widetilde g
=-\frac12\widetilde F_{,\psi}.
\]

Thus the structural clock can be used as the scalar coordinate of the same
Horndeski theory; this is a field redefinition, not a new degree of freedom.

## 2. Direct structural field and constant-kinetic form

Now define the nonlogarithmic structural field

\[
\sigma\equiv e^\psi=\frac{S}{S_0}.
\]

Let

\[
Z\equiv-\frac12\nabla_\mu\sigma\nabla^\mu\sigma.
\]

Since \(d\psi/d\sigma=1/\sigma\), the action again retains the same
linear-\(G_3\), quadratic-\(G_2\) form,

\[
\widehat G_4=\frac{\widehat F}{2},
\qquad
\widehat G_3=\widehat g\,Z,
\qquad
\widehat G_2=\widehat k_1Z+\widehat k_2Z^2-\widehat V,
\]

where

\[
\widehat F=\widetilde F,
\qquad
\widehat g=\frac{\widetilde g}{\sigma^3},
\qquad
\widehat k_1=\frac{\widetilde k_1}{\sigma^2},
\]

\[
\widehat k_2
=\frac{\widetilde k_2-2\widetilde g}{\sigma^4},
\qquad
\widehat V=\widetilde V.
\]

The transformed No-Slip identity remains

\[
Z\widehat G_{3Z}
=Z\widehat g
=-\frac12\widehat F_{,\sigma}.
\]

For a homogeneous FLRW background,

\[
Z=\frac12\dot\sigma^2.
\]

The chronological SDMC trajectory

\[
\sigma=\frac{t}{t_0}
\]

is therefore exactly the constant-kinetic orbit

\[
\dot\sigma=\frac1{t_0},
\qquad
Z=Z_\star\equiv\frac1{2t_0^2}.
\]

This converts the previous chronological closure into a precise dynamical
question: does the structural-field action admit and select the constant-\(Z\)
orbit?

## 3. Structural lapse as the scalar kinetic invariant

Because \(S=S_0\sigma\),

\[
\mathcal N\equiv t_P\dot S
=t_PS_0\dot\sigma.
\]

Hence

\[
\boxed{\mathcal N=t_PS_0\sqrt{2Z}}
\]

on the expanding branch, or equivalently

\[
\boxed{Z=\frac{\mathcal N^2}{2t_P^2S_0^2}}.
\]

Therefore a constant structural lapse is not an extra independent coupling in
the structural-field representation. It is equivalent to a constant
homogeneous kinetic density of the same scalar already present in the accepted
covariant action.

The stretch-rate index becomes

\[
p\equiv\frac{d\ln S}{d\ln a}
=\frac{\dot\sigma}{H\sigma}
=\frac{\sqrt{2Z}}{H\sigma}.
\]

For the chronological orbit this reduces to

\[
p=\frac1{Ht}.
\]

## 4. Matter and baryon mappers from the same field

For conserved pressureless matter, \(\rho_m\propto a^{-3}\), so the matter
reference scale obeys \(S_m/S_{m0}=a\). Therefore

\[
\frac{K_m}{K_{m0}}
=\frac{S/S_0}{a}
=\frac{\sigma}{a},
\]

and

\[
\frac{d\ln K_m}{d\ln a}=p-1.
\]

For a constant baryon fraction,

\[
\frac{K_p}{K_{p0}}=\frac{\sigma}{a},
\qquad
\frac{K_p}{K_m}=f_b^{1/3}.
\]

Thus \(\mathcal N,K_m,K_p\) are all composites of the same structural scalar
and conserved matter densities.

## 5. Activation factor as an action-derived composite

For the structural-field action, the positive homogeneous active source has the
same form as in the accepted linear-\(G_3\) reconstruction,

\[
\rho_X=
\widehat k_1 Z+3\widehat k_2 Z^2+\widehat V
+6H\dot\sigma Z\widehat g
-2Z^2\widehat g_{,\sigma}
-3H\dot{\widehat F}.
\]

The bare inverse-square SDMC density is

\[
\rho_v=\frac{\rho_P}{S_0^2\sigma^2}.
\]

Hence the activation factor can be written without introducing a new field:

\[
\boxed{
\chi(\sigma,Z,H)
=\frac{\rho_X}{\rho_v}
=\frac{S_0^2\sigma^2}{\rho_P}\,\rho_X
}.
\]

Its logarithmic rate is

\[
A_\chi
=\frac{1}{H}\frac{d\ln\chi}{dt}
=\frac{1}{H}\frac{d\ln\rho_X}{dt}+2p
=2(p-p_X),
\]

where

\[
p_X\equiv-\frac1{2H}\frac{d\ln\rho_X}{dt}.
\]

This recovers the previously derived active/bare slope relation, but now
expresses \(\chi\) explicitly as a composite of the structural-field action
and its background solution.

## 6. Density route to the structural lapse

Define

\[
\Omega_X\equiv\frac{8\pi G\rho_X}{3H^2},
\qquad
\Xi_P^2\equiv\frac{8\pi G}{3}\rho_Pt_P^2.
\]

Using \(\rho_X=\chi\rho_P/S^2\) gives

\[
Ht_PS
=\Xi_P\sqrt{\frac{\chi}{\Omega_X}}.
\]

Since \(\mathcal N=t_PHS\,p\),

\[
\boxed{
\mathcal N
=p\,\Xi_P\sqrt{\frac{\chi}{\Omega_X}}
}.
\]

This is the density representation of the same lapse. It is algebraically
identical to \(\mathcal N=t_P\dot S\) once the active/bare split is imposed.

## 7. Constant-kinetic equation from the homogeneous action

For the structural field, write

\[
\widehat G_2=k_1(\sigma)Z+k_2(\sigma)Z^2-V(\sigma),
\qquad
\widehat G_3=g(\sigma)Z,
\qquad
\widehat G_4=F(\sigma)/2.
\]

On homogeneous FLRW, after removing total derivatives, the scalar-dependent
minisuperspace Lagrangian per unit \(a^3\) can be written

\[
\ell=
\frac12k_1v^2
+Bv^4
+Hg v^3
-V
-3HF_{,\sigma}v
-3FH^2,
\]

where

\[
v\equiv\dot\sigma,
\qquad
B\equiv\frac14k_2-\frac16g_{,\sigma}.
\]

The exact homogeneous scalar equation is

\[
\frac{d}{dt}\left(\ell_{,v}\right)
+3H\ell_{,v}
-\ell_{,\sigma}=0.
\]

For the chronological constant-kinetic orbit \(v=v_\star=1/t_0\) and
\(\dot v=0\), this reduces to

\[
\boxed{
\begin{aligned}
0={}&
V_{,\sigma}
+\frac12 k_{1,\sigma}v_\star^2
+3B_{,\sigma}v_\star^4
+3Hk_1v_\star
+12HBv_\star^3 \\
&+2Hg_{,\sigma}v_\star^3
+3(\dot H+3H^2)gv_\star^2
-3(\dot H+2H^2)F_{,\sigma}.
\end{aligned}
}
\]

This is the action-level constant-kinetic closure condition. The remaining
question is not whether the chronological trajectory can be represented by the
action—the field-redefined free replay already shows that it can—but whether
this orbit is dynamically selected rather than requiring finely tuned scalar
initial data.

## 8. New GitHub attractor audit

Branch: sdmc-structural-clock-attractor-audit

The new audit leaves \(G_2,G_3,G_4\) unchanged and perturbs only the initial
structural scalar position and velocity. It measures

\[
q\equiv\frac{d\psi}{d\ln t},
\]

for which the chronological solution requires \(q=1\), together with

\[
\Delta_{\rm chrono}
=(\psi-\psi_0)-\ln(t/t_0),
\]

and the normalized structural speed

\[
\frac{u}{u_0},
\qquad
u\equiv\frac{d(e^\psi)}{dt}.
\]

If small perturbations decay while the background returns to the nominal
trajectory, the chronological orbit has a local basin of attraction. If they
do not, then the current reconstruction demonstrates an exact solution but
not a dynamically selected closure.

## Scientific status

The structural-field rewrite is an exact field redefinition of the accepted
covariant action on the monotonic cosmological branch. It introduces no new
propagating degree of freedom. The identities for \(\mathcal N,K_m,K_p\) and
\(\chi\) then follow from the definitions stated above.

What is not yet established is that the constant-\(Z\) trajectory is a unique
or generic prediction of an independently specified microscopic action. The
new attractor audit is designed to test the first necessary condition:
dynamical stability of that trajectory under scalar initial-condition
perturbations.


## 9. Epoch-dependent matter and baryon mappers

The recent chronological closure sharpens an important point that was only
conditional in Part II.  The mapper need not be constant at all epochs.

For a general matter reference density,

\[
S_m=\left(\frac{\rho_{m,P}}{\rho_m}\right)^{1/3},
\qquad
K_m=\frac{S}{S_m},
\]

hence directly

\[
\boxed{
K_m(S,\rho_m)
=
S\left(\frac{\rho_m}{\rho_{m,P}}\right)^{1/3}
}.
\]

With the working common primordial reference \(\rho_{m,P}=\rho_P\),

\[
\boxed{
K_m
=
S\left(\frac{\rho_m}{\rho_P}\right)^{1/3}
}.
\]

Combining this with the inverse-square bare density
\(\rho_v=\rho_P/S^2\) gives the equivalent density-only form

\[
\boxed{
K_m=
\left(
\frac{\rho_P\rho_m^2}{\rho_v^3}
\right)^{1/6}
}.
\]

For conserved pressureless matter, \(\rho_m=\rho_{m0}a^{-3}\), so

\[
\boxed{
\frac{K_m(a)}{K_{m0}}
=
\frac{S(a)/S_0}{a}
}
\]

and therefore

\[
\boxed{
\frac{d\ln K_m}{d\ln a}=p-1.
}
\]

Likewise, with constant baryon fraction,

\[
K_p=K_m f_b^{1/3},
\]

so

\[
\boxed{
\frac{K_p(a)}{K_{p0}}
=
\frac{K_m(a)}{K_{m0}}
=
\frac{S(a)/S_0}{a}
},
\qquad
\boxed{
\frac{K_p}{K_m}=f_b^{1/3}.
}
\]

Thus the baryon and total-matter mappers have the same fractional epoch
evolution; their ratio remains constant.

This also resolves the older Part-II conditional statement about an
epoch-invariant \(K_m\).  Part II showed that \(K_m\) is constant if the global
and matter radii have identical fractional expansion rates.  The covariant
chronological closure instead gives

\[
\frac{\dot K_m}{K_m}=H(p-1).
\]

Hence the constant-\(K_m\) result is recovered only in the synchronized mature
limit \(p\to1\).  During radiation and matter tracking, \(p>1\), so \(K_m\) and
\(K_p\) evolve.

### Last-scattering value in the accepted late266 closure

The structural-closure artifact samples the accepted covariant history near
last scattering and gives, at \(z=1091.31\),

\[
\frac{S_*}{S_0}=2.69548\times10^{-5},
\qquad
\frac{K_{m,*}}{K_{m0}}
=
\frac{K_{p,*}}{K_{p0}}
=0.02944299.
\]

Interpolating the same free covariant background to \(z_*=1090\) gives

\[
\frac{S_*}{S_0}=2.70096\times10^{-5},
\qquad
\frac{K_{m,*}}{K_{m0}}
=
\frac{K_{p,*}}{K_{p0}}
=0.02946748.
\]

Using the current action-derived present anchors,

\[
K_{m0}=2.19926192\times10^{20},
\qquad
K_{p0}=1.17431363\times10^{20},
\]

therefore gives

\[
\boxed{
K_m(z_*=1090)\simeq6.481\times10^{18}
}
\]

and

\[
\boxed{
K_p(z_*=1090)\simeq3.460\times10^{18}.
}
\]

The ratio is unchanged,

\[
\frac{K_p(z_*)}{K_m(z_*)}
=0.53395806=f_b^{1/3}.
\]

This is the epoch-dependent \(K_p\) value that must be distinguished from the
present-day anchor \(K_{p0}\simeq1.1743\times10^{20}\).  It should also be kept
separate from the historical acoustic use of \(K_p\): once the photon-baryon
sound horizon has already been computed in the local thermal domain, the
absolute mapper must not automatically be applied to that ruler a second time.

## 10. Source-map conservation constraint

The older source hierarchy is

\[
M_L=M_GK_m^{-1/2}{\cal N}^{-3},
\qquad
\lambda=K_m^{-1/6}{\cal N}^{-1}.
\]

The accepted action has universal Jordan-frame matter conservation.  Therefore
this relation cannot mean that the rest mass of a fixed collection of ordinary
particles literally changes with epoch.

With the chronological closure \({\cal N}=\) constant,

\[
\frac{d}{d\ln a}
\ln\left(K_m^{-1/2}{\cal N}^{-3}\right)
=
-\frac12(p-1).
\]

The factor varies during radiation and matter tracking and freezes only as
\(p\to1\).  Its viable interpretation is therefore a map between coarse-grained
source representations in different structural domains, not literal particle
mass evolution.

## 11. Placement of the activation features

The reconstructed activation history is not independent of the already fitted
background transitions.  In the accepted late266 solution the short negative
activation-rate precursor near \(z\simeq18.3\) and the minimum of \(\chi\) near
\(z\simeq12.7\) bracket the tracker handoff centered at \(z_t\simeq16.17\).
The maximum activation rate near \(z\simeq3.01\) lies within the effective
Planck-mass transition centered at \(z_c\simeq3.43\).

This favors treating \(\chi\) as a composite diagnostic of the existing early
handoff plus late scalar-tensor response, rather than adding a new propagating
activation field solely to reproduce the reconstructed \(\chi(a)\).


## 12. Linearized stability of the constant-kinetic orbit

The minisuperspace equation can be written in a form that isolates the
acceleration of the structural field.  With

\[
v\equiv\dot\sigma,
\qquad
B(\sigma)\equiv\frac14k_2-\frac16g_{,\sigma},
\]

the reduced scalar Lagrangian is

\[
\ell=
\frac12k_1v^2+Bv^4+Hgv^3-V
-3HF_{,\sigma}v-3FH^2.
\]

Define

\[
\boxed{
{\cal K}_{\rm hom}\equiv
\frac{\partial^2\ell}{\partial v^2}
=
k_1+12Bv^2+6Hgv
}.
\]

The homogeneous scalar equation can then be arranged exactly as

\[
\boxed{
{\cal K}_{\rm hom}\,\dot v+{\cal C}(\sigma,v;t)=0
}
\]

with

\[
\boxed{
\begin{aligned}
{\cal C}={}&
V_{,\sigma}
+\frac12 k_{1,\sigma}v^2
+3B_{,\sigma}v^4
+3Hk_1v
+12HBv^3\\
&+2Hg_{,\sigma}v^3
+3(\dot H+3H^2)gv^2
-3(\dot H+2H^2)F_{,\sigma}.
\end{aligned}
}
\]

The chronological orbit is the constant-velocity solution

\[
v=v_\star=\frac1{t_0},
\qquad
{\cal C}(\bar\sigma,v_\star;t)=0.
\]

For a scalar-only local stability diagnostic, temporarily holding the metric
background on the nominal solution, write

\[
\sigma=\bar\sigma+\delta\sigma,
\qquad
v=v_\star+\delta v.
\]

Since \(\dot v_\star=0\), the linearized system is

\[
\frac{d}{dt}
\begin{pmatrix}
\delta\sigma\\
\delta v
\end{pmatrix}
=
\begin{pmatrix}
0 & 1\\
-{\cal M}_\star^2 & -\Gamma_\star
\end{pmatrix}
\begin{pmatrix}
\delta\sigma\\
\delta v
\end{pmatrix},
\]

where

\[
\boxed{
\Gamma_\star=
\frac{{\cal C}_{,v}}{{\cal K}_{\rm hom}}
\bigg|_\star,
\qquad
{\cal M}_\star^2=
\frac{{\cal C}_{,\sigma}}{{\cal K}_{\rm hom}}
\bigg|_\star.
}
\]

If the coefficients are regarded as locally constant, the instantaneous
eigenvalues are

\[
\boxed{
\lambda_\pm=
\frac{-\Gamma_\star
\pm\sqrt{\Gamma_\star^2-4{\cal M}_\star^2}}{2}.
}
\]

Positive effective damping and a non-tachyonic local restoring term,

\[
\Gamma_\star>0,
\qquad
{\cal M}_\star^2>0,
\]

are sufficient local indicators in this frozen-background approximation.
They are not the final cosmological attractor criterion because the scalar
perturbation backreacts on \(H\), \(F\), and the matter fractions.

That is why the GitHub basin test perturbs the actual hi_class scalar initial
conditions and re-integrates the full homogeneous system rather than declaring
the orbit stable from the reduced equation alone.

A convenient numerical variable is

\[
q\equiv\frac{d\psi}{d\ln t}
=t\dot\psi
=\frac{t\dot\sigma}{\sigma}.
\]

Since

\[
p=\frac{\dot\sigma}{H\sigma},
\]

we also have

\[
\boxed{q=Htp}.
\]

The chronological solution has

\[
q=1,
\]

so convergence of \(q\to1\), together with decay of the phase-space
separation from the nominal solution, is a direct numerical attractor test.

## 13. Covariant composite hierarchy in the structural field

On the reconstructed monotonic cosmological branch
\(\phi=f(\psi)\) and \(a=e^\phi\).  Therefore the matter mapper can be written
directly as a composite of the same scalar coordinate,

\[
\boxed{
\frac{K_m}{K_{m0}}
=e^{\psi-\phi}
=e^{\psi-f(\psi)}
}
\]

and similarly

\[
\boxed{
\frac{K_p}{K_{p0}}
=e^{\psi-f(\psi)}.
}
\]

The lapse is

\[
\boxed{
{\cal N}
=t_PS_0e^\psi\dot\psi
=t_PS_0\dot\sigma.
}
\]

Thus the structural hierarchy may be summarized as

\[
\boxed{
\{G_2,G_3,G_4\}
\longrightarrow
\{\sigma,Z,\rho_X\}
\longrightarrow
\{\chi,{\cal N},K_m,K_p\},
}
\]

subject to the bare-density identification

\[
\rho_v=\frac{\rho_P}{S_0^2\sigma^2}.
\]

This is stronger than treating \({\cal N}\), \(K_m\), and \(K_p\) as separate
phenomenological parameters: once the structural scalar trajectory is fixed,
all three are derived composites.  The unresolved step is whether an
independently motivated microscopic form of the same covariant functions
selects the constant-\(Z\) orbit without reconstruction from the desired
background.


## 14. Shift-current route toward a dynamical chronological attractor

The reduced homogeneous equation admits a useful current form. Define

\[
\boxed{
J_\sigma\equiv\ell_{,v}
=
k_1v+4Bv^3+3Hgv^2-3HF_{,\sigma}.
}
\]

Then

\[
\boxed{
\dot J_\sigma+3HJ_\sigma=\ell_{,\sigma},
}
\]

where

\[
\ell_{,\sigma}
=
\frac12k_{1,\sigma}v^2
+B_{,\sigma}v^4
+Hg_{,\sigma}v^3
-V_{,\sigma}
-3HF_{,\sigma\sigma}v
-3H^2F_{,\sigma}.
\]

This gives a possible microscopic route to the chronological closure. If the
mature structural action approaches an approximately shift-symmetric regime,

\[
k_{1,\sigma},\ B_{,\sigma},\ g_{,\sigma},\ V_{,\sigma},\
F_{,\sigma},\ F_{,\sigma\sigma}\longrightarrow0,
\]

then

\[
\ell_{,\sigma}\longrightarrow0
\]

and the current obeys

\[
\dot J_\sigma+3HJ_\sigma\simeq0,
\qquad
J_\sigma\propto a^{-3}.
\]

Expansion therefore drives the system toward an algebraic current root,

\[
\boxed{
4Bv_\star^3+3Hgv_\star^2+k_1v_\star-3HF_{,\sigma}=0.
}
\]

In the further asymptotic limit in which \(H\to0\) and
\(F_{,\sigma}\to0\), the nonzero root becomes

\[
\boxed{
v_\star^2=-\frac{k_{1,\infty}}{4B_\infty}.
}
\]

Because \(v_\star=\dot\sigma\), such a nonzero constant root gives

\[
\sigma\propto t,\qquad
S\propto t,\qquad
{\cal N}=t_PS_0v_\star=\text{constant}.
\]

The corresponding lapse would be predicted directly by the asymptotic
Lagrangian coefficients,

\[
\boxed{
{\cal N}_\infty
=
t_PS_0
\sqrt{-\frac{k_{1,\infty}}{4B_\infty}}
}.
\]

At the same reduced level,

\[
{\cal K}_{{\rm hom},\star}
=
k_{1,\infty}+12B_\infty v_\star^2
=
-2k_{1,\infty}.
\]

Thus a real nonzero root with positive homogeneous kinetic curvature requires

\[
\boxed{
k_{1,\infty}<0,\qquad B_\infty>0
}
\]

in this asymptotic approximation. This is a ghost-condensate-like
constant-velocity mechanism, but the statement is deliberately limited to the
homogeneous reduced dynamics: the full Horndeski no-ghost and sound-speed
conditions remain the \(D>0\) and \(c_s^2>0\) tests already used in the
covariant runs.

This route is important because it changes the outstanding question from
"why impose \(S\propto t\)?" to a concrete action-level test: do the reconstructed
functions approach a regime in which the scalar current is Hubble-diluted
toward a stable nonzero kinetic root?

## 15. Reduced numerical stability screen of the accepted action

A new repository audit named run_structural_clock_linear_stability_audit.py
evaluates the quantities in Section 12 directly from the accepted free
linear-\(G_3\) background after the two structural field redefinitions.

A local reproduction on the accepted artifact gives, over \(z<100\),

\[
{\cal K}_{\rm hom}>0,
\]

with

\[
0.36998\lesssim{\cal K}_{\rm hom}\lesssim1.5042\times10^4,
\]

and

\[
0.1802\lesssim\frac{{\cal M}_\star^2}{H^2}\lesssim1.6913.
\]

Thus the reduced restoring term remains positive throughout this interval.
The effective damping changes sign once, near

\[
\boxed{z\simeq57.46.}
\]

For \(z\lesssim57.5\),

\[
\Gamma_\star>0,
\]

and the instantaneous reduced eigenvalues have negative real parts. At the
present epoch the diagnostic gives approximately

\[
\frac{\Gamma_\star}{H_0}\simeq2.1133,
\qquad
\frac{{\cal M}_\star^2}{H_0^2}\simeq0.18777,
\]

with

\[
\frac{\lambda_+}{H_0}\simeq-0.09294,
\qquad
\frac{\lambda_-}{H_0}\simeq-2.0204.
\]

This is a positive late-time indication: in the frozen-background scalar
subsystem the chronological constant-kinetic orbit is locally damped today.

The same reduced test is not uniformly attractive in the early universe. Near
last scattering it gives roughly

\[
\frac{\Gamma_\star}{H}\simeq-0.240,
\qquad
\frac{{\cal M}_\star^2}{H^2}\simeq2.105,
\]

so the local eigenvalues have a small positive real part,
approximately \(0.120H\). This does not by itself reject the chronological
history, because the frozen-background approximation omits the simultaneous
response of the metric, matter fractions and Planck mass. It does show that a
claim of an all-epoch scalar-only attractor would be too strong.

The decisive next test remains the full hi_class basin integration with
perturbed structural-scalar initial conditions. The reduced calculation now
provides a sharper expectation for that run: late perturbations should decay,
while the early radiation/matter branch may require the tracker handoff and
metric backreaction to carry the solution onto the mature structural
attractor.


## 16. Mature inverse-square coasting normal form

The shift-current discussion above identifies an important limitation.  A
strictly shift-symmetric constant-velocity root would generically drive
\(P_{,Z}\to0\), which is cosmological-constant-like rather than the required
mature SDMC coasting source.  The mature inverse-square branch instead suggests
a scale-covariant, explicitly \(\sigma\)-dependent normal form.

Take the asymptotic No-Slip/GR limit

\[
F\to F_\infty=\text{constant},
\qquad
F_{,\sigma}\to0,
\qquad
g\to0,
\]

and let the leading scalar Lagrangian be

\[
\boxed{
G_2(\sigma,Z)
=
\frac{1}{\sigma^2}
\left(
\kappa_1 Z+\kappa_2 Z^2-U_0
\right).
}
\]

Assume the synchronized constant-velocity orbit

\[
\sigma=\frac{t-t_B}{t_\star},
\qquad
v_\star\equiv\dot\sigma=\text{constant},
\qquad
H=\frac{v_\star}{\sigma},
\qquad
Z_\star=\frac12v_\star^2.
\]

The homogeneous scalar equation then gives

\[
\boxed{
U_0
=
2\kappa_1 Z_\star
+
3\kappa_2 Z_\star^2.
}
\]

The scalar energy density and pressure are

\[
\rho_\sigma
=
\frac{
\kappa_1 Z_\star
+3\kappa_2Z_\star^2
+U_0
}{\sigma^2},
\]

\[
p_\sigma
=
\frac{
\kappa_1 Z_\star
+\kappa_2Z_\star^2
-U_0
}{\sigma^2}.
\]

Using the constant-velocity field equation,

\[
\boxed{
\rho_\sigma
=
\frac{
3Z_\star(\kappa_1+2\kappa_2Z_\star)
}{\sigma^2},
}
\]

\[
\boxed{
p_\sigma
=
-\frac{
Z_\star(\kappa_1+2\kappa_2Z_\star)
}{\sigma^2}.
}
\]

Therefore

\[
\boxed{
w_\sigma=-\frac13.
}
\]

This is a major structural closure: the same inverse-square dependence that
defines the mature SDMC density produces the coasting equation of state
directly from the scalar field equation.  The relation \(w=-1/3\) is not
inserted separately.

If this scalar sector dominates a spatially flat future and \(F_\infty\) is
constant, the Friedmann equation gives

\[
\boxed{
\kappa_1+2\kappa_2Z_\star=2F_\infty.
}
\]

Consequently

\[
\rho_\sigma
=
\frac{6F_\infty Z_\star}{\sigma^2},
\qquad
p_\sigma
=
-\frac{2F_\infty Z_\star}{\sigma^2}.
\]

For the asymptotic k-essence normal form,

\[
P_{,Z}+2ZP_{,ZZ}
=
\frac{\kappa_1+6\kappa_2Z_\star}{\sigma^2},
\]

and the scalar sound speed is

\[
\boxed{
c_s^2
=
\frac{\kappa_1+2\kappa_2Z_\star}
{\kappa_1+6\kappa_2Z_\star}
=
\frac{F_\infty}
{F_\infty+2\kappa_2Z_\star}.
}
\]

Thus \(F_\infty>0\) and \(\kappa_2Z_\star>0\) give

\[
0<c_s^2<1.
\]

Unlike the pure shift-current \(P_{,Z}=0\) root, the inverse-square normal form
can therefore support a finite positive propagation speed while retaining a
constant structural velocity and exact coasting.

Defining

\[
r\equiv\frac{\kappa_2Z_\star}{F_\infty},
\]

one may write

\[
c_s^2=\frac1{1+2r},
\qquad
\kappa_1=2F_\infty(1-r),
\]

\[
U_0=F_\infty Z_\star(4-r).
\]

The region

\[
0<r<4
\]

has subluminal positive \(c_s^2\) and positive \(U_0\), while

\[
r>1
\]

also gives \(\kappa_1<0\).  The present reconstructed structural coefficient
\(k_1\) becomes negative near \(z\simeq17.41\), so the sign pattern required by
this candidate late normal form is already entered during the tracker-to-late
handoff, although the full scale-covariant asymptote has not yet been reached.

## 17. Present acceleration is not yet the mature normal form

The accepted \(z\ge0\) reconstruction can be tested for proximity to the
inverse-square normal form by examining

\[
\sigma^2k_1(\sigma),
\qquad
\sigma^2k_2(\sigma),
\qquad
\sigma^2V(\sigma).
\]

These quantities would tend to constants if the \(G_2\propto\sigma^{-2}\)
normal form had already been reached.

At the present endpoint the reconstructed logarithmic slopes are approximately

\[
\frac{d\ln|\sigma^2k_1|}{d\ln\sigma}
\simeq-0.640,
\]

\[
\frac{d\ln|\sigma^2k_2|}{d\ln\sigma}
\simeq+0.681,
\]

\[
\frac{d\ln|\sigma^2V|}{d\ln\sigma}
\simeq+1.931.
\]

Meanwhile

\[
\frac{d\ln F}{d\ln\sigma}
\simeq1.05\times10^{-3}.
\]

Thus the effective Planck mass is already close to frozen, but the reconstructed
\(G_2\) sector is not yet in its mature \(\sigma^{-2}\) form.  In particular,
the potential is still behaving much more like the nearly constant
activation-supported source required by the present accelerating phase.

This is consistent with the independently reconstructed activation rate,

\[
A_{\chi,0}\simeq1.99837,
\]

which is close to temporary cosmological-constant mimicry rather than the
future coasting requirement

\[
A_\chi\to0.
\]

The present accepted action therefore describes the transient accelerating
handoff.  The manuscript's mature coasting normal form is a future completion,
not a property that should already hold at \(z=0\).

## 18. Future lapse consistency and the role of the ideal structural coefficient

The density identity derived earlier is

\[
{\cal N}
=
p\,\Xi_P
\sqrt{\frac{\chi}{\Omega_X}},
\qquad
\Xi_P^2
=
\frac{8\pi G}{3}\rho_Pt_P^2.
\]

With the numerical Planck references used in the current structural audit,

\[
\boxed{
\Xi_P\simeq2.89445,
}
\]

which is essentially the manuscript's ideal flat-vacuum structural coefficient
\(\Xi_v=\sqrt{8\pi/3}\).

At the present accepted-action point,

\[
p_0=1.03421566,
\qquad
\chi_0=0.94108957,
\qquad
\Omega_{X0}=0.72374015,
\]

and therefore

\[
p_0\Xi_P
\sqrt{\frac{\chi_0}{\Omega_{X0}}}
=
3.41350846,
\]

exactly reproducing the chronological lapse.

Now consider the mature flat scalar-dominated normal form.  Then

\[
p\to1,
\]

and for constant \(F_\infty\),

\[
\Omega_X\to F_\infty.
\]

Hence

\[
\boxed{
{\cal N}_\infty
=
\Xi_P
\sqrt{\frac{\chi_\infty}{F_\infty}}.
}
\]

This creates a useful consistency fork.

If the activation factor is a literal bounded fraction and saturates at

\[
\chi_\infty=1,
\]

while ordinary gravity is recovered with

\[
F_\infty\simeq1,
\]

then

\[
\boxed{
{\cal N}_\infty\simeq2.89445,
}
\]

rather than the present chronological value \(3.41351\).

Keeping

\[
{\cal N}_\infty=3.41351
\]

with \(F_\infty=1\) instead requires

\[
\boxed{
\chi_\infty\simeq1.39082.
}
\]

If \(F_\infty\) remained at its present value
\(F_0\simeq1.02347\), the required value would be

\[
\chi_\infty\simeq1.42346.
\]

Therefore the three assumptions

\[
{\cal N}=\text{present constant forever},
\qquad
\chi_\infty\le1,
\qquad
F_\infty\simeq1
\]

cannot all hold simultaneously in the simplest flat scalar-dominated mature
limit.

This is not a contradiction in the current \(z\ge0\) reconstruction.  It says
that the future completion must choose among several physically distinct
possibilities:

1. the structural lapse relaxes from its present value toward the ideal
   \(\Xi_P\) value as activation saturates;
2. \(\chi\) is an effective coupling rather than a strictly bounded fraction
   and may exceed unity;
3. the asymptotic modified-gravity bookkeeping keeps
   \(\Omega_X\ne F_\infty\simeq1\);
4. the simple chronological relation \(S/S_0=t/t_0\) is only an excellent
   classical-history closure and is replaced by a slightly different future
   linear branch.

The first possibility has a particularly clean interpretation.  If
\(\chi_\infty=F_\infty=1\),

\[
\frac{{\cal N}_\infty}{{\cal N}_0}
\simeq0.84794.
\]

Thus the mature structural speed would be about \(15.2\%\) below the present
chronological slope.  This is compatible with the manuscript's earlier
distinction between the present structural coefficient near \(3.37\) and the
ideal flat-vacuum coefficient near \(2.894\).

The future asymptote therefore provides a new discriminator for the microscopic
Lagrangian: it must determine whether \({\cal N}\) remains fixed, relaxes toward
\(\Xi_P\), or is supported by a non-bounded effective activation coupling.
