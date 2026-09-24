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

For the structural-field action, two closely related density conventions must
be kept separate. The **positive structural source** used to define the SDMC
activation factor is the KGB-side density before the non-minimal Planck-mass
term is folded into an effective dark source,

\[
\boxed{
\rho_X^{\rm struct}=
\widehat k_1 Z+3\widehat k_2 Z^2+\widehat V
+6H\dot\sigma Z\widehat g
-2Z^2\widehat g_{,\sigma}
}.
\]

The background equation may equivalently be written

\[
3\widehat F H^2
=
\rho_m+\rho_r
+\rho_X^{\rm struct}
-3H\dot{\widehat F}.
\]

Hence the quantity that appears as a single effective right-hand-side source is

\[
\boxed{
\rho_X^{\rm eff}
=
\rho_X^{\rm struct}-3H\dot{\widehat F}.
}
\]

The positive source used in the Part-V/Part-IX structural-density audit is

\[
\boxed{
\frac{\rho_X^{\rm struct}}{3M_{\rm Pl}^2H_0^2}
=
E^2(F+F')
-\Omega_{m0}a^{-3}
-\Omega_{r0}a^{-4}.
}
\]

This distinction resolves an otherwise confusing notation collision in the
earlier continuation: the coefficient-level expression containing
\(-3H\dot F\) is the metric-RHS effective source, whereas the source whose
slope defines \(p_X\) and whose positivity is audited is
\(\rho_X^{\rm struct}\). The two coincide only when \(\dot F=0\).

In the remainder of the structural-clock discussion, \(\rho_X\) means
\(\rho_X^{\rm struct}\) unless an eff superscript is written explicitly.

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


## 19. Activation as the running of the rescaled action density

The mature normal form also suggests a cleaner interpretation of the activation
factor itself.  Define the rescaled active structural energy

\[
{\cal R}\equiv\sigma^2\rho_X.
\]

Because

\[
\rho_v=\frac{\rho_P}{S_0^2\sigma^2},
\qquad
\chi=\frac{\rho_X}{\rho_v},
\]

we have the exact identity

\[
\boxed{
\chi
=
\frac{S_0^2}{\rho_P}\,{\cal R}.
}
\]

Therefore the activation history is simply the running of the rescaled
action-level energy:

\[
\boxed{
A_\chi
=
\frac{d\ln{\cal R}}{d\ln a}.
}
\]

Since

\[
p=\frac{d\ln\sigma}{d\ln a},
\]

one may also define a structural beta function,

\[
\boxed{
\beta_\chi
\equiv
\frac{d\ln\chi}{d\ln\sigma}
=
\frac{A_\chi}{p}.
}
\]

The mature coasting fixed point is

\[
\beta_\chi\to0,
\]

because \({\cal R}\) becomes constant when
\(\rho_X\propto\sigma^{-2}\).

At the accepted present point,

\[
A_{\chi,0}=1.99836914,
\qquad
p_0=1.03421566,
\]

so

\[
\boxed{
\beta_{\chi,0}
\simeq1.93226.
}
\]

The reconstructed potential provides a striking cross-check.  At the same
endpoint,

\[
\frac{d\ln(\sigma^2V)}{d\ln\sigma}
\simeq1.93130,
\]

differing from \(\beta_{\chi,0}\) by only about

\[
9.5\times10^{-4}.
\]

Thus the present activation is almost entirely aligned with the running of the
rescaled potential term.  This is exactly what should occur during a
cosmological-constant-like phase: \(V\) is nearly constant while multiplication
by \(\sigma^2\) makes the rescaled structural energy grow with a slope near
two.

The future activation problem can therefore be reframed without introducing an
independent activation field:

\[
\boxed{
\text{late activation}
\quad\Longleftrightarrow\quad
\text{RG-like flow of the reconstructed structural functions toward }
G_2\propto\sigma^{-2}.
}
\]

In this picture the microscopic theory must explain why the present
near-constant-potential regime flows into the inverse-square normal form.  Once
that flow completes,

\[
A_\chi\to0,
\qquad
\chi\to\chi_\infty,
\qquad
w_X\to-\frac13.
\]

## 20. Exact structural-lapse evolution equation

The structural lapse admits two equivalent representations,

\[
{\cal N}=t_PS H p
\]

and

\[
{\cal N}
=
p\,\Xi_P
\sqrt{\frac{\chi}{\Omega_X}}.
\]

Taking a logarithmic derivative with respect to \(a\) gives

\[
\boxed{
A_{\cal N}
\equiv
\frac{d\ln{\cal N}}{d\ln a}
=
p-1-q
+
\frac{d\ln p}{d\ln a},
}
\]

where

\[
q=-1-\frac{\dot H}{H^2}.
\]

The density representation gives the equivalent relation

\[
\boxed{
A_{\cal N}
=
\frac{d\ln p}{d\ln a}
+
\frac12
\left(
A_\chi
-
\frac{d\ln\Omega_X}{d\ln a}
\right).
}
\]

These are identical because

\[
\frac{d\ln\Omega_X}{d\ln a}
=
A_\chi-2p+2(1+q).
\]

This equation clarifies the status of the chronological constant-lapse
closure.  If

\[
A_{\cal N}=0,
\]

then

\[
\boxed{
q
=
p-1
+
\frac{d\ln p}{d\ln a}.
}
\]

The present and past chronological branch approximately satisfies this
relation.  A future relaxation of \({\cal N}\), however, can occur without
breaking the structural framework; it corresponds to a controlled departure
from this special constant-lapse trajectory during the final activation
saturation.

For the bounded-activation future example

\[
\chi_\infty=1,
\qquad
F_\infty=1,
\qquad
{\cal N}_\infty=\Xi_P,
\]

the required integrated change is

\[
\ln\frac{{\cal N}_\infty}{{\cal N}_0}
\simeq-0.16495.
\]

Because \(p_\infty=1\),

\[
\boxed{
\int_{a_0}^{\infty}
\left(p-1-q\right)d\ln a
\simeq-0.13130.
}
\]

If \(p\) settles close to unity quickly, this becomes approximately

\[
\int_{a_0}^{\infty}q\,d\ln a
\simeq0.1313.
\]

Thus the bounded-activation route naturally prefers some net future
decelerating area before the universe settles to \(q\to0\).  This connects
directly with the manuscript's earlier observation that the final coasting
state may be approached from the weakly decelerating side rather than remaining
eternally accelerating.

The future asymptotic sign of \(q\) is therefore no longer merely qualitative:
once the microscopic flow of \(\chi\), \(p\), and \(F\) is specified, the lapse
evolution equation fixes the integrated approach to the mature coasting branch.


## 21. Reconciliation with the earlier 2.8944 asymptotic endpoint

The future-lapse result obtained above is not an isolated new number.  Part III
already derived, under its flat canonical scalar closure,

\[
{\cal N}
=
\Xi_v\,\frac{p}{\sqrt{\Omega_\phi}},
\qquad
\Xi_v=\sqrt{\frac{8\pi}{3}}\simeq2.894405,
\]

and therefore obtained the late scalar-dominated endpoint

\[
\boxed{
{\cal N}_\infty=\Xi_v\simeq2.8944.
}
\]

The current covariant density identity,

\[
{\cal N}
=
p\,\Xi_P
\sqrt{\frac{\chi}{\Omega_X}},
\]

is the active/bare generalization of that earlier relation.  In the canonical
limit \(\chi\to1\), \(F\to1\), and \(\Omega_X\to1\), the two formulas coincide.

Thus the present action-level development has recovered a result that was
already present in the earlier manuscript from a different route.  The current
value \({\cal N}_0=3.4135\) belongs to the matter-containing, activation-driven
present branch; the ideal \(2.8944\) value belongs to the asymptotic
scalar-dominated structural branch.

This distinction should not be confused with the unrelated numerical value

\[
\Delta\chi^2_{\rm fair,Planck+DESI}=+2.84892452
\]

that appeared at an earlier Part IX likelihood milestone before the
Planck+DESI residual was subsequently driven through zero.

## 22. General scale-covariant fixed-point theorem

The polynomial normal form of Section 16 is a special case of a more general
result.

Let the asymptotic scalar sector have

\[
\boxed{
G_2(\sigma,Z)=\sigma^{-m}f(Z),
}
\]

with constant \(F\), negligible \(G_3\), and a homogeneous constant-velocity
orbit

\[
\dot\sigma=v_\star=\text{constant},
\qquad
Z_\star=\frac12v_\star^2.
\]

Suppose simultaneously that

\[
a(t)\propto t^h.
\]

Since \(\sigma\propto t\),

\[
H=\frac{h}{t}=\frac{hv_\star}{\sigma}.
\]

The scalar equation

\[
\frac{d}{dt}(G_{2,Z}\dot\sigma)
+3HG_{2,Z}\dot\sigma
-G_{2,\sigma}=0
\]

reduces exactly to

\[
\boxed{
2Z_\star(3h-m)f_{,Z}
+
m f
=
0.
}
\]

Meanwhile,

\[
\rho_\sigma
=
\sigma^{-m}
\left(
2Z_\star f_{,Z}-f
\right).
\]

But the Friedmann side satisfies

\[
H^2\propto\sigma^{-2}.
\]

For a scalar-dominated constant-\(F\) solution with nonzero finite amplitude,
the powers can agree only if

\[
\boxed{
m=2.
}
\]

Thus the inverse-square structural dependence is not merely a convenient
choice for a constant-velocity power-law cosmology: it is the unique monomial
field scaling compatible with the flat Friedmann equation under those
assumptions.

With \(m=2\), the scalar equation becomes

\[
\boxed{
f+(3h-2)Z_\star f_{,Z}=0.
}
\]

The equation of state is then

\[
w_\sigma
=
\frac{f}{2Z_\star f_{,Z}-f}
=
-\frac{3h-2}{3h},
\]

hence

\[
\boxed{
w_\sigma=-1+\frac{2}{3h}.
}
\]

Because

\[
p
=
\frac{d\ln\sigma}{d\ln a}
=
\frac1h,
\]

we obtain

\[
\boxed{
w_\sigma=-1+\frac{2p}{3}.
}
\]

This is precisely the SDMC inverse-square structural equation-of-state relation
previously obtained from density scaling.  It now follows directly from the
scale-covariant Lagrangian fixed point.

The flat Friedmann amplitude gives a second exact relation.  Using

\[
3F H^2=\rho_\sigma
\]

together with the scalar equation yields

\[
\boxed{
f_{,Z}(Z_\star)=2Fh=\frac{2F}{p}.
}
\]

Therefore the structural exponent itself is encoded in the kinetic slope:

\[
\boxed{
p=\frac{2F}{f_{,Z}(Z_\star)}.
}
\]

For mature coasting,

\[
h=1,\qquad p=1,
\]

so the two fixed-point conditions reduce to

\[
\boxed{
f+Z_\star f_{,Z}=0,
\qquad
f_{,Z}=2F.
}
\]

For

\[
f(Z)=\kappa_1Z+\kappa_2Z^2-U_0,
\]

these reproduce exactly

\[
U_0=2\kappa_1Z_\star+3\kappa_2Z_\star^2,
\]

and

\[
\kappa_1+2\kappa_2Z_\star=2F.
\]

## 23. Running structural exponent and activation flow

The accepted action provides a direct variable that measures how far the
active source is from the inverse-square fixed point. Define

\[
\boxed{
m_X
\equiv
-\frac{d\ln\rho_X}{d\ln\sigma}.
}
\]

Since

\[
p_X
=
-\frac12\frac{d\ln\rho_X}{d\ln a},
\qquad
p=\frac{d\ln\sigma}{d\ln a},
\]

we have

\[
\boxed{
m_X=\frac{2p_X}{p}.
}
\]

The activation beta function satisfies

\[
\boxed{
\beta_\chi
\equiv
\frac{d\ln\chi}{d\ln\sigma}
=
\frac{A_\chi}{p}
=
2-m_X.
}
\]

Thus

\[
\boxed{
m_X+\beta_\chi=2
}
\]

is an exact consequence of \(\rho_X=\chi\rho_P/(S_0^2\sigma^2)\).

At the accepted present point,

\[
p_0=1.03421566,
\qquad
p_{X0}=0.03503109,
\]

so

\[
\boxed{
m_{X0}\simeq0.0677443,
}
\]

and

\[
\boxed{
\beta_{\chi0}\simeq1.932256.
}
\]

The active source is therefore still much closer to the
\(m_X=0\) cosmological-constant-like side than to the mature inverse-square
fixed point \(m_X=2\).

The effective equation of state can be written

\[
\boxed{
w_X
=
-1+\frac{2p_X}{3}
=
-1+\frac{p\,m_X}{3}.
}
\]

At the present point this gives

\[
w_{X0}\simeq-0.97665,
\]

while the mature limit

\[
p\to1,\qquad m_X\to2
\]

gives

\[
w_X\to-\frac13.
\]

A minimal autonomous fixed-point model for the future action flow is

\[
\boxed{
\frac{dm_X}{d\ln\sigma}
=
\gamma\,m_X(2-m_X),
\qquad
\gamma>0.
}
\]

It has

\[
m_X=0
\]

as an unstable de-Sitter-like fixed point and

\[
m_X=2
\]

as a stable inverse-square fixed point.  Equivalently,

\[
\frac{d\beta_\chi}{d\ln\sigma}
=
-\gamma\,\beta_\chi(2-\beta_\chi).
\]

The exact solution starting from \(m_{X0}\) at \(\sigma=1\) is

\[
\boxed{
m_X(\sigma)
=
\frac{2}
{1+C\sigma^{-2\gamma}},
\qquad
C=\frac{2-m_{X0}}{m_{X0}}.
}
\]

The corresponding activation growth is finite:

\[
\boxed{
\ln\frac{\chi_\infty}{\chi_0}
=
\frac{1}{\gamma}\ln(1+C).
}
\]

Therefore specifying the asymptotic activation normalization fixes the flow
rate,

\[
\boxed{
\gamma
=
\frac{\ln(1+C)}
{\ln(\chi_\infty/\chi_0)}.
}
\]

For a literal bounded activation with

\[
\chi_\infty=1,
\]

the present values imply approximately

\[
\gamma\simeq55.75.
\]

If instead the effective activation saturates near the current accepted
Planck-mass asymptote,

\[
F_\infty=e^{A_F}\simeq1.02389,
\qquad
\chi_\infty=F_\infty,
\]

then

\[
\gamma\simeq40.15.
\]

These large values quantify the statement already made in the manuscript:
because the present activation is close to saturation while
\(A_{\chi0}\simeq2\), the activation rate must fall rapidly in the future.

This autonomous equation is not yet claimed as the microscopic law.  It is a
compact normal-form candidate whose fixed points exactly match the two
action-level regimes now identified:

\[
\text{present nearly constant active density}
\quad
m_X\simeq0
\]

flowing toward

\[
\text{mature inverse-square structural density}
\quad
m_X=2.
\]

The next covariant task is to realize this \(m_X\) flow through smooth
future extensions of \(k_1(\sigma)\), \(k_2(\sigma)\), \(V(\sigma)\),
\(g(\sigma)\), and \(F(\sigma)\), and then release the extension under native
hi_class evolution rather than imposing \(m_X(\sigma)\) directly.


## 24. The 2.860--2.894 mature-lapse bracket

The earlier manuscript contains several nearby numbers that should remain
separate.

The structural one is

\[
\Xi_v=\sqrt{\frac{8\pi}{3}}=2.894405018\ldots,
\]

and Part III already identified

\[
{\cal N}_\infty=\Xi_v\simeq2.8944
\]

for the ideal late scalar-dominated canonical endpoint.

Separately, Part VI contains an unrelated \(2.86\%\) growth-suppression entry,
and Part IX contains an unrelated likelihood value
\(\Delta\chi^2\simeq2.86998\).  Neither of those is a structural lapse.

The current action-level completion nevertheless produces a genuine structural
interval close to the remembered \(2.86\)--\(2.89\) range.

For the accepted Planck-mass amplitude,

\[
F_\infty=e^{A_F},
\qquad
A_F=0.02360463334,
\]

so

\[
F_\infty=1.02388543.
\]

In the mature scalar-dominated branch,

\[
{\cal N}_\infty
=
\Xi_v\sqrt{\frac{\chi_\infty}{F_\infty}}.
\]

If the asymptotic activation lies in the natural interval

\[
1\le\chi_\infty\le F_\infty,
\]

then

\[
\boxed{
2.860445
\le
{\cal N}_\infty
\le
2.894405.
}
\]

The lower endpoint corresponds to

\[
\chi_\infty=1,
\qquad
F_\infty=1.02388543,
\]

while the upper endpoint corresponds to

\[
\chi_\infty=F_\infty
\]

or to the canonical limit

\[
\chi_\infty=F_\infty=1.
\]

Thus a structural range of approximately

\[
\boxed{
{\cal N}_\infty\simeq2.86\text{--}2.89
}
\]

arises naturally in the present covariant completion, even though the older
manuscript's exact structural reference was the upper endpoint \(2.8944\).

Relative to the accepted present lapse,

\[
{\cal N}_0=3.41350846,
\]

the corresponding mature structural-speed range is

\[
0.83798
\lesssim
\frac{{\cal N}_\infty}{{\cal N}_0}
\lesssim
0.84793.
\]

Since

\[
Z\propto{\cal N}^2,
\]

the mature kinetic density is therefore predicted to lie in the interval

\[
\boxed{
0.70221
\lesssim
\frac{Z_\infty}{Z_0}
\lesssim
0.71898.
}
\]

This sharpens the future completion problem: the present structural field need
not stop; its homogeneous speed relaxes by only about \(15\)--\(16\%\), while
its kinetic density relaxes by about \(28\)--\(30\%\), before settling on the
mature coasting fixed point.

## 25. Universal homogeneous attraction of the inverse-square fixed point

The scale-covariant fixed point admits a stronger stability statement than the
earlier frozen-background diagnostic.

Take the asymptotic scalar-dominated action

\[
G_2(\sigma,Z)=\sigma^{-2}f(Z),
\qquad
G_3\to0,
\qquad
G_4\to\frac{F}{2},
\]

with constant positive \(F\).  The flat Friedmann equation gives

\[
H\sigma
=
\sqrt{
\frac{2Zf_{,Z}-f}{3F}
}.
\]

Define

\[
{\cal K}(Z)
=
f_{,Z}+2Zf_{,ZZ}.
\]

Using \(\ln\sigma\) as the evolution variable, the exact homogeneous scalar
equation can be written

\[
\boxed{
\frac{dZ}{d\ln\sigma}
=
-\frac{{\cal B}(Z)}{{\cal K}(Z)},
}
\]

where

\[
\boxed{
{\cal B}(Z)
=
3(H\sigma)f_{,Z}\sqrt{2Z}
+
2f
-
4Zf_{,Z}.
}
\]

At the mature coasting fixed point,

\[
f+Z_\star f_{,Z}=0,
\]

and the Friedmann amplitude requires

\[
f_{,Z}(Z_\star)=2F.
\]

These conditions imply

\[
H\sigma=\sqrt{2Z_\star}.
\]

Now perturb

\[
Z=Z_\star+\delta Z.
\]

Because

\[
\frac{d}{dZ}
\left(2Zf_{,Z}-f\right)
=
{\cal K},
\]

one finds exactly at the fixed point

\[
{\cal B}_{,Z}\big|_\star
=
2{\cal K}_\star.
\]

Therefore

\[
\boxed{
\frac{d\,\delta Z}{d\ln\sigma}
=
-2\,\delta Z.
}
\]

Hence

\[
\boxed{
\delta Z\propto\sigma^{-2}.
}
\]

This decay exponent is independent of the detailed shape of \(f(Z)\).  The only
local requirements are that the fixed point exists and that

\[
{\cal K}_\star>0.
\]

Since

\[
{\cal N}
=
t_PS_0\sqrt{2Z},
\]

small lapse perturbations satisfy

\[
\boxed{
\frac{\delta{\cal N}}{{\cal N}}
\simeq
\frac12\frac{\delta Z}{Z_\star}
\propto\sigma^{-2}.
}
\]

Thus the mature structural lapse is not merely a consistent endpoint; within
the scalar-dominated inverse-square normal form it is a genuine homogeneous
attractor.

For the quadratic family

\[
f(Z)
=
\kappa_1Z+\kappa_2Z^2-U_0,
\]

define

\[
r=\frac{\kappa_2Z_\star}{F}.
\]

Then

\[
c_{s,\star}^2=\frac1{1+2r}.
\]

If, only as a continuity reference, the asymptotic sound speed is chosen equal
to the accepted present scalar value

\[
c_{s,0}^2\simeq0.138335,
\]

then

\[
\boxed{
r\simeq3.11441.
}
\]

For this value the positive-gradient region is

\[
\frac{Z}{Z_\star}
>
1-\frac1r
\simeq0.67891,
\]

while the homogeneous kinetic coefficient remains positive already for

\[
\frac{Z}{Z_\star}
>
\frac{r-1}{3r}
\simeq0.22630.
\]

A nonlinear basin integration of the exact scalar-dominated equations confirms
convergence toward \(Z/Z_\star=1\) from representative initial values between
\(0.70\) and \(2.0\).  For example, after
\(\Delta\ln\sigma=3\), the \(r\simeq3.1144\) cases give approximately

\[
Z/Z_\star=0.99913
\]

from an initial value \(0.70\), and

\[
Z/Z_\star=1.00169
\]

from an initial value \(2.0\).

The nonlinear result is consistent with the universal local law
\(\delta Z\propto\sigma^{-2}\).

This is stronger than the earlier late-time frozen-background stability test:
the metric response is included through the scalar-dominated Friedmann
constraint.  What remains to be demonstrated is that a smooth continuation of
the accepted late266 Horndeski functions enters this normal-form basin without
crossing \(D=0\), \(c_s^2=0\), or generating an unacceptable future singularity.

## 26. Smooth covariant stitching conditions

A future completion should leave the accepted \(z\ge0\) action unchanged and
modify only the unobserved \(\sigma>1\) branch.

A convenient variable is

\[
x=\ln\sigma.
\]

For any coefficient \(Q\) that must approach an inverse-square mature form,

\[
Q(\sigma)\longrightarrow Q_\infty\sigma^{-2},
\]

define the rescaled quantity

\[
\bar Q(x)=e^{2x}Q(x).
\]

A \(C^2\) future tail can be written

\[
\boxed{
\bar Q(x)
=
Q_\infty
+
e^{-\mu x}
\left(
A_0+A_1x+\frac12A_2x^2
\right),
\qquad x\ge0.
}
\]

If the accepted action supplies the present value and first two derivatives

\[
\bar Q_0,\qquad
\bar Q'_0,\qquad
\bar Q''_0,
\]

then exact \(C^2\) matching at \(x=0\) requires

\[
\boxed{
A_0=\bar Q_0-Q_\infty,
}
\]

\[
\boxed{
A_1=\bar Q'_0+\mu A_0,
}
\]

\[
\boxed{
A_2=\bar Q''_0+2\mu\bar Q'_0+\mu^2A_0.
}
\]

This construction can be applied separately to the future
\(k_1,k_2,V\) coefficients.

For the Planck mass, which approaches a constant rather than an inverse-square
law, the same formula is used directly on \(F(x)\) with \(F_\infty\) replacing
\(Q_\infty\).

The No-Slip relation should not be independently interpolated.  Once a smooth
future \(F(\sigma)\) and a target kinetic branch \(Z(\sigma)\) are chosen,
define

\[
\boxed{
g(\sigma)
=
-\frac{F_{,\sigma}}{2Z(\sigma)}.
}
\]

Then

\[
ZG_{3,Z}
=
Zg
=
-\frac12F_{,\sigma}
\]

remains exact by construction and \(g\to0\) automatically as
\(F_{,\sigma}\to0\).

Likewise, a \(C^2\) structural-speed relaxation from the present chronological
kinetic density \(Z_0\) to the mature value \(Z_\infty\) can be represented by

\[
\boxed{
Z(x)
=
Z_\infty
+
(Z_0-Z_\infty)
e^{-\mu_Zx}
\left(
1+\mu_Zx+\frac12\mu_Z^2x^2
\right).
}
\]

This satisfies

\[
Z(0)=Z_0,
\qquad
Z'(0)=0,
\qquad
Z''(0)=0,
\]

so it joins smoothly onto the presently constant-\(Z\) chronological branch,
while

\[
Z\to Z_\infty
\]

in the future.

The next numerical task is therefore well defined: construct these future
coefficient tails from the actual accepted endpoint derivatives, release the
extended action under native evolution, and test whether the solution reaches
the universal inverse-square basin without imposing the target trajectory.


## 27. Exact homogeneous covariant evolution without imposing the clock

The future stitch can be released under the homogeneous equations of the same
linear-\(G_3\) Horndeski action rather than evolving a prescribed
\(Z(\sigma)\).

Write

\[
v\equiv\dot\sigma,
\qquad
Z=\frac12v^2,
\qquad
B\equiv\frac{k_2}{4}-\frac{g_{,\sigma}}{6}.
\]

The action-level energy density is

\[
\rho_X
=
\frac12k_1v^2
+
\left(
\frac34k_2-\frac12g_{,\sigma}
\right)v^4
+
V
+
3Hgv^3
-
3HF_{,\sigma}v.
\]

Hence the flat Friedmann constraint is a quadratic equation for \(H\),

\[
\boxed{
3FH^2
-
3\left(gv^3-F_{,\sigma}v\right)H
-
{\cal E}
=
0,
}
\]

where

\[
{\cal E}
=
\rho_m+\rho_r
+
\frac12k_1v^2
+
\left(
\frac34k_2-\frac12g_{,\sigma}
\right)v^4
+
V.
\]

The expanding solution is

\[
\boxed{
H
=
\frac{
3A_H+\sqrt{9A_H^2+12F{\cal E}}
}{6F},
\qquad
A_H\equiv gv^3-F_{,\sigma}v.
}
\]

Thus \(H\) is determined algebraically once
\(\sigma,v,\rho_m,\rho_r\) are known.

The scalar equation and Raychaudhuri equation are both linear in
\(\dot v\) and \(\dot H\).  Define

\[
{\cal K}
=
k_1+12Bv^2+6Hgv,
\]

\[
Q
=
F_{,\sigma}-gv^2,
\]

\[
{\cal R}
=
\rho_m+\frac43\rho_r
+k_1v^2
+\left(k_2-g_{,\sigma}\right)v^4
+3Hgv^3
-HF_{,\sigma}v
+F_{,\sigma\sigma}v^2,
\]

and

\[
\begin{aligned}
{\cal C}
={}&
V_{,\sigma}
+\frac12k_{1,\sigma}v^2
+3B_{,\sigma}v^4
+3Hk_1v
+12HBv^3
+2Hg_{,\sigma}v^3\\
&+9H^2gv^2
-6H^2F_{,\sigma}.
\end{aligned}
\]

The exact homogeneous system becomes

\[
\boxed{
\begin{pmatrix}
{\cal K} & -3Q\\
Q & 2F
\end{pmatrix}
\begin{pmatrix}
\dot v\\
\dot H
\end{pmatrix}
=
-
\begin{pmatrix}
{\cal C}\\
{\cal R}
\end{pmatrix}.
}
\]

Its determinant is

\[
\boxed{
\Delta_{\rm hom}
=
2F{\cal K}+3Q^2.
}
\]

Therefore

\[
\boxed{
\dot v
=
-\frac{
2F{\cal C}+3Q{\cal R}
}{
2F{\cal K}+3Q^2
},
}
\]

and

\[
\boxed{
\dot H
=
\frac{
Q{\cal C}-{\cal K}{\cal R}
}{
2F{\cal K}+3Q^2
}.
}
\]

Using the cosmological e-fold \(N_a=\ln a\),

\[
\boxed{
\frac{d\sigma}{dN_a}
=
\frac{v}{H},
}
\]

\[
\boxed{
\frac{dv}{dN_a}
=
\frac{\dot v}{H}.
}
\]

This closes the future homogeneous system with ordinary matter and radiation
dilution,

\[
\rho_m\propto a^{-3},
\qquad
\rho_r\propto a^{-4}.
\]

No chronological relation such as
\(S/S_0=t/t_0\) is imposed in this integration.  The structural lapse is
read out afterwards from

\[
{\cal N}
=
t_PS_0v.
\]

The same is true for

\[
p
=
\frac{v}{H\sigma}.
\]

This is the first stage of the present investigation in which the future
structural clock is evolved directly from the stitched covariant action rather
than inserted as a closure condition.

## 28. Present-endpoint consistency of the released equations

Before using the future extension, the coupled system was evaluated at the
accepted present endpoint using the reconstructed action itself.

The algebraic Friedmann solution reproduces the accepted free background with

\[
\left|
\frac{H_{\rm EOM}}{H_{\rm free}}-1
\right|
\simeq
1.8\times10^{-8}.
\]

The released system gives

\[
p_{\rm EOM}
=
1.03421568,
\]

against

\[
p_{\rm free}
=
1.03421566,
\]

and

\[
q_{\rm EOM}
\simeq
-0.53565562,
\]

against the direct background value

\[
q_{\rm free}
\simeq
-0.53565646.
\]

Thus the structural-field homogeneous equations reproduce the present accepted
solution before any future asymptotic behavior is tested.

For the linear-\(G_3\) action, the physical No-Slip combination is

\[
\boxed{
\alpha_B+2\alpha_M
=
\frac{
2v
}{
HF
}
\left(
Zg+\frac12F_{,\sigma}
\right).
}
\]

This is a better future diagnostic than dividing the two No-Slip terms by
their own instantaneous magnitudes, because both terms become very small as
\(F_{,\sigma}\to0\) and \(g\to0\).

## 29. First freely evolved future coasting candidates

Two \(C^2\) future tails were released under the exact homogeneous system.

### 29.1 Bounded-activation endpoint

For

\[
\chi_\infty=1,
\]

the structural target is

\[
\boxed{
{\cal N}_\infty
=
2.86044513,
}
\]

with

\[
\frac{Z_\infty}{Z_0}
=
0.70220719.
\]

A representative stable stitch uses

\[
r=3.43684,
\qquad
\mu_Q=3.75024,
\qquad
\mu_F=4.47766.
\]

After ten future e-folds the released solution gives

\[
\boxed{
{\cal N}=2.86044025,
}
\]

\[
\boxed{
p=0.99999826,
}
\]

and

\[
\boxed{
q=-1.68\times10^{-6}.
}
\]

Thus the freely evolved field reaches the intended mature structural speed and
coasting kinematics without directly imposing \(v=v_\infty\).

The maximum physical No-Slip departure in this future run is only

\[
\boxed{
\max|\alpha_B+2\alpha_M|
\simeq
5.7\times10^{-5}.
}
\]

The homogeneous kinetic coefficient remains positive in the tested interval.
The k-essence-sector sound-speed proxy remains in the approximate range

\[
0.127
\lesssim
c_{s,\rm proxy}^2
\lesssim
0.215.
\]

This is not yet the full Horndeski \(c_s^2\) test.

### 29.2 \(F_\infty\)-matched activation endpoint

For

\[
\chi_\infty=F_\infty,
\]

the ideal manuscript endpoint is recovered,

\[
\boxed{
{\cal N}_\infty
=
2.89440502.
}
\]

Here

\[
\frac{Z_\infty}{Z_0}
=
0.71897971.
\]

A representative released stitch uses

\[
r=3.75612,
\qquad
\mu_Q=3.75390,
\qquad
\mu_F=4.49987.
\]

After ten future e-folds,

\[
\boxed{
{\cal N}=2.89440058,
}
\]

\[
\boxed{
p=0.99999844,
}
\]

and

\[
\boxed{
q=-1.50\times10^{-6}.
}
\]

The maximum physical No-Slip departure is

\[
\boxed{
\max|\alpha_B+2\alpha_M|
\simeq
5.3\times10^{-5}.
}
\]

The k-essence-sector sound-speed proxy remains between approximately

\[
0.117
\lesssim
c_{s,\rm proxy}^2
\lesssim
0.207.
\]

Again, this proxy is not a substitute for the complete Horndeski scalar-health
calculation.

Both future candidates therefore exhibit the same qualitative release:

\[
{\cal N}_0\simeq3.4135
\]

relaxes dynamically toward the mature interval,

\[
{\cal N}_\infty\simeq2.86\text{--}2.89,
\]

while

\[
p\to1,
\qquad
q\to0,
\qquad
F\to F_\infty,
\qquad
g\to0.
\]

The current acceleration ends naturally in these candidate continuations.
For the \(F_\infty\)-matched example, the first
\(q=0\) crossing occurs near

\[
\ln a\simeq1.155,
\qquad
a\simeq3.17,
\]

about \(23.5\) Gyr after the present endpoint in this particular tail design.
The solution then enters a weakly decelerating interval, with a maximum

\[
q\simeq0.088
\]

around \(58\) Gyr after the present, before tending back toward

\[
q\to0.
\]

These future times are candidate-tail diagnostics rather than observationally
fixed predictions.

The most important result is narrower: the accepted action can be extended
smoothly into a future inverse-square branch whose freely evolved homogeneous
solution selects the same \(2.86\)--\(2.89\) structural interval identified
independently from the density closure and the old
\(\Xi_v\simeq2.8944\) endpoint.

The remaining promotion gate is now sharper.  The future extension must be
implemented in a genuinely future-capable hi_class background so that the full
Horndeski

\[
D>0,
\qquad
c_s^2>0,
\qquad
c_s^2\le1
\]

conditions can be checked along the complete released trajectory rather than
through the reduced homogeneous and k-essence diagnostics alone.


## 30. Full background-level Horndeski stability of the released future tails

The reduced k-essence proxy in Section 29 can now be replaced by the exact
Bellini-Sawicki background stability functions for the accepted Horndeski
subclass

\[
G_2=k_1(\sigma)Z+k_2(\sigma)Z^2-V(\sigma),
\qquad
G_3=g(\sigma)Z,
\qquad
G_4=\frac12F(\sigma),
\]

with \(G_5=0\) and \(c_T^2=1\).

Along the released homogeneous trajectory,

\[
\boxed{
\alpha_M
=
\frac{vF_{,\sigma}}{HF},
}
\]

\[
\boxed{
\alpha_B
=
\frac{2v}{HF}
\left(
Zg-\frac12F_{,\sigma}
\right),
}
\]

and

\[
\boxed{
\alpha_K
=
\frac{
2Z\left(k_1+6k_2Z-4Zg_{,\sigma}\right)
+
12HvZg
}{
H^2F
}.
}
\]

The exact scalar kinetic combination is therefore

\[
\boxed{
D
=
\alpha_K
+
\frac32\alpha_B^2.
}
\]

For this \(c_T=1\) Horndeski subclass, the scalar sound speed may be written

\[
\boxed{
c_s^2
=
\frac{
(2-\alpha_B)
\left[
-\frac{\dot H}{H^2}
+
\frac12\alpha_B
+
\alpha_M
\right]
-
\frac{\rho_m+\frac43\rho_r}{H^2F}
+
\frac{d\alpha_B}{d\ln a}
}{
D
}.
}
\]

This expression was cross-checked at the accepted present endpoint against the
native hi_class background columns.  The reconstructed present value is

\[
D_0\simeq0.38666113,
\]

against the direct hi_class value

\[
D_0^{\rm hi\_class}=0.38666113,
\]

while

\[
c_{s,0}^2\simeq0.13836
\]

is consistent with the direct native value

\[
c_{s,0}^{2\,{\rm hi\_class}}
=
0.13833497
\]

at the few-\(10^{-5}\) level.  The residual difference is dominated by the
endpoint derivative used for \(d\alpha_B/d\ln a\) in the future continuation.

### 30.1 Bounded-activation future tail

For the released

\[
\chi_\infty=1
\]

candidate, the complete background-level Horndeski audit over

\[
0\le\ln a\le10
\]

gives

\[
\boxed{
D_{\min}
=
0.38666113>0,
}
\]

with the minimum occurring at the present endpoint.

The sound speed remains positive and subluminal throughout the tested future:

\[
\boxed{
0.126971
\lesssim
c_s^2
\lesssim
0.214469.
}
\]

At ten future e-folds,

\[
D\simeq15.74725,
\qquad
c_s^2\simeq0.12700435.
\]

### 30.2 \(F_\infty\)-matched activation tail

For

\[
\chi_\infty=F_\infty,
\]

the full background-level stability audit gives

\[
\boxed{
D_{\min}
=
0.38666113>0,
}
\]

and

\[
\boxed{
0.117446
\lesssim
c_s^2
\lesssim
0.206273.
}
\]

At ten future e-folds,

\[
D\simeq17.02435,
\qquad
c_s^2\simeq0.11747705.
\]

Thus neither candidate encounters a ghost crossing or a scalar gradient
instability in the entire released homogeneous interval.

The physical No-Slip combination remains extremely small:

\[
\max|\alpha_B+2\alpha_M|
\simeq5.7\times10^{-5}
\]

for the bounded-\(\chi\) tail and

\[
\max|\alpha_B+2\alpha_M|
\simeq5.3\times10^{-5}
\]

for the \(F_\infty\)-matched tail.

These values clarify the previously quoted few-percent term-balance residual:
that ratio becomes misleading when both \(Zg\) and
\(F_{,\sigma}/2\) have already become individually tiny.  The physical
dimensionless No-Slip departure remains at the \(10^{-5}\) level and decays
rapidly toward zero.

### 30.3 Exact mature stability invariant

At the mature fixed point,

\[
\alpha_M\to0,
\qquad
\alpha_B\to0,
\qquad
g\to0,
\qquad
F\to F_\infty.
\]

For the quadratic inverse-square normal form, with

\[
r\equiv\frac{\kappa_2Z_\infty}{F_\infty},
\]

the fixed-point relations give

\[
\boxed{
D_\infty
=
2(1+2r),
}
\]

whereas

\[
\boxed{
c_{s,\infty}^2
=
\frac1{1+2r}.
}
\]

Hence the mature branch obeys the exact invariant

\[
\boxed{
D_\infty c_{s,\infty}^2=2.
}
\]

For the bounded-activation candidate,

\[
r=3.43684444,
\]

so

\[
D_\infty=15.74737778,
\qquad
c_{s,\infty}^2=0.12700527.
\]

For the \(F_\infty\)-matched candidate,

\[
r=3.75611884,
\]

so

\[
D_\infty=17.02447535,
\qquad
c_{s,\infty}^2=0.11747792.
\]

The numerical future integrations approach these analytic limits.

This closes the background-level ghost and scalar-gradient gate for the two
released future tails.  It does not yet replace an actual future perturbation
integration inside hi_class: the stock background solver terminates at
\(\ln(a/a_0)=0\), so a separate future-capable numerical branch is still needed
before the perturbation module itself can be propagated beyond the present
epoch.


## 31. C2-preserving iterative No-Slip refinement

The first freely evolved future tails already gave a small physical No-Slip
departure,

\[
|\alpha_B+2\alpha_M|\lesssim6\times10^{-5},
\]

but the residual can be reduced further without changing the accepted present
action or its first two field derivatives.

After releasing a candidate future tail, regard the resulting monotonic
trajectory as

\[
Z=Z_{\rm rel}(\sigma).
\]

The exact on-trajectory No-Slip value of the cubic coefficient is

\[
\boxed{
g_{\rm NS}(\sigma)
=
-\frac{F_{,\sigma}}
{2Z_{\rm rel}(\sigma)}.
}
\]

Directly replacing \(g\) by this expression would slightly change the future
side of the \(C^2\) jet at \(\sigma=1\). Instead define

\[
x=\ln\sigma,
\]

and the blending function

\[
\boxed{
W(x)
=
1-e^{-u}
\left(
1+u+\frac12u^2
\right),
\qquad
u=\lambda_{\rm NS}x.
}
\]

It satisfies

\[
W(0)=W'(0)=W''(0)=0,
\]

while

\[
W(x)\to1
\]

rapidly for \(x>0\).

The refined cubic coefficient is therefore

\[
\boxed{
g_{\rm ref}
=
g_{\rm base}
+
W(x)
\left(
g_{\rm NS}-g_{\rm base}
\right).
}
\]

Consequently,

\[
g_{\rm ref}(1)=g_{\rm base}(1),
\]

\[
g_{{\rm ref},\sigma}(1)
=
g_{{\rm base},\sigma}(1),
\]

and

\[
g_{{\rm ref},\sigma\sigma}(1)
=
g_{{\rm base},\sigma\sigma}(1).
\]

Thus the present accepted action and the complete \(C^2\) matching data are
unchanged.

A fixed-point iteration was then performed:

\[
g^{(n)}
\longrightarrow
Z^{(n)}_{\rm rel}(\sigma)
\longrightarrow
g^{(n+1)}_{\rm ref}.
\]

Using two refinement iterations and

\[
\lambda_{\rm NS}=40,
\]

the bounded-activation branch improves from

\[
\max|\alpha_B+2\alpha_M|
\simeq5.7\times10^{-5}
\]

to approximately

\[
\boxed{
9.3\times10^{-6},
}
\]

while the \(F_\infty\)-matched branch improves from

\[
5.3\times10^{-5}
\]

to approximately

\[
\boxed{
8.8\times10^{-6}.
}
\]

In the independent refinement audit the numerical changes of the endpoint
\(g\), \(g_{,\sigma}\), and \(g_{,\sigma\sigma}\) were zero at the quoted
precision, as required by the construction.

The full background-level Horndeski health conditions remain intact. For the
refined bounded-\(\chi\) trajectory,

\[
D_{\min}\simeq0.38666113,
\]

and

\[
0.126971
\lesssim
c_s^2
\lesssim
0.21449.
\]

For the refined \(F_\infty\)-matched trajectory,

\[
D_{\min}\simeq0.38666113,
\]

and

\[
0.117446
\lesssim
c_s^2
\lesssim
0.20630.
\]

The present deceleration parameter is unchanged by the \(C^2\) refinement,

\[
q_0\simeq-0.53565562,
\]

and the mature structural endpoint is unchanged to the numerical accuracy of
the test.

This demonstrates that the small transient No-Slip error is not tied to the
future coasting dynamics. It is mainly a choice of interpolation for
\(g(\sigma)\), and it can be systematically driven downward while preserving
the accepted endpoint jet.

The next numerical promotion criterion is therefore stronger than approximate
No-Slip. A future production tail can demand

\[
|\alpha_B+2\alpha_M|<10^{-5}
\]

throughout the released background, together with

\[
D>0,
\qquad
0<c_s^2<1,
\]

before the tail is admitted to a future-capable perturbation run.


## 32. Future scalar perturbations of the mature coasting fixed point

The background-level conditions

\[
D>0,
\qquad
c_s^2>0
\]

exclude ghost and gradient instabilities, but one can go further and ask
whether the mature inverse-square branch admits any growing scalar curvature
mode.

For the Bellini-Sawicki curvature variable \(\zeta\), the quadratic action may
be written schematically as

\[
S_\zeta^{(2)}
=
\int dt\,d^3x\,
a^3 Q_s
\left[
\dot\zeta^2
-
c_s^2\frac{(\nabla\zeta)^2}{a^2}
\right],
\]

with

\[
\boxed{
Q_s
=
\frac{2FD}{(2-\alpha_B)^2}.
}
\]

The Fourier-mode equation is therefore

\[
\boxed{
\ddot\zeta_k
+
\left(
3H+\frac{\dot Q_s}{Q_s}
\right)
\dot\zeta_k
+
c_s^2\frac{k^2}{a^2}\zeta_k
=
0.
}
\]

At the mature inverse-square fixed point,

\[
\alpha_M\to0,
\qquad
\alpha_B\to0,
\qquad
F\to F_\infty,
\]

and Section 30 gave

\[
D_\infty=2(1+2r),
\qquad
c_{s,\infty}^2=\frac1{1+2r}.
\]

Hence

\[
\boxed{
Q_{s,\infty}
=
F_\infty(1+2r)
}
\]

is constant.

Because mature coasting has

\[
a\propto t,
\qquad
H=\frac1t,
\qquad
aH=\text{constant},
\]

the dimensionless physical wavenumber

\[
\nu^2
\equiv
c_s^2\frac{k^2}{a^2H^2}
\]

is itself constant.  Writing

\[
\zeta_k\propto a^s
\]

gives

\[
\boxed{
s^2+2s+\nu^2=0,
}
\]

so

\[
\boxed{
s_\pm
=
-1\pm\sqrt{1-\nu^2}.
}
\]

This immediately separates the mature perturbations into two regimes.

For

\[
\nu\ll1,
\]

the roots approach

\[
s_+\to0,
\qquad
s_-\to-2,
\]

so the super-horizon solution is a constant curvature mode plus a decaying

\[
a^{-2}
\]

mode.

For

\[
\nu>1,
\]

the roots are

\[
s_\pm
=
-1
\pm
i\sqrt{\nu^2-1},
\]

so sub-horizon curvature perturbations oscillate with envelope

\[
\boxed{
|\zeta_k|\propto a^{-1}.
}
\]

There is therefore no growing scalar curvature mode at the mature fixed point.

Tensor modes satisfy the same qualitative coasting equation because

\[
F\to\text{constant},
\qquad
c_T^2=1.
\]

Their long-wavelength branch consists of a constant mode plus a decaying mode,
while finite-wavelength modes oscillate with a decaying envelope.

Ordinary matter also becomes asymptotically subdominant,

\[
\frac{\rho_m}{\rho_X}
\propto
a^{-1}
\to0.
\]

Its leading growth equation therefore tends to

\[
\ddot\delta_m+2H\dot\delta_m\simeq0,
\]

with asymptotic solutions

\[
\boxed{
\delta_m
=
C_1+\frac{C_2}{t}.
}
\]

Thus the mature branch does not generate a divergent future matter-growth
mode either.

## 33. Perturbation envelope through the finite future transition

The mature analytic result does not by itself exclude transient amplification
during the finite transition from the present accelerating branch to the
coasting fixed point.  The released backgrounds were therefore used to evolve
the curvature-mode equation directly over

\[
0\le\ln a\le6.
\]

Using e-fold time,

\[
\boxed{
\zeta_{,NN}
+
\left[
3+\frac{d\ln H}{dN}
+\frac{d\ln Q_s}{dN}
\right]
\zeta_{,N}
+
c_s^2
\left(
\frac{k}{aH}
\right)^2
\zeta
=
0.
}
\]

For the bounded-activation tail, the exact background gives

\[
Q_s>0
\]

throughout the tested interval, with

\[
Q_{s,\min}\simeq0.19744,
\]

and the total damping coefficient obeys approximately

\[
\boxed{
1.9994
\lesssim
3+\frac{d\ln H}{dN}
+\frac{d\ln Q_s}{dN}
\lesssim
11.63.
}
\]

For the \(F_\infty\)-matched tail,

\[
Q_{s,\min}\simeq0.19744,
\]

and

\[
\boxed{
2.0000
\lesssim
3+\frac{d\ln H}{dN}
+\frac{d\ln Q_s}{dN}
\lesssim
11.80.
}
\]

The friction term is therefore positive over the full released transition in
both candidates.

Representative modes were initialized at the present endpoint with

\[
\zeta(0)=1,
\qquad
\frac{d\zeta}{dN}(0)=0.
\]

Using the present-horizon ratios

\[
\nu_0\equiv\frac{k}{H_0}
=
0,\ 0.1,\ 1,\ 10,
\]

no tested mode exceeded its initial amplitude.

For the bounded-\(\chi\) future tail, the amplitudes at
\(\ln a=6\) are approximately

\[
\zeta=
1,\quad
0.99776,\quad
0.79573,\quad
-1.14\times10^{-3},
\]

for the four representative \(\nu_0\) values respectively.

For the \(F_\infty\)-matched tail, the corresponding values are

\[
\zeta=
1,\quad
0.99796,\quad
0.81278,\quad
-7.02\times10^{-4}.
\]

The maximum absolute amplitude in every tested case is

\[
\boxed{
\max|\zeta|=1.
}
\]

Thus the finite future transition shows no scalar-curvature transient
amplification in this representative mode set, and the numerical solutions
join smoothly onto the analytic mature behavior: a conserved long-wavelength
mode and decaying finite-wavelength modes.

Scientific status.  This closes a stronger future stability gate than the
background \(D\) and \(c_s^2\) conditions alone.  It still does not replace a
complete future multi-species Boltzmann integration, because the stock
hi_class background table terminates at \(a=1\).  What has now been shown is
that the released covariant background is ghost-free, gradient-stable, has a
positive scalar kinetic normalization, and does not exhibit a growing
curvature mode either during the tested transition or at the exact mature
fixed point.


## 34. Full Bellini-Sawicki health of the released future branch

The reduced homogeneous kinetic proxy can now be replaced by the exact
Bellini-Sawicki functions for the accepted Horndeski subclass

\[
G_2=k_1(\sigma)Z+k_2(\sigma)Z^2-V(\sigma),
\qquad
G_3=g(\sigma)Z,
\qquad
G_4=\frac{F(\sigma)}{2}.
\]

With

\[
v=\dot\sigma,
\qquad
Z=\frac12v^2,
\]

the three nonzero alpha functions are

\[
\boxed{
\alpha_M
=
\frac{vF_{,\sigma}}{HF},
}
\]

\[
\boxed{
\alpha_B
=
\frac{v\left(-F_{,\sigma}+2Zg\right)}{HF},
}
\]

and

\[
\boxed{
\alpha_K
=
\frac{2Z}{H^2F}
\left(
k_1+6k_2Z-4Zg_{,\sigma}+6Hgv
\right).
}
\]

Hence

\[
\boxed{
D
=
\alpha_K+\frac32\alpha_B^2.
}
\]

The structural No-Slip relation

\[
Zg=-\frac12F_{,\sigma}
\]

immediately yields

\[
\boxed{
\alpha_B=-2\alpha_M.
}
\]

The future trajectory can therefore be tested with exactly the same scalar
health quantities used by hi_class, without replacing them by a k-essence
proxy.

For this subclass the effective hi_class density and pressure entering the
Bellini-Sawicki sound-speed numerator are

\[
\boxed{
\rho_{\rm smg}
=
\frac{
k_1Z+3k_2Z^2+V-2g_{,\sigma}Z^2
}{3}
-
Hv\left(F_{,\sigma}-2Zg\right)
-
(F-1)H^2,
}
\]

and

\[
\boxed{
\begin{aligned}
p_{\rm smg}
={1\over3}\Big[
&
k_1Z+k_2Z^2-V
-2g_{,\sigma}Z^2
+2F_{,\sigma\sigma}Z\\
&
+3(F-1)H^2
+2(F-1)\dot H
+2F_{,\sigma}Hv\\
&
+\left(F_{,\sigma}-2Zg\right)\dot v
\Big].
\end{aligned}
}
\]

These are bookkeeping variables internal to the hi_class background equations
and should not be confused with the positive action-level source
\(\rho_X\) used in the structural density interpretation.

For Horndeski with

\[
\alpha_T=\alpha_H=0,
\]

the hi_class scalar sound-speed numerator reduces to

\[
\begin{aligned}
c_{s,\rm num}^2
={}&
\frac12(2-\alpha_B)(\alpha_B+2\alpha_M)\\
&+\frac32(2-\alpha_B)
\frac{\rho_{\rm smg}+p_{\rm smg}}{H^2}\\
&-\frac32
\frac{
2-2F+\alpha_BF
}{F}
\frac{\rho_{\rm m,CLASS}+p_{\rm m,CLASS}}{H^2}\\
&+
\frac{d\alpha_B}{d\ln a},
\end{aligned}
\]

with

\[
\boxed{
c_s^2=\frac{c_{s,\rm num}^2}{D}.
}
\]

An algebraically equivalent background form is

\[
\boxed{
c_s^2
=
\frac{
(2-\alpha_B)
\left(
-\frac{d\ln H}{d\ln a}
+\frac{\alpha_B}{2}
+\alpha_M
\right)
-
\frac{
\rho_m+\frac43\rho_r
}{H^2F}
+
\frac{d\alpha_B}{d\ln a}
}{D},
}
\]

where the last expression uses the action normalization for
\(\rho_m,\rho_r\).

The two independent implementations agree along the released future
trajectory to approximately

\[
\boxed{
\max|\Delta c_s^2|
\lesssim2\times10^{-5},
}
\]

the remaining difference being dominated by independent spline derivatives.

At the present endpoint, the structural reconstruction gives

\[
\boxed{
D_0\simeq0.38666113,
}
\]

in agreement with the stored hi_class value at the few-\(10^{-9}\) level, and

\[
\boxed{
c_{s,0}^2\simeq0.138335,
}
\]

reproducing the accepted hi_class scalar sound speed.

For the bounded-activation future candidate,

\[
{\cal N}_\infty=2.860445,
\]

the released trajectory gives

\[
\boxed{
D_{\min}\simeq0.386661,
}
\]

and

\[
\boxed{
0.12697
\lesssim
c_s^2
\lesssim
0.21447.
}
\]

For the \(F_\infty\)-matched candidate,

\[
{\cal N}_\infty=2.894405,
\]

the corresponding range is

\[
\boxed{
D_{\min}\simeq0.386661,
}
\]

and

\[
\boxed{
0.11745
\lesssim
c_s^2
\lesssim
0.20627.
}
\]

Thus both candidate completions remain ghost-free and gradient-stable in the
full Bellini-Sawicki scalar sector of the released homogeneous trajectory, and
the scalar propagation speed remains subluminal throughout the tested future.

The mature quadratic fixed point gives an especially simple analytic closure.
With

\[
r=\frac{\kappa_2Z_\star}{F_\infty},
\]

one has

\[
\alpha_M,\alpha_B\to0,
\]

and

\[
\boxed{
D_\infty
=
2(1+2r).
}
\]

Since

\[
c_{s,\infty}^2
=
\frac1{1+2r},
\]

the mature branch obeys the exact invariant

\[
\boxed{
D_\infty c_{s,\infty}^2=2.
}
\]

For the two released candidates this predicts

\[
D_\infty\simeq15.7474
\]

and

\[
D_\infty\simeq17.0245,
\]

respectively, exactly matching the numerical approach of the future
trajectories.

The remaining future gate is therefore no longer the sign of \(D\) or
\(c_s^2\) at the homogeneous action level.  It is an implementation-level
cross-check: extend the native hi_class background beyond \(a=1\), keep the
shooting target anchored at the actual present epoch, and verify that hi_class
reproduces the same future \(D\), \(c_s^2\), No-Slip, and structural-lapse
trajectory directly from the log-structural action table.

## 35. Isolated native-hi_class future test

A dedicated future-only infrastructure has now been added without altering any
of the accepted \(z\ge0\) likelihood products.

The experiment performs four controlled operations after all ordinary
present-day and basin tests are complete:

1. extend the log-structural coefficient table into \(\sigma>1\) with the
   \(C^2\) inverse-square tail;
2. extend the background integration endpoint from \(\ln a=0\) to
   \(\ln a=5\);
3. keep the \(\Omega_{\rm smg}\) shooting condition explicitly evaluated at
   \(\ln a=0\), rather than at the new future endpoint;
4. run hi_class in background mode and compare its future
   \(D,c_s^2,\alpha_B+2\alpha_M,p,q,\mathcal N\) against the independent
   homogeneous-action integration.

Both structural endpoint candidates are included:

\[
{\cal N}_\infty=2.860445
\]

and

\[
{\cal N}_\infty=2.894405.
\]

This native future gate is deliberately isolated at the end of the
experimental workflow.  It cannot alter the already generated CMB, matter
power, DESI, supernova, or present-background products.

At the time of this derivation the native future workflow has been launched by
the branch update, but its GitHub Actions result has not yet been independently
retrieved.  Until that run is inspected, the action-level future health result
above is the established result and the native hi_class result remains a
pending implementation cross-check.


## 36. Matter-loaded coasting bridge and an interior 2.869 reference

The asymptotic structural interval

\[
2.860445
\lesssim
{\cal N}_\infty
\lesssim
2.894405
\]

can be viewed from another angle once the mature source has reached

\[
\rho_X\propto a^{-2},
\qquad
w_X=-\frac13,
\qquad
p=1.
\]

At the fluid level, with constant \(F\), conserved matter and radiation give

\[
H^2
=
H_X^2 a^{-2}
\left(
1+\frac{\mu_m}{a}+\frac{\mu_r}{a^2}
\right),
\]

where

\[
\mu_m=\frac{\rho_{m,\rm ref}}{\rho_{X,\rm ref}},
\qquad
\mu_r=\frac{\rho_{r,\rm ref}}{\rho_{X,\rm ref}}.
\]

If the structural mapper has synchronized so that

\[
\sigma\propto a,
\]

then

\[
v=\dot\sigma=\sigma H,
\]

and consequently

\[
\boxed{
\frac{{\cal N}(a)}{{\cal N}_\infty}
=
\sqrt{
1+\frac{\mu_m}{a}+\frac{\mu_r}{a^2}
}.
}
\]

Residual matter and radiation therefore load the structural speed above the
pure-scalar endpoint.

The associated deceleration parameter is

\[
\boxed{
q(a)
=
\frac12
\frac{
\mu_m/a+2\mu_r/a^2
}{
1+\mu_m/a+\mu_r/a^2
}
>0.
}
\]

Equivalently,

\[
\boxed{
\frac{d\ln{\cal N}}{d\ln a}=-q.
}
\]

Thus an exact mature coasting source mixed with residual conserved matter does
not approach \(q=0\) from the accelerating side.  It approaches from the
weakly decelerating side while the matter loading disappears.

This is a fluid-level mature-manifold relation, not yet an exact statement
about the present accepted action.  The current universe still has

\[
w_{X0}\simeq-0.977
\]

rather than \(-1/3\), so applying the mature formula at \(a=1\) is only a
diagnostic extrapolation.

Using the accepted physical densities,

\[
\Omega_{m0}
=
\frac{\omega_b+\omega_c}{h^2}
\simeq0.300756,
\]

\[
\Omega_{r0}
\simeq8.67\times10^{-5},
\]

and

\[
\Omega_{X0}=0.72374015,
\]

gives

\[
\mu_m\simeq0.415557,
\qquad
\mu_r\simeq1.198\times10^{-4}.
\]

If one asks which pure-scalar lapse would reproduce the present
\({\cal N}_0=3.41350846\) through this mature matter-loading formula, the result
is

\[
\boxed{
{\cal N}_{\infty,\rm load}
=
\frac{{\cal N}_0}
{\sqrt{1+\mu_m+\mu_r}}
\simeq2.868923.
}
\]

This value lies naturally inside the independently derived structural bracket,

\[
2.860445
<
2.868923
<
2.894405.
\]

It is about one quarter of the way from the bounded-\(\chi\) endpoint to the
ideal \(\Xi_v\) endpoint.

At the accepted

\[
F_\infty=1.02388543,
\]

this interior lapse corresponds through

\[
{\cal N}_\infty
=
\Xi_v
\sqrt{\frac{\chi_\infty}{F_\infty}}
\]

to

\[
\boxed{
\chi_\infty\simeq1.00594.
}
\]

That is noteworthy because it is only about \(0.6\%\) above unity.  A future
activation interpreted as an effective action-level coupling therefore does
not have to rise all the way to \(F_\infty\simeq1.0239\) in order to land in
the structural range suggested by the present matter loading.

The comparison should not be overinterpreted.  The present action is not yet
on the mature \(p=1,w=-1/3\) manifold, and \(F\), \(\chi\), and the scalar
equation are still evolving.  The \(2.868923\) value is therefore a useful
interior reference, not a newly selected fundamental constant.

It should also be kept distinct from the historically unrelated likelihood
number

\[
\Delta\chi^2_{\rm fair,Planck+DESI}
\simeq2.86998,
\]

whose numerical proximity is coincidental.

The matter-loaded bridge gives a new physical interpretation of why a
\(2.86\)--\(2.89\) pure-scalar endpoint can coexist naturally with a present
structural lapse near \(3.41\): part of the excess structural speed is simply
the finite matter/radiation loading of a universe that has not yet reached its
pure-scalar asymptote.

The sign of the late approach should not, however, be promoted from this
frozen-fluid construction to a fundamental action-level requirement.  In the
actual inverse-square scalar theory, residual matter perturbs the scalar
kinetic state itself.  Section 41 derives that forced response and shows that
the healthy covariant fixed point instead has a universal leading matter
correction with \(q\to0^-\).  The \(2.868923\) value therefore remains a
useful interior normalization diagnostic, while the sign of the asymptotic
approach must be taken from the scalar field equation rather than from the
frozen-source bridge.


## 37. Perturbation-ready C3 stitching and near-exact future No-Slip

The previous future tails matched the accepted action through second field
derivatives at the present boundary.  That is sufficient for the homogeneous
equations, but it is not the correct final smoothness requirement for a native
hi_class perturbation implementation.

The hi_class gravity-function module contains explicit
\(G_{3,\phi\phi\phi}\) and \(G_{4,\phi\phi\phi}\) terms.  Therefore the future
continuation of the present action should preserve the complete third-order
field jet of \(G_3\) and \(G_4\).  In the structural variable
\(x=\ln\sigma\), the future \(F\) and \(g\) tails are consequently upgraded to

\[
Q(x)
=
Q_\infty
+
e^{-\mu x}
\left(
A_0+A_1x+\frac{A_2}{2}x^2+\frac{A_3}{6}x^3
\right).
\]

For a prescribed present jet

\[
q_n
\equiv
\left.
\frac{d^nQ}{dx^n}
\right|_{x=0},
\qquad n=0,1,2,3,
\]

the coefficients are fixed recursively by

\[
\boxed{
A_n
=
q_n-\delta_{n0}Q_\infty
-
\sum_{j=0}^{n-1}
{n\choose j}
(-\mu)^{\,n-j}A_j.
}
\]

This guarantees continuity of

\[
Q,\quad Q_{,x},\quad Q_{,xx},\quad Q_{,xxx}
\]

through the present boundary.

For the accepted structural reconstruction the present \(x\)-derivatives are
approximately

\[
F(0)=1.02346657,
\]

\[
F_{,x}(0)=1.07891\times10^{-3},
\]

\[
F_{,xx}(0)=-2.13347\times10^{-3},
\]

\[
F_{,xxx}(0)=3.47914\times10^{-3}.
\]

The transformed cubic-braiding coefficient likewise has a finite third-order
jet and is matched at the same order.

The asymptotic No-Slip matching condition is modified only in its leading
polynomial power.  With the C3 tails,

\[
F-F_\infty
\sim
\frac{A_3^{(F)}}{6}x^3e^{-\mu_Fx},
\]

while

\[
g
\sim
\frac{A_3^{(g)}}{6}x^3e^{-\mu_gx}.
\]

Since

\[
F_{,\sigma}
=
e^{-x}F_{,x},
\]

asymptotic No-Slip requires

\[
\boxed{
\mu_g=\mu_F+1,
}
\]

and

\[
\boxed{
Z_\infty
=
\frac{
\mu_F A_3^{(F)}
}{
2A_3^{(g)}
}.
}
\]

The fast positive C3 branches selected by the two mature structural endpoints
are

\[
\boxed{
\mu_F\simeq5.55458
}
\]

for the bounded-\(\chi\) candidate and

\[
\boxed{
\mu_F\simeq5.57323
}
\]

for the \(F_\infty\)-matched candidate.

The released trajectory can then be used to refine \(g\) toward the exact
on-trajectory relation

\[
g_{\rm NS}(\sigma)
=
-\frac{F_{,\sigma}}{2Z(\sigma)}.
\]

To preserve the complete C3 jet of the accepted branch, the refinement is
multiplied by

\[
\boxed{
W(u)
=
1-e^{-u}
\left(
1+u+\frac{u^2}{2}+\frac{u^3}{6}
\right),
}
\]

with

\[
u=\lambda_{\rm NS}\ln\sigma.
\]

This satisfies

\[
W(0)=W'(0)=W''(0)=W'''(0)=0.
\]

A two-iteration release/refinement audit with
\(\lambda_{\rm NS}=500\) gives

\[
\boxed{
\max|\alpha_B+2\alpha_M|
\simeq7.6\times10^{-8}
}
\]

for the bounded-activation branch and

\[
\boxed{
\max|\alpha_B+2\alpha_M|
\simeq7.1\times10^{-8}
}
\]

for the \(F_\infty\)-matched branch.

The improved No-Slip matching leaves the background dynamics and scalar-health
intervals essentially unchanged.  In the same local audit the bounded branch
retains

\[
D_{\min}\simeq0.38666113,
\]

\[
0.12697117
\lesssim
c_s^2
\lesssim
0.21447760,
\]

while the \(F_\infty\)-matched branch retains

\[
D_{\min}\simeq0.38666113,
\]

\[
0.11744554
\lesssim
c_s^2
\lesssim
0.20628612.
\]

The C3 upgrade therefore removes a perturbation-level smoothness concern
without changing the mature structural endpoint.

## 38. Recovery of the No-Slip gravitational sector

The hi_class Horndeski effective-gravity variables make the physical meaning of
the refined future branch especially transparent.

For \(c_T^2=1\), define

\[
\beta_1
=
\alpha_B+2\alpha_M.
\]

The quasi-static effective Newton coupling used by hi_class is

\[
G_{\rm eff}
=
\frac1F
\left[
1-
\frac{
\alpha_B\beta_1
}{
\alpha_B\beta_1-\beta_2
}
\right],
\]

with \(\beta_2\) the corresponding background combination.

Thus exact No-Slip,

\[
\beta_1=0,
\]

immediately implies

\[
\boxed{
G_{\rm eff}=\frac1F.
}
\]

The gravitational slip simultaneously tends to

\[
\boxed{
\eta_{\rm slip}=1.
}
\]

The C3-refined future audit gives

\[
\max|F G_{\rm eff}-1|
\lesssim1.45\times10^{-9}
\]

for the bounded branch and

\[
\max|F G_{\rm eff}-1|
\lesssim1.35\times10^{-9}
\]

for the \(F_\infty\)-matched branch.

The maximum slip departures are of the same order,

\[
\boxed{
\max|\eta_{\rm slip}-1|
\sim10^{-9}.
}
\]

At the mature endpoint,

\[
F_\infty
=
1.02388543,
\]

so

\[
\boxed{
G_{{\rm eff},\infty}
=
F_\infty^{-1}
=
0.97667178
}
\]

in the current hi_class normalization.

At the accepted present endpoint,

\[
F_0=1.02346657,
\]

which corresponds to

\[
F_0^{-1}=0.97707148.
\]

Therefore the remaining evolution of the effective gravitational normalization
from the present branch to the mature branch is only

\[
\boxed{
\frac{
G_{{\rm eff},\infty}
}{
G_{{\rm eff},0}
}
-1
\simeq
-4.09\times10^{-4},
}
\]

or about

\[
\boxed{
-0.0409\%.
}
\]

This is a useful structural result.  The late transition in
\({\cal N}\) from \(3.4135\) toward \(2.86\)--\(2.89\) is not accompanied by a
large residual modification of gravitational slip or the effective Newton
sector.  The scalar becomes dynamically important through its structural
energy while the metric coupling itself approaches a nearly frozen No-Slip
limit.

## 39. Future geometric regularity of the mature coasting attractor

The inverse-square fixed point also gives a simple global future geometry.

On the mature branch,

\[
p\to1,
\qquad
a(t)\propto t-t_B,
\]

so

\[
H
=
\frac1{t-t_B},
\]

and

\[
\dot H
=
-\frac1{(t-t_B)^2}.
\]

Therefore

\[
\boxed{
q
=
-1-\frac{\dot H}{H^2}
=
0.
}
\]

For a spatially flat FLRW metric the Ricci scalar is

\[
R_{\rm Ricci}
=
6(2H^2+\dot H),
\]

so the mature solution gives

\[
\boxed{
R_{\rm Ricci}
=
\frac{6}{(t-t_B)^2}
\rightarrow0.
}
\]

The Kretschmann scalar is

\[
{\cal K}
=
12
\left[
(H^2+\dot H)^2+H^4
\right].
\]

Because \(H^2+\dot H=0\) in exact coasting,

\[
\boxed{
{\cal K}
=
\frac{12}{(t-t_B)^4}
\rightarrow0.
}
\]

Thus the mature fixed point contains no finite-future curvature blow-up; the
curvature invariants decay toward zero.

The future conformal interval is

\[
\Delta\eta_{\rm future}
=
\int_t^\infty
\frac{dt'}{a(t')}.
\]

For

\[
a(t')\propto t'-t_B,
\]

this becomes

\[
\Delta\eta_{\rm future}
\propto
\int_t^\infty
\frac{dt'}{t'-t_B},
\]

which diverges logarithmically.  Hence the exact mature coasting limit has

\[
\boxed{
\Delta\eta_{\rm future}=\infty.
}
\]

In particular, the asymptotic coasting branch does not possess the finite
future event horizon characteristic of an eternal de-Sitter phase.

The assumed scalar domination is also self-consistent.  In the mature branch,

\[
\rho_X\propto a^{-2},
\]

whereas

\[
\rho_m\propto a^{-3},
\qquad
\rho_r\propto a^{-4}.
\]

Therefore

\[
\frac{\rho_m}{\rho_X}\propto a^{-1}\to0,
\]

and

\[
\frac{\rho_r}{\rho_X}\propto a^{-2}\to0.
\]

At the same time,

\[
\frac{K_m}{K_{m0}}
=
\frac{K_p}{K_{p0}}
=
\frac{\sigma}{a}
\]

tends to a constant because the mature attractor has
\(\sigma\propto a\).

The future completion is therefore internally coherent at the background
level:

\[
\boxed{
\text{stable inverse-square scalar domination}
\;\Longrightarrow\;
q\to0,\;
D>0,\;
c_s^2>0,\;
\eta_{\rm slip}\to1,\;
G_{\rm eff}\to F_\infty^{-1},
}
\]

with decaying curvature, frozen matter/baryon mapping ratios, and no finite
future event horizon.

The remaining numerical problem is no longer the existence or background
health of the mature endpoint.  It is to propagate the actual hi_class scalar
perturbation variables beyond the present boundary using this C3 future action,
without redefining the observationally calibrated \(z\ge0\) branch.


## 40. Literal bounded-activation stress test

The smooth future tail used in Sections 29--35 was designed primarily to test
whether the accepted present action could reach the mature inverse-square
coasting basin.  It succeeds dynamically, but the action-level activation
history exposes an additional distinction between two interpretations of
\(\chi\).

Using

\[
\chi(\sigma)
=
\chi_0
\frac{\rho_X(\sigma)}{\rho_{X0}}
\sigma^2,
\]

the original bounded-endpoint tail with

\[
\mu_Q\simeq3.75024
\]

reaches a transient maximum

\[
\boxed{
\chi_{\max}\simeq1.96856
}
\]

near

\[
\ln a\simeq0.778.
\]

It later relaxes to the intended asymptotic normalization.  Therefore this
slow tail is consistent if \(\chi\) is interpreted as an effective
action-level coupling, but it is not compatible with treating \(\chi\) as
a literal fraction constrained to remain below unity at every future epoch.

A dedicated stress test was therefore performed with the same accepted
present endpoint and the same mature coasting target, but with progressively
faster \(C^2\) flow of the \(G_2\) coefficients toward their inverse-square
normal form.

At fixed

\[
r=3.43684444,
\]

representative results are

\[
\begin{array}{c|c|c}
\mu_Q & \chi_{\max} & \chi_{\max}/\chi_{\infty,\rm rec}-1\\
\hline
3.75024 & 1.96856 & 0.96560\\
20 & 1.15519 & 0.15346\\
40 & 1.05219 & 0.05062\\
60 & 1.01644 & 0.01492\\
68 & 1.00789 & 0.00638\\
72 & 1.00430 & 0.00280\\
76 & 1.00150 & 0\\
80 & 1.00150 & 0
\end{array}
\]

where \(\chi_{\infty,\rm rec}\simeq1.00150\) is the common late
normalization reached by the reconstructed action.  The approximately
\(1.5\times10^{-3}\) offset from the nominal target \(\chi_\infty=1\)
is at the same level as the endpoint reconstruction normalization mismatch and
is therefore kept separate from a resolved transient overshoot.

The important result is that by

\[
\boxed{
\mu_Q\simeq76
}
\]

the future history no longer develops a resolved activation peak above its own
asymptotic level.

This fast-flow branch remains healthy in the exact Bellini-Sawicki scalar
sector.  For \(\mu_Q=76\),

\[
\boxed{
D_{\min}\simeq0.386661>0,
}
\]

and

\[
\boxed{
0.11943
\lesssim
c_s^2
\lesssim
0.20866.
}
\]

Thus the bounded-activation requirement and scalar stability are not mutually
exclusive.

The price is a much sharper future handoff.  The first acceleration-to-
deceleration crossing moves from approximately

\[
\ln a\simeq1.13
\]

for the slow tail to approximately

\[
\boxed{
\ln a\simeq0.063
}
\]

for the \(\mu_Q=76\) branch, while the second crossing occurs near

\[
\ln a\simeq1.50.
\]

The matter mapper also freezes at a different future normalization.  Since

\[
\frac{K_m}{K_{m0}}
=
\frac{K_p}{K_{p0}}
=
\frac{\sigma}{a},
\]

the slow tail approaches approximately

\[
\frac{K_{m,\infty}}{K_{m0}}
\simeq0.713,
\]

whereas the fast bounded-activation tail approaches approximately

\[
\boxed{
\frac{K_{m,\infty}}{K_{m0}}
\simeq0.893.
}
\]

The corresponding baryonic mapper has the same fractional evolution.

There is also a No-Slip tradeoff.  The original smooth tail has

\[
\max|\alpha_B+2\alpha_M|
\simeq5.7\times10^{-5},
\]

whereas the rapid \(\mu_Q=76\) branch reaches approximately

\[
\boxed{
2.45\times10^{-4}.
}
\]

This remains numerically small, but it is larger because the kinetic trajectory
changes before the independently stitched \(F\) and \(g\) tails can fully
readjust.  A future production model that insists simultaneously on

\[
\chi\le1,
\qquad
|\alpha_B+2\alpha_M|<10^{-5},
\]

will therefore require more than a single common \(G_2\) transition rate.
The natural next step is a coupled future reconstruction in which
\(k_1,k_2,V,F\), and \(g\) are solved together against the activation and
No-Slip conditions rather than stitched independently.

This stress test therefore sharpens the interpretation fork:

- a smooth future continuation naturally favors \(\chi\) as an effective
  coupling that may temporarily exceed unity;
- a literal bounded fraction remains possible, but it requires a rapid
  future action flow and a jointly reconstructed No-Slip sector.

The branch now records both possibilities rather than assuming either one in
advance.


## 41. Action-level matter-forced approach to the mature fixed point

The fluid-level construction in Section 36 deliberately froze the mature
structural source while adding conserved matter.  The covariant scalar does
not remain frozen in this way.  Once the inverse-square normal form has been
reached, residual matter perturbs the kinetic state \(Z\), and the scalar
responds dynamically.

Take

\[
G_2(\sigma,Z)=\sigma^{-2}f(Z),
\qquad
G_3\to0,
\qquad
F\to F_\infty=\text{constant}.
\]

Define

\[
{\cal R}(Z)
\equiv
2Zf_{,Z}-f,
\]

so that

\[
\rho_X=\frac{{\cal R}(Z)}{\sigma^2}.
\]

The homogeneous scalar equation is

\[
\frac{d}{dt}
\left(
\sigma^{-2}f_{,Z}\dot\sigma
\right)
+
3H\sigma^{-2}f_{,Z}\dot\sigma
+
2\sigma^{-3}f
=
0.
\]

Writing

\[
p\equiv\frac{d\ln\sigma}{d\ln a}
=
\frac{\dot\sigma}{H\sigma},
\]

and using

\[
Z=\frac12\dot\sigma^2,
\]

this becomes the exact first-order relation

\[
\boxed{
{\cal K}\frac{dZ}{d\ln a}
+
6Zf_{,Z}
-
2p{\cal R}
=
0,
}
\]

where

\[
{\cal K}
=
f_{,Z}+2Zf_{,ZZ}.
\]

The flat Friedmann equation gives independently

\[
\boxed{
\frac{6FZ}{p^2}
=
\rho_m\sigma^2
+
\rho_r\sigma^2
+
{\cal R}(Z).
}
\]

These two equations show why the frozen-fluid bridge is not the exact
action-level late solution.  If one imposed simultaneously

\[
p=1,
\qquad
\rho_X\propto a^{-2},
\]

then \(\sigma\propto a\), and therefore
\({\cal R}(Z)\) would have to remain constant.  But

\[
\frac{d{\cal R}}{dZ}
=
{\cal K}.
\]

For a healthy scalar,

\[
{\cal K}>0,
\]

so constant \({\cal R}\) forces

\[
Z=\text{constant}.
\]

Then \(p=1\) implies exact coasting,

\[
H\propto a^{-1},
\]

which cannot simultaneously contain a finite separately conserved
\(\rho_m\propto a^{-3}\) contribution.  Thus the exact
\(p=1,\ w_X=-1/3\) state is reached only asymptotically as matter disappears.

The leading forced correction can be derived analytically.  Let

\[
Z=Z_\star(1+z),
\qquad
p=1+u,
\]

and define the dimensionless matter loading

\[
m
\equiv
\frac{\rho_m\sigma^2}{{\cal R}_\star}.
\]

At the pure-scalar fixed point,

\[
{\cal R}_\star=6FZ_\star,
\qquad
f_{,Z}(Z_\star)=2F.
\]

Introduce

\[
\boxed{
\kappa
\equiv
\frac{{\cal K}_\star}{2F}.
}
\]

For the quadratic mature family,

\[
\kappa=1+2r
=
\frac{1}{c_{s,\star}^2}.
\]

To first order in \(m\), the Friedmann constraint gives

\[
-2u
=
m
+
\left(
\frac{\kappa}{3}-1
\right)z.
\]

The scalar equation gives

\[
\frac{\kappa}{3}\frac{dz}{d\ln a}
+
\left(
1+\frac{\kappa}{3}
\right)z
-
2u
=
0.
\]

Eliminating \(u\),

\[
\boxed{
\frac{dz}{d\ln a}
+
2z
=
-\frac{3}{\kappa}m.
}
\]

At late times,

\[
\frac{dm}{d\ln a}
=
-m+O(m^2),
\]

so the solution is

\[
\boxed{
z
=
-\frac{3}{\kappa}m
+
C a^{-2}.
}
\]

This separates two physical relaxation modes.  The intrinsic scalar
perturbation found in Section 25 decays as

\[
a^{-2},
\]

whereas the matter-forced correction decays only as

\[
a^{-1}.
\]

Consequently the matter-forced term eventually dominates the final approach
to the fixed point even though matter itself becomes negligible.

Substituting the forced solution back into the Friedmann relation gives

\[
\boxed{
p-1
=
-\frac{3}{2\kappa}m
+
O(a^{-2}),
}
\]

and because

\[
{\cal N}\propto\sqrt{Z},
\]

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-\frac{3}{2\kappa}m
+
O(a^{-2}).
}
\]

The rescaled scalar density responds as

\[
\frac{{\cal R}}{{\cal R}_\star}-1
=
\frac{\kappa}{3}z
=
-m+O(a^{-2}).
\]

Thus, to first order, the scalar density decreases by precisely the amount
introduced by the residual matter term.  The direct matter loading and the
scalar response cancel in the Friedmann amplitude at \(O(m)\).

As a result,

\[
H\sigma
=
\text{constant}
+
O(m^2),
\]

and therefore

\[
\boxed{
q
=
-\frac{3}{2\kappa}m
+
O(a^{-2}).
}
\]

The exact inverse-square scalar consequently approaches coasting from the
very weakly accelerating side,

\[
\boxed{
q\to0^-,
}
\]

rather than the \(q\to0^+\) behavior of the frozen-fluid diagnostic.

The total effective equation of state is correspondingly

\[
\boxed{
w_{\rm eff}
=
-\frac13
-
\frac{m}{\kappa}
+
O(a^{-2}),
}
\]

while the scalar itself obeys

\[
\boxed{
w_X
=
-\frac13
-
\left(
\frac13+\frac1\kappa
\right)m
+
O(a^{-2}).
}
\]

These corrections remain safely above the phantom boundary for sufficiently
late \(m\ll1\).

The freely evolved future tails provide a direct numerical check.  For the
bounded-\(\chi\) candidate,

\[
r=3.43684444,
\qquad
\kappa=7.87368889.
\]

At \(\ln a=10\), the released action has approximately

\[
m=9.02\times10^{-6},
\]

for which the asymptotic formula predicts

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
\simeq
-1.72\times10^{-6}.
\]

The numerical trajectory gives

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
\simeq
-1.71\times10^{-6},
\]

and

\[
q\simeq-1.68\times10^{-6}.
\]

For the \(F_\infty\)-matched candidate,

\[
r=3.75611884,
\qquad
\kappa=8.51223767,
\]

and at the same epoch

\[
m=8.76\times10^{-6}.
\]

The analytic prediction is

\[
-\frac{3m}{2\kappa}
\simeq
-1.54\times10^{-6},
\]

while the released action gives approximately

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
\simeq
-1.53\times10^{-6},
\]

and

\[
q\simeq-1.50\times10^{-6}.
\]

The agreement is at the percent level by \(\ln a=10\).

This explains a feature that previously looked like a possible interpolation
artifact: after the finite future tail passes through a weakly decelerating
phase, \(q\) crosses back below zero at very large scale factor.  The final
negative sign is in fact the expected matter-forced asymptotic response of the
healthy inverse-square scalar.

The mature endpoint itself is unchanged,

\[
q\to0,
\qquad
p\to1,
\qquad
{\cal N}\to{\cal N}_\infty,
\]

but the direction from which the action approaches that endpoint is now
derived rather than guessed.

The interior value

\[
{\cal N}_{\infty,\rm load}\simeq2.868923
\]

from Section 36 remains useful as a normalization diagnostic because it shows
that the present \(3.4135\) lapse and the \(2.86\)--\(2.89\) pure-scalar bracket
are of the right relative scale.  It should not, however, be interpreted as
the exact dynamical matter correction to the mature scalar action.


## 34. Scalar-curvature envelope through the future transition

The asymptotic fixed-point argument can be extended through the finite
present-to-coasting transition without replacing the full hi_class hierarchy.

For Horndeski scalar perturbations define

\[
\boxed{
Q_s
=
\frac{2FD}{(2-\alpha_B)^2}.
}
\]

The curvature-mode envelope obeys the standard action-level equation

\[
\boxed{
\zeta_{,NN}
+
\left(
3+\frac{d\ln H}{dN}
+\frac{d\ln Q_s}{dN}
\right)\zeta_{,N}
+
c_s^2
\left(
\frac{k}{aH}
\right)^2
\zeta
=
0,
}
\]

where here \(N=\ln a\).

For the C3 No-Slip-refined future branches, the exact background
Bellini-Sawicki quantities remain positive through the tested interval.  The
bounded-activation branch gives approximately

\[
0.19744
\lesssim
Q_s
\lesssim
8.0581,
\]

with

\[
\min
\left[
3+\frac{d\ln H}{dN}
+\frac{d\ln Q_s}{dN}
\right]
\simeq
1.9994.
\]

The \(F_\infty\)-matched branch gives

\[
0.19744
\lesssim
Q_s
\lesssim
8.7119,
\]

with minimum friction approximately

\[
2.0004.
\]

Thus neither branch develops a negative scalar kinetic normalization or a
negative curvature-mode friction interval in this reduced envelope test.

Representative modes were initialized at the present boundary with

\[
\zeta(0)=1,
\qquad
\zeta_{,N}(0)=0,
\]

and present-horizon ratios

\[
\nu_0\equiv\frac{k}{H_0}
=
0,\ 0.1,\ 1,\ 10.
\]

Across six future e-folds the maximum absolute amplitude did not exceed its
initial value in any tested case.

For the bounded branch the final absolute amplitudes are approximately

\[
|\zeta_f|
=
1,\quad
0.99776,\quad
0.79572,\quad
1.14\times10^{-3},
\]

for \(\nu_0=0,0.1,1,10\), respectively.

For the \(F_\infty\)-matched branch they are approximately

\[
|\zeta_f|
=
1,\quad
0.99796,\quad
0.81277,\quad
7.02\times10^{-4}.
\]

These numbers are not a substitute for the multi-species Boltzmann hierarchy.
They show that the scalar-curvature degree of freedom of the released action
does not acquire a transient growing envelope in the representative band that
was tested.

At the exact mature fixed point,

\[
\alpha_B=\alpha_M=0,
\]

\[
D=2(1+2r),
\]

\[
c_s^2=\frac1{1+2r},
\]

and

\[
Q_s=F_\infty(1+2r)=\text{constant}.
\]

Exact coasting has

\[
\frac{d\ln H}{dN}=-1
\]

and

\[
aH=\text{constant}.
\]

Therefore the mode equation becomes

\[
\boxed{
\zeta_{,NN}
+
2\zeta_{,N}
+
c_s^2
\left(
\frac{k}{aH}
\right)^2
\zeta
=
0.
}
\]

Writing

\[
\zeta\propto a^s
\]

gives

\[
\boxed{
s
=
-1
\pm
\sqrt{
1-c_s^2
\left(
\frac{k}{aH}
\right)^2
}.
}
\]

For \(k/(aH)\to0\),

\[
s=0,-2,
\]

so the two mature super-horizon modes are a constant mode and an
\(a^{-2}\) decaying mode.

For

\[
c_s^2
\left(
\frac{k}{aH}
\right)^2>1,
\]

the roots are complex with real part \(-1\), giving an oscillatory envelope

\[
\boxed{
|\zeta|\propto a^{-1}.
}
\]

Hence the mature inverse-square attractor possesses no growing
scalar-curvature mode in this sector.

## 35. Native future perturbation propagation as a separate code gate

The standard hi_class time bookkeeping contains an important distinction that
must be preserved.

For the observational calculation, the quantity

\[
\texttt{conformal\_age}
\]

means the conformal time at the physical present epoch.  The isolated future
background patch deliberately keeps that quantity anchored at \(a=1\), even
while the background integration table is extended to \(a>1\).  This prevents
the future audit from corrupting present distances, growth normalization, and
shooting targets.

The same choice means that the stock perturbation module still stops its source
sampling at the present epoch.  Therefore a future perturbation audit cannot
simply reuse the extended background table and call itself complete.

A separate future-only perturbation endpoint has now been introduced for the
experimental branch:

\[
\boxed{
\tau_{\rm end}^{\rm pert}
=
\tau_{\rm bg}(a_{\rm max})
}
\]

for non-CMB diagnostic runs, while CMB calculations remain anchored to the
ordinary observational conformal age.

There is a second implementation issue.  The standard thermodynamics
interpolation table is defined for

\[
z\ge0.
\]

Future perturbation evolution reaches

\[
-1<z<0.
\]

For the isolated diagnostic run only, a controlled post-today continuation is
therefore used.  It loads the \(z=0\) ionization state and extrapolates

\[
\dot\kappa
\propto
a^{-2},
\]

\[
T_b\propto a^{-2},
\]

\[
w_b\propto a^{-2},
\qquad
c_b^2\propto a^{-2},
\]

while future visibility-source terms are set to zero.

This continuation is not being promoted as a prediction of future atomic or
astrophysical thermodynamics.  Its purpose is narrower: after photons have
decoupled and baryon pressure is dynamically negligible, it prevents an
unphysical interpolation failure at negative redshift while allowing the
metric, matter, radiation, and Horndeski scalar perturbation equations to be
propagated through the future structural transition.

The native test is configured at representative wavenumbers

\[
k
=
10^{-4},\ 10^{-3},\ 10^{-2}\ {\rm Mpc}^{-1},
\]

with the background extended to

\[
\ln a=5.
\]

The resulting perturbation tables expose the native hi_class scalar variables,
including the dynamical scalar-field perturbation and its first two time
derivatives, together with the metric and standard species perturbations.

This creates a clean hierarchy of future tests:

\[
\boxed{
\text{analytic fixed point}
\rightarrow
\text{exact homogeneous action}
\rightarrow
\text{Bellini-Sawicki health}
\rightarrow
\text{curvature envelope}
\rightarrow
\text{native hi_class perturbation propagation}.
}
\]

The first four stages are already closed for the tested future candidates.  The
last stage is now implemented as an isolated workflow gate and should be
treated as pending until its native run artifacts are available and inspected.


## 42. Exact Planck identity behind the 2.894405 coefficient

The density-route coefficient introduced above is not merely numerically close
to the earlier manuscript value \(\Xi_v=\sqrt{8\pi/3}\).  With the standard
Planck definitions used by SDMC,

\[
t_P^2=\frac{\hbar G}{c^5},
\qquad
\rho_P=\frac{c^5}{\hbar G^2},
\]

where \(\rho_P\) is the Planck mass density.  Their product is exactly

\[
\boxed{
\rho_P t_P^2=\frac1G.
}
\]

Therefore

\[
\Xi_P^2
=
\frac{8\pi G}{3}\rho_Pt_P^2
=
\frac{8\pi}{3},
\]

and hence

\[
\boxed{
\Xi_P
=
\Xi_v
=
\sqrt{\frac{8\pi}{3}}
=
2.894405018\ldots
}
\]

identically.

The small \(2.89445\)-level discrepancy that appeared when rounded numerical
values of \(G,t_P,\rho_P\) were inserted separately is consequently only a
rounding artifact.  The exact structural coefficient is \(2.894405018\ldots\).

The mature lapse formula can therefore be written without any independent
Planck-density normalization constant:

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
\sqrt{\frac{\chi_\infty}{F_\infty}}.
}
\]

For the accepted

\[
F_\infty=e^{0.02360463334}=1.02388543,
\]

this gives the two endpoint normalizations already found numerically,

\[
\boxed{
{\cal N}_\infty=2.86044513
\quad
(\chi_\infty=1),
}
\]

and

\[
\boxed{
{\cal N}_\infty=2.89440502
\quad
(\chi_\infty=F_\infty).
}
\]

Thus the remembered \(2.86\)--\(2.89\) interval is anchored at its upper end by
an exact Planck identity, not by a fit parameter.  What the future covariant
action must determine is the ratio \(\chi_\infty/F_\infty\), which selects the
location inside that interval.


## 43. Exact Lagrangian evolution equation for the structural lapse

The structural lapse no longer needs to be treated as an external closure once
the structural field is used directly.

With

\[
\sigma=\frac{S}{S_0},
\qquad
v\equiv\dot\sigma,
\]

the lapse is

\[
\boxed{
{\cal N}=t_P S_0 v.
}
\]

Hence

\[
\boxed{
\frac{d\ln{\cal N}}{d\ln a}
=
\frac{\dot v}{Hv}.
}
\]

Using the exact homogeneous system of Section 27,

\[
\Delta_{\rm hom}
=
2F{\cal K}+3Q^2,
\]

\[
\dot v
=
-\frac{
2F{\cal C}+3Q{\cal R}
}{
\Delta_{\rm hom}
},
\]

the structural-clock equation becomes

\[
\boxed{
\frac{d\ln{\cal N}}{d\ln a}
=
-
\frac{
2F{\cal C}+3Q{\cal R}
}{
Hv\left(2F{\cal K}+3Q^2\right)
}.
}
\]

Equivalently,

\[
\boxed{
\frac{d{\cal N}}{d\ln a}
=
-
\frac{
t_P S_0
\left(2F{\cal C}+3Q{\cal R}\right)
}{
H\left(2F{\cal K}+3Q^2\right)
}.
}
\]

This is the direct Lagrangian evolution law for the SDMC lapse.  No relation of
the form

\[
S/S_0=t/t_0
\]

is required after the initial endpoint is specified.

The bridge exponent

\[
p
=
\frac{d\ln\sigma}{d\ln a}
=
\frac{v}{H\sigma}
\]

is also a dynamical output.  Its exact first-order evolution is

\[
\boxed{
\frac{d\ln p}{d\ln a}
=
\frac{d\ln{\cal N}}{d\ln a}
+
1+q-p.
}
\]

The deceleration parameter follows from the same homogeneous matrix,

\[
\boxed{
q
=
-1
-
\frac{
Q{\cal C}-{\cal K}{\cal R}
}{
H^2\left(2F{\cal K}+3Q^2\right)
}.
}
\]

Thus

\[
{\cal N},\qquad p,\qquad q
\]

are not three independent closures.  They are three projections of the same
two-dimensional homogeneous scalar-metric system.

The positivity of

\[
\Delta_{\rm hom}=2F{\cal K}+3Q^2
\]

along the released healthy branch is particularly useful.  It shows that the
clock equation remains nonsingular whenever the effective Planck mass and
homogeneous scalar kinetic sector remain healthy.

## 44. Exact lapse-mapper identity and the late-time cancellation

The matter and baryon mappers derived earlier obey

\[
\boxed{
M
\equiv
\frac{K_m}{K_{m0}}
=
\frac{K_p}{K_{p0}}
=
\frac{\sigma}{a}.
}
\]

Therefore

\[
\boxed{
\frac{d\ln M}{d\ln a}=p-1.
}
\]

Combining

\[
{\cal N}=t_PS_0v,
\qquad
p=\frac{v}{H\sigma},
\qquad
M=\frac{\sigma}{a},
\]

gives the exact identity

\[
\boxed{
{\cal N}
=
t_PS_0\,p\,M\,aH.
}
\]

This relation is useful because it unifies the structural clock and the domain
mapper in one equation.  It also provides an immediate differential
consistency condition,

\[
\frac{d\ln{\cal N}}{d\ln a}
=
\frac{d\ln p}{d\ln a}
+
(p-1)
-q.
\]

The right-hand side is identical to the direct lapse equation after using the
exact \(p\)-flow above.

At the mature matter-forced fixed point of Section 41 define

\[
A\equiv\frac{3}{2\kappa},
\qquad
\kappa=1+2r.
\]

The action gives

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-Am+O(a^{-2}),
\]

\[
p-1
=
-Am+O(a^{-2}),
\]

and

\[
q
=
-Am+O(a^{-2}),
\]

where

\[
m\propto a^{-1}
\]

is the residual matter loading.

Since

\[
\frac{d\ln M}{d\ln a}
=
-Am+O(a^{-2})
\]

and

\[
\frac{dm}{d\ln a}
=
-m+O(m^2),
\]

integration toward the fixed point gives

\[
\boxed{
\ln\frac{M}{M_\infty}
=
Am+O(a^{-2}),
}
\]

or

\[
\boxed{
\frac{M}{M_\infty}
=
1+Am+O(a^{-2}).
}
\]

Thus the matter mapper approaches its future constant from the opposite side
to the structural lapse:

\[
\boxed{
\frac{K_m}{K_{m,\infty}}-1
=
\frac{K_p}{K_{p,\infty}}-1
=
+\frac{3m}{2\kappa}
+O(a^{-2}),
}
\]

whereas

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-\frac{3m}{2\kappa}
+O(a^{-2}).
}
\]

The leading matter correction therefore cancels in their product,

\[
\boxed{
\frac{{\cal N}K_m}
{{\cal N}_\infty K_{m,\infty}}
=
1+O(a^{-2}),
}
\]

and identically

\[
\boxed{
\frac{{\cal N}K_p}
{{\cal N}_\infty K_{p,\infty}}
=
1+O(a^{-2}).
}
\]

This is a new structural invariant of the mature action: the lapse and the
matter-domain mapper each retain an \(O(a^{-1})\) matter correction, but with
equal magnitude and opposite sign.

The C3 No-Slip-refined numerical trajectories reproduce this cancellation.

For the bounded-\(\chi\) branch,

\[
\kappa=7.87368889.
\]

At

\[
\ln a=10
\]

the matter loading gives

\[
Am\simeq1.71818\times10^{-6}.
\]

The numerical lapse satisfies

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
\simeq
-1.70644\times10^{-6},
\]

while the inferred mapper correction is

\[
\frac{M}{M_\infty}-1
\simeq
+1.71818\times10^{-6}.
\]

Their normalized product differs from unity by only

\[
\boxed{
1.17\times10^{-8}.
}
\]

For the \(F_\infty\)-matched branch,

\[
\kappa=8.51223767,
\]

with

\[
Am\simeq1.54360\times10^{-6}
\]

at the same epoch.  The numerical lapse gives

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
\simeq
-1.53304\times10^{-6},
\]

while

\[
\frac{M}{M_\infty}-1
\simeq
+1.54360\times10^{-6}.
\]

The normalized lapse-mapper product differs from unity by only

\[
\boxed{
1.06\times10^{-8}.
}
\]

The cancellation is already visible at \(\ln a=8\), but the larger
\(a^{-2}\) intrinsic scalar transient is still measurable there.  By
\(\ln a=10\), the first-order matter duality is the dominant description.

This result also explains why

\[
H\sigma
\]

is asymptotically constant at first order.  Since

\[
\frac{{\cal N}}{p}
=
t_PS_0H\sigma,
\]

and both \({\cal N}/{\cal N}_\infty\) and \(p\) contain the same
\(-Am\) correction, their ratio has no \(O(a^{-1})\) term.

## 45. Native future perturbation gate is now closed

The native future perturbation gate described earlier as pending has now been
run successfully on the isolated future-capable hi_class branch.

For both mature endpoint candidates, the native background reaches

\[
\ln a=5,
\qquad
a=148.413159\ldots,
\]

and the scalar perturbation hierarchy is propagated at

\[
k=
10^{-4},
\quad
10^{-3},
\quad
10^{-2}
\ {\rm Mpc}^{-1}.
\]

All stored perturbation columns remain finite through the future endpoint and
all six candidate-mode endpoint gates pass.  The workflow returns

\[
\boxed{
\texttt{NATIVE\_FUTURE\_GATE\_PASS}.
}
\]

The metric potentials remain bounded in the tested interval.  For the
\(F_\infty\)-matched branch, for example, the \(k=10^{-2}\,{\rm Mpc}^{-1}\)
mode has

\[
|\phi(a=148.4)|
\simeq3.66\times10^{-3},
\]

compared with approximately

\[
0.312
\]

at the first future sample.  The corresponding \(\psi\) behavior is
essentially identical.

For the longest tested wavelength,

\[
k=10^{-4}\ {\rm Mpc}^{-1},
\]

the metric potentials remain of order unity relative to their present
normalization rather than developing a runaway mode.  The scalar field
perturbation itself remains finite in every table.

The same statement holds for the bounded-\(\chi\) endpoint.

The matter density contrasts continue to evolve in synchronous gauge and can
grow in absolute coordinate amplitude.  That behavior is not by itself a
future instability criterion.  The decisive numerical statements from this
gate are narrower:

\[
\boxed{
\text{native scalar hierarchy integrates successfully},
}
\]

\[
\boxed{
\text{all stored variables remain finite},
}
\]

\[
\boxed{
\text{the future background and perturbation endpoint gates both pass}.
}
\]

This closes the implementation gap left by the earlier background-only
Horndeski audit.  The remaining theory question is no longer whether a healthy
future structural-clock branch can be propagated.  It is which structural
closure, if any, uniquely selects the endpoint ratio

\[
\frac{\chi_\infty}{F_\infty}
\]

and therefore fixes one value inside the currently allowed

\[
2.860445
\lesssim
{\cal N}_\infty
\le
2.894405
\]

interval.


## 46. Why the mature action does not yet select a unique 2.86--2.89 endpoint

The remaining endpoint ambiguity can now be isolated mathematically.

The mature inverse-square normal form is

\[
G_2
=
\frac1{\sigma^2}
\left(
\kappa_1Z+\kappa_2Z^2-U_0
\right),
\qquad
F\to F_\infty,
\qquad
g\to0.
\]

At the coasting fixed point define

\[
r
\equiv
\frac{\kappa_2Z_\star}{F_\infty}.
\]

The exact fixed-point relations are

\[
\boxed{
\kappa_1
=
2F_\infty(1-r),
}
\]

\[
\boxed{
\kappa_2
=
\frac{rF_\infty}{Z_\star},
}
\]

and

\[
\boxed{
U_0
=
F_\infty Z_\star(4-r).
}
\]

The perturbative health of the fixed point depends on \(r\),

\[
D_\infty=2(1+2r),
\]

\[
c_{s,\infty}^2=\frac1{1+2r},
\]

but not on the absolute value of \(Z_\star\).

This exposes an exact normalization freedom.  For any positive constant
\(\lambda\), the transformation

\[
\boxed{
Z_\star\rightarrow\lambda Z_\star,
}
\]

together with

\[
\boxed{
\kappa_2\rightarrow\frac{\kappa_2}{\lambda},
}
\]

\[
\boxed{
U_0\rightarrow\lambda U_0,
}
\]

and unchanged

\[
F_\infty,\qquad r,\qquad\kappa_1
\]

leaves

\[
w_\sigma=-\frac13,
\qquad
D_\infty,
\qquad
c_{s,\infty}^2
\]

unchanged.

But the structural lapse is

\[
{\cal N}_\infty
=
t_PS_0\sqrt{2Z_\star},
\]

so the same transformation gives

\[
\boxed{
{\cal N}_\infty
\rightarrow
\sqrt{\lambda}\,{\cal N}_\infty.
}
\]

Therefore coasting, No-Slip, ghost freedom, gradient stability, and the mature
sound speed determine the **shape** of the fixed point but do not determine its
absolute structural-clock normalization.

This is why both

\[
{\cal N}_\infty=2.86044513
\]

and

\[
{\cal N}_\infty=2.89440502
\]

can belong to healthy mature actions.

It is useful to isolate the normalization variable

\[
\boxed{
\Gamma_\infty
\equiv
\frac{\chi_\infty}{F_\infty}.
}
\]

The exact Planck identity gives

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}\,
\sqrt{\Gamma_\infty}.
}
\]

Thus the unresolved endpoint question is now one-dimensional:

\[
\boxed{
\text{what fixes }\Gamma_\infty?
}
\]

The mature normal form itself cannot answer this because \(Z_\star\) is a
normalization zero mode of the fixed-point family.

Two closures already studied correspond to two distinct normalization
statements:

\[
\chi_\infty=1
\]

gives

\[
\Gamma_\infty=\frac1{F_\infty}
\]

and therefore

\[
{\cal N}_\infty=2.86044513,
\]

whereas

\[
\chi_\infty=F_\infty
\]

gives

\[
\Gamma_\infty=1
\]

and therefore

\[
{\cal N}_\infty=\sqrt{\frac{8\pi}{3}}
=2.89440502.
\]

Neither follows from fixed-point stability alone.

A value strictly between them would correspond to

\[
1<\chi_\infty<F_\infty
\]

or, equivalently,

\[
\frac1{F_\infty}
<
\Gamma_\infty
<
1.
\]

The remembered \(2.86\)--\(2.89\) band is therefore not a numerical accident.
It is the image, under the exact Planck relation, of the remaining structural
normalization freedom.

The next derivation should consequently not perform another stability scan in
\(r\).  That direction cannot remove the degeneracy.  The missing condition
must come from physics that fixes the absolute normalization of the structural
source relative to the effective Planck mass: for example the microscopic
activation definition, a matching condition inherited from the
tracker-to-structural handoff, or another independently stated SDMC identity.

Until such a condition is derived, the scientifically correct statement is

\[
\boxed{
2.86044513
\le
{\cal N}_\infty
\le
2.89440502
}
\]

for the two currently motivated endpoint normalizations, with no unique
preferred value supplied by the mature covariant action alone.


## 47. Covariant extension of the Unified Balanced Identity

The remaining normalization zero mode can be tested against the oldest
structural identity in the manuscript.

The Part-I Unified Balanced Identity contains the bare inverse-square density

\[
\boxed{
\rho_v=\frac{\rho_P}{S^2}
}
\]

and the corresponding structural radius

\[
\boxed{
R=\ell_P S
=
\frac{c}{\sqrt{G\rho_v}}.
}
\]

The modern covariant action contains instead the positive active structural
source

\[
\boxed{
\rho_X=\chi\,\rho_v
}
\]

and a non-minimal Planck-mass factor \(F\).  In the mature No-Slip limit,

\[
F\rightarrow F_\infty=\text{constant},
\qquad
g\rightarrow0,
\]

the background Friedmann weighting and the No-Slip effective Newton sector
both carry the factor

\[
\boxed{
\frac{G_{\rm eff}}{G}=\frac1F.
}
\]

The active source therefore gravitates relative to the original Balanced
Identity source with the dimensionless weight

\[
\boxed{
\Gamma
\equiv
\frac{G_{\rm eff}\rho_X}
{G\rho_v}
=
\frac{\chi}{F}.
}
\]

This is exactly the same normalization variable isolated in Section 46.

Define the gravitational radius associated with the active covariant source by

\[
R_X
\equiv
\frac{c}
{\sqrt{G_{\rm eff}\rho_X}}.
\]

Using the original structural radius gives

\[
\boxed{
\frac{R_X}{R}
=
\sqrt{\frac{F}{\chi}}
=
\Gamma^{-1/2}.
}
\]

The mature structural coefficient is correspondingly

\[
\boxed{
\frac{HR}{c}
=
\sqrt{\frac{8\pi}{3}}
\sqrt{\Gamma}.
}
\]

Since \(p\to1\) on the mature coasting attractor,

\[
\frac{HR}{c}
\to
{\cal N}_\infty.
\]

This produces a new candidate closure condition.

If the covariant completion is required to preserve the Part-I statement that
the radius column and density column describe the **same structural state
after gravitational weighting is included**, then one should require

\[
\boxed{
R_X=R.
}
\]

Equivalently,

\[
\boxed{
G_{\rm eff}\rho_X
=
G\rho_v.
}
\]

This immediately gives

\[
\boxed{
\Gamma_\infty=1,
}
\]

hence

\[
\boxed{
\chi_\infty=F_\infty
}
\]

and therefore

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
=
2.894405018\ldots .
}
\]

This condition may be called the **covariant Balanced-Identity matching
condition**.

It should be distinguished from the alternative source-amplitude condition

\[
\boxed{
\rho_X\to\rho_v,
}
\]

which gives

\[
\chi_\infty=1
\]

and hence

\[
\boxed{
{\cal N}_\infty
=
\frac{\sqrt{8\pi/3}}{\sqrt{F_\infty}}
=
2.86044513.
}
\]

The difference between the two closures is now transparent.

The lower endpoint restores the bare **source amplitude**,

\[
\rho_X=\rho_v,
\]

but because \(F_\infty>1\), its gravitational radius is

\[
\boxed{
\frac{R_X}{R}
=
\sqrt{F_\infty}
\simeq1.01187224.
}
\]

Thus the active gravitational radius remains about

\[
\boxed{
1.1872\%
}
\]

larger than the original Balanced-Identity radius.

The upper endpoint instead restores the bare **gravitational strength**,

\[
\frac{G}{F_\infty}\rho_X
=
G\rho_v,
\]

so the two radii coincide exactly.

The released future trajectories show the distinction directly.  At
\(\ln a=10\), the C3 No-Slip-refined \(F_\infty\)-matched branch has

\[
\frac{\chi}{F}
\simeq0.99996176
\]

and therefore

\[
\frac{R_X}{R}
\simeq1.00001912,
\]

already within about \(1.9\times10^{-5}\) of the covariant Balanced-Identity
condition.

The bounded-\(\chi\) branch instead has

\[
\frac{\chi}{F}
\simeq0.97663417,
\]

with

\[
\frac{R_X}{R}
\simeq1.01189172,
\]

and tends to the non-unit asymptotic value
\(\sqrt{F_\infty}\).

This result does **not** prove that \(\chi_\infty=F_\infty\).  It identifies the
precise additional principle that selects the upper endpoint: the original
Unified Balanced Identity must remain a single density-radius state after the
covariant theory's effective gravitational coupling is included.

If that principle is accepted as part of the covariant completion, the
normalization zero mode is removed and the mature lapse is uniquely

\[
\boxed{
{\cal N}_\infty=2.894405018\ldots .
}
\]

If only the active source amplitude is required to return to the bare
inverse-square density, the lower endpoint remains the appropriate closure.

## 48. Why the Planck-to-radiation handoff does not yet fix the endpoint

The early manuscript supplies a second possible source of normalization:
the Planck boundary.

The structural boundary is

\[
S_P=1,
\qquad
t_{s,P}=t_P,
\qquad
{\cal N}_P=1,
\qquad
\rho_{v,P}=\rho_P.
\]

However, Part III already showed that this boundary cannot be identified
directly with the classical radiation tracker.

The radiation tracker has

\[
p_r=2,
\]

and, for the accepted early exponential slope,

\[
{\cal N}_r
=
\Xi_v\lambda_e.
\]

For

\[
\lambda_e=20,
\]

this is approximately

\[
{\cal N}_r\simeq57.89,
\]

which is incompatible with imposing the classical tracker directly at

\[
{\cal N}_P=1.
\]

The manuscript therefore inserts a separate Planck-to-radiation transition
with structural-to-radiation energy transfer \(J\),

\[
\dot\rho_\phi
+
3H(1+w_\phi)\rho_\phi
=
-J,
\]

\[
\dot\rho_r+4H\rho_r=+J.
\]

The structural density equation becomes

\[
2p
=
3(1+w_\phi)
+
\frac{J}{H\rho_\phi},
\]

and the radiation fraction evolves according to

\[
\frac{d\ln(\rho_r/\rho_\phi)}{d\ln a}
=
2p-4+\frac{J}{H\rho_r}.
\]

Consequently, the Planck boundary fixes the initial structural normalization,
but the accumulated transfer history

\[
\int J\,a^4\,dt
\]

controls how the classical thermal branch is reached.

The exact structural reconstruction

\[
S(t_c)
=
1+
\frac1{t_P}
\int_{t_P}^{t_c}
{\cal N}(t')\,dt'
\]

must simultaneously satisfy the density column of the Balanced Identity.
Without a microscopic law for \(J\), or an action valid continuously through
the nonclassical Planck-to-thermal transition, these equations do not supply a
unique map from

\[
{\cal N}_P=1
\]

to the late normalization \(\Gamma_\infty\).

Therefore the existing Planck-to-radiation handoff does **not** presently
select either

\[
\chi_\infty=1
\]

or

\[
\chi_\infty=F_\infty.
\]

It constrains the topology of the history -- the classical radiation tracker
must emerge only after a separate transition -- but it leaves an integrated
normalization freedom because \(J(t)\) is not yet derived.

This is useful because it isolates the remaining theoretical fork very
cleanly:

\[
\boxed{
\text{Planck boundary alone}
\;\not\Rightarrow\;
\Gamma_\infty.
}
\]

The strongest currently available candidate for removing the mature
normalization zero mode is instead the covariant Balanced-Identity condition

\[
\boxed{
G_{\rm eff}\rho_X=G\rho_v,
}
\]

which selects

\[
\boxed{
\Gamma_\infty=1
}
\]

and the exact upper endpoint

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
=
2.894405018\ldots .
}
\]

The next theoretical task is to determine whether that matching condition can
itself be derived from a variational or microscopic statement, rather than
adopted as an additional structural closure.


## 49. Renormalized activation flow and the endpoint sum rule

The covariant Balanced-Identity analysis shows that the physically relevant
normalization variable is not \(\chi\) or \(F\) separately, but

\[
\boxed{
\Gamma
\equiv
\frac{\chi}{F}.
}
\]

It is useful to promote this combination to its own structural beta function.

Define

\[
\beta_\Gamma
\equiv
\frac{d\ln\Gamma}{d\ln\sigma}.
\]

Since

\[
\beta_\chi
=
\frac{d\ln\chi}{d\ln\sigma}
=
2-m_X
\]

and

\[
\alpha_M
=
\frac{d\ln F}{d\ln a}
=
p\frac{d\ln F}{d\ln\sigma},
\]

we obtain the exact identity

\[
\boxed{
\beta_\Gamma
=
2-m_X-\frac{\alpha_M}{p}.
}
\]

Equivalently, using \(d\ln\sigma=p\,d\ln a\),

\[
\boxed{
\frac{d\ln\Gamma}{d\ln a}
=
A_\chi-\alpha_M.
}
\]

Thus the endpoint normalization is controlled by the integrated difference
between activation growth and Planck-mass running,

\[
\boxed{
\ln\frac{\Gamma_\infty}{\Gamma_0}
=
\int_{\sigma_0}^{\infty}
\left(
2-m_X-\frac{\alpha_M}{p}
\right)d\ln\sigma.
}
\]

This is the endpoint **sum rule**.

At the accepted present point,

\[
\chi_0=0.94108957,
\]

\[
F_0=1.02346657,
\]

so

\[
\boxed{
\Gamma_0
=
\frac{\chi_0}{F_0}
\simeq0.91951178.
}
\]

The covariant Balanced-Identity endpoint requires

\[
\Gamma_\infty=1,
\]

and therefore the required integrated beta-function area is

\[
\boxed{
\int\beta_\Gamma\,d\ln\sigma
=
-\ln\Gamma_0
\simeq0.08391242.
}
\]

The bare-source restoration endpoint instead has

\[
\Gamma_\infty=\frac1{F_\infty}
\simeq0.97667178,
\]

which requires

\[
\boxed{
\int\beta_\Gamma\,d\ln\sigma
\simeq0.06030779.
}
\]

The difference between the two endpoint areas is exactly

\[
\boxed{
0.08391242-0.06030779
=
0.02360463334
=
\ln F_\infty.
}
\]

This equality is not numerical coincidence.  The two closures differ only by
whether the asymptotic activation acquires the additional factor \(F_\infty\).

The released future actions satisfy the corresponding sum rules directly.

At \(\ln a=10\), the \(F_\infty\)-matched branch has accumulated approximately

\[
\ln\frac{\Gamma}{\Gamma_0}
\simeq0.08387418,
\]

already within about

\[
3.8\times10^{-5}
\]

of the exact covariant-matching area.

The bounded-\(\chi\) branch gives approximately

\[
\ln\frac{\Gamma}{\Gamma_0}
\simeq0.06026929,
\]

again within about

\[
3.9\times10^{-5}
\]

of its asymptotic target.

The endpoint problem can therefore be restated in a particularly compact
form:

\[
\boxed{
\text{What microscopic dynamics fixes }
\int\beta_\Gamma\,d\ln\sigma?
}
\]

The local fixed-point equations determine

\[
\beta_\Gamma\to0
\]

but do not determine the total accumulated area.  That area is a global
boundary-value datum.

This clarifies why neither the mature stability conditions nor the local
inverse-square normal form can choose the endpoint.  A local fixed point fixes
the **slope** at the end of the flow; the Balanced-Identity closure fixes the
**integrated normalization** of the entire flow.

If the earlier Part-III ideal endpoint

\[
{\cal N}_\infty=\Xi_v
\]

is required to survive in the covariant completion, then the sum rule must take
the \(0.08391242\) branch and therefore

\[
\boxed{
\Gamma_\infty=1,
\qquad
\chi_\infty=F_\infty.
}
\]

If instead the defining late condition is

\[
\rho_X\to\rho_v,
\]

then the \(0.06030779\) branch is selected.

The next microscopic derivation should therefore target
\(\beta_\Gamma\), rather than \(\chi\) and \(F\) independently.


## 50. Legacy Part-III endpoint continuity as a boundary condition

The covariant matching condition has an independent interpretation when the
new action is compared directly with the earlier temporal reconstruction.

Part III derived, under its flat canonical scalar closure,

\[
{\cal N}
=
\Xi_v\frac{p}{\sqrt{\Omega_\phi}},
\qquad
\Xi_v=\sqrt{\frac{8\pi}{3}},
\]

and therefore obtained the ideal late scalar-dominated endpoint

\[
\boxed{
{\cal N}_\infty=\Xi_v.
}
\]

The modern covariant action generalizes the same relation to

\[
{\cal N}
=
p\,\Xi_v
\sqrt{
\frac{\Gamma}{\widehat\Omega_X}
},
\]

where

\[
\Gamma=\frac{\chi}{F}
\]

and

\[
\widehat\Omega_X
\equiv
\frac{8\pi G\rho_X}{3FH^2}
\]

is the structural-source fraction measured relative to the actual
non-minimal Friedmann normalization.

On the mature scalar-dominated branch,

\[
p\to1,
\qquad
\widehat\Omega_X\to1,
\]

so

\[
\boxed{
{\cal N}_\infty
=
\Xi_v\sqrt{\Gamma_\infty}.
}
\]

Consequently, preserving the older Part-III endpoint in the modern covariant
completion requires

\[
\Xi_v\sqrt{\Gamma_\infty}
=
\Xi_v,
\]

hence

\[
\boxed{
\Gamma_\infty=1.
}
\]

Therefore

\[
\boxed{
\chi_\infty=F_\infty.
}
\]

This is exactly the covariant Balanced-Identity matching condition of
Section 47.

The two arguments are logically distinct but algebraically equivalent:

\[
\boxed{
R_X=R
}
\]

preserves the original density-radius column of the Unified Balanced Identity,
while

\[
\boxed{
{\cal N}_\infty=\Xi_v
}
\]

preserves the earlier temporal endpoint.  Both lead to

\[
\boxed{
\frac{\chi_\infty}{F_\infty}=1.
}
\]

The bounded-source condition

\[
\chi_\infty=1
\]

instead gives

\[
{\cal N}_\infty
=
\frac{\Xi_v}{\sqrt{F_\infty}}
=
2.86044513,
\]

so it should be interpreted as a deliberate modification of the older ideal
endpoint caused by the frozen non-minimal Planck mass.

The fractional shift is

\[
\frac{
2.86044513
}{
2.89440502
}
-1
\simeq
-1.1733\%.
\]

Thus the endpoint fork can be stated operationally:

- **legacy-preserving covariant completion:**  
  \(\Gamma_\infty=1\), \(\chi_\infty=F_\infty\),
  \({\cal N}_\infty=2.89440502\);

- **bare-source-amplitude completion:**  
  \(\chi_\infty=1\),
  \(\Gamma_\infty=1/F_\infty\),
  \({\cal N}_\infty=2.86044513\).

The present derivation therefore gives a structural reason to treat the upper
endpoint as more than an arbitrary member of the stable interval: it is the
unique member that leaves the original Part-III ideal lapse unchanged after
the covariant Planck-mass dressing is introduced.

This remains a conditional selection, because it assumes that the older ideal
endpoint is a defining SDMC boundary condition rather than a result that may be
renormalized by the later covariant completion.  The microscopic theory must
decide which interpretation is fundamental.
