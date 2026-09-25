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


## 51. Frame-invariant form of the Balanced-Identity closure

The ratio

\[
\Gamma=\frac{\chi}{F}
\]

has a further useful property: in the mature constant-\(F\) limit it is the
conformal-frame-invariant gravitational strength of the structural source.

In the Jordan-frame description,

\[
G_{\rm eff}=\frac{G}{F},
\]

so define the dimensionless structural-strength invariant

\[
\boxed{
{\cal I}_{\rm BI}
\equiv
\frac{
G_{\rm eff}\rho_X R^2
}{c^2}.
}
\]

Using

\[
\rho_X=\chi\rho_v
\]

and the original Unified Balanced Identity

\[
\frac{G\rho_vR^2}{c^2}=1,
\]

one obtains immediately

\[
\boxed{
{\cal I}_{\rm BI}
=
\frac{\chi}{F}
=
\Gamma.
}
\]

Now perform the constant conformal transformation appropriate to the mature
endpoint,

\[
g_{\mu\nu}^{E}
=
F_\infty g_{\mu\nu}^{J}.
\]

For constant \(F_\infty\),

\[
R_E=\sqrt{F_\infty}\,R_J
\]

and an energy density transforms as

\[
\rho_E=\frac{\rho_J}{F_\infty^2}.
\]

Therefore

\[
\frac{G\rho_ER_E^2}{c^2}
=
\frac{
G\rho_JR_J^2
}{
F_\infty c^2
}.
\]

Hence

\[
\boxed{
{\cal I}_{\rm BI}^{E}
=
{\cal I}_{\rm BI}^{J}
=
\Gamma_\infty.
}
\]

The endpoint variable \(\Gamma_\infty\) is therefore not an artifact of the
Jordan-frame normalization used by hi_class.

The covariant Balanced-Identity condition can be written as the frame-neutral
statement

\[
\boxed{
{\cal I}_{\rm BI}=1.
}
\]

This again gives

\[
\Gamma_\infty=1,
\qquad
\chi_\infty=F_\infty.
\]

The same condition can be stated directly as an action-level asymptotic
boundary condition.

Define the rescaled active energy

\[
{\cal R}(Z)
\equiv
\sigma^2\rho_X.
\]

Since

\[
\chi
=
\frac{S_0^2}{\rho_P}{\cal R},
\]

the covariant Balanced-Identity condition requires

\[
\boxed{
{\cal R}_\star
=
\frac{\rho_P}{S_0^2}\,F_\infty.
}
\]

By contrast, bare-source restoration requires

\[
\boxed{
{\cal R}_\star
=
\frac{\rho_P}{S_0^2}.
}
\]

Thus the endpoint distinction can be imposed at the level of the asymptotic
Lagrangian energy coefficient, before the lapse is calculated.

For the mature normal form

\[
G_2=\sigma^{-2}f(Z),
\]

\[
{\cal R}(Z)=2Zf_{,Z}-f.
\]

The coasting fixed point already satisfies

\[
f+Z_\star f_{,Z}=0,
\qquad
f_{,Z}(Z_\star)=2F_\infty.
\]

The additional Balanced-Identity equation

\[
\boxed{
\frac{S_0^2}{\rho_PF_\infty}
\left[
2Z_\star f_{,Z}(Z_\star)-f(Z_\star)
\right]
=1
}
\]

then fixes the normalization of \(Z_\star\), and therefore fixes

\[
{\cal N}_\infty
=
t_PS_0\sqrt{2Z_\star}.
\]

In this formulation the value \(2.894405\) is not inserted as a lapse target.
The boundary condition is instead placed on the rescaled energy carried by the
covariant Lagrangian, and the structural lapse follows from the scalar kinetic
solution.

This is the closest current formulation to a microscopic endpoint condition.
What is still missing is a derivation of

\[
{\cal I}_{\rm BI}=1
\]

from a deeper symmetry or variational constraint, rather than taking the
Unified Balanced Identity itself as the defining structural principle.


## 52. Kinetic form of the covariant Balanced-Identity condition

The endpoint matching condition can be expressed directly in the scalar kinetic
variable, without reference to the future lapse as an input.

At the mature inverse-square fixed point,

\[
G_2=\sigma^{-2}f(Z),
\]

and the coasting equations give

\[
f+Z_\star f_{,Z}=0,
\qquad
f_{,Z}=2F_\infty.
\]

The rescaled active energy is therefore

\[
{\cal R}_\star
=
2Z_\star f_{,Z}-f
=
3Z_\star f_{,Z}
=
6F_\infty Z_\star.
\]

Because

\[
\Gamma_\infty
=
\frac{\chi_\infty}{F_\infty}
\]

and

\[
{\cal N}_\infty
=
\Xi_v\sqrt{\Gamma_\infty},
\]

while independently

\[
{\cal N}_\infty
=
t_PS_0\sqrt{2Z_\star},
\]

define the Balanced-Identity kinetic normalization

\[
\boxed{
Z_{\rm BI}
\equiv
\frac{\Xi_v^2}
{2t_P^2S_0^2}.
}
\]

Since

\[
\Xi_v^2=\frac{8\pi}{3},
\]

this is

\[
\boxed{
Z_{\rm BI}
=
\frac{4\pi}
{3t_P^2S_0^2}.
}
\]

The mature endpoint then satisfies the exact identity

\[
\boxed{
\Gamma_\infty
=
\frac{Z_\star}{Z_{\rm BI}}.
}
\]

Thus the covariant Balanced-Identity condition

\[
\Gamma_\infty=1
\]

is simply

\[
\boxed{
Z_\star=Z_{\rm BI}.
}
\]

The effective Planck-mass factor cancels from this kinetic normalization.

This explains why the legacy-preserving endpoint remains

\[
{\cal N}_\infty=\Xi_v
\]

even when

\[
F_\infty\ne1.
\]

The \(F_\infty\) factor is compensated by the asymptotic source amplitude

\[
\chi_\infty=F_\infty,
\]

leaving the structural kinetic speed fixed at the original Balanced-Identity
value.

By contrast, the bare-source condition

\[
\chi_\infty=1
\]

gives

\[
\Gamma_\infty=\frac1{F_\infty}
\]

and therefore

\[
\boxed{
Z_\star
=
\frac{Z_{\rm BI}}{F_\infty}.
}
\]

This is the lower \(2.860445\) endpoint.

Relative to the accepted present kinetic density,

\[
\frac{Z_{\rm BI}}{Z_0}
=
\left(
\frac{\Xi_v}{{\cal N}_0}
\right)^2
\simeq0.71897971.
\]

The matched future action uses precisely this mature kinetic ratio.

The bounded-source endpoint instead has

\[
\frac{Z_\star}{Z_0}
=
\frac{0.71897971}{F_\infty}
\simeq0.70220719.
\]

The endpoint ambiguity is therefore equivalently a question about whether the
non-minimal Planck mass is allowed to renormalize the mature **structural
kinetic speed**.

The covariant Balanced-Identity branch says no:

\[
\boxed{
Z_\star=Z_{\rm BI}
}
\]

and lets the active source amplitude absorb \(F_\infty\).

The bare-source branch says yes:

\[
\boxed{
Z_\star=Z_{\rm BI}/F_\infty
}
\]

so that the active source amplitude itself returns to the bare inverse-square
normalization.

This kinetic formulation is useful for future action building because the
boundary condition can be imposed directly on \(Z_\star\), after which
\({\cal N}_\infty\) is derived from

\[
{\cal N}=t_PS_0\sqrt{2Z}
\]

rather than supplied numerically.


## 53. Residual mature shape modulus after the lapse normalization is fixed

Fixing the covariant Balanced-Identity normalization removes the freedom in
\(Z_\star\), but it does not determine the entire mature scalar Lagrangian.

For the quadratic normal form,

\[
f(Z)=\kappa_1Z+\kappa_2Z^2-U_0,
\]

the coasting conditions give

\[
\kappa_1=2F_\infty(1-r),
\]

\[
\kappa_2=\frac{rF_\infty}{Z_\star},
\]

\[
U_0=F_\infty Z_\star(4-r),
\]

where

\[
r\equiv\frac{\kappa_2Z_\star}{F_\infty}.
\]

Once

\[
Z_\star=Z_{\rm BI}
\]

is fixed by the covariant Balanced Identity, the mature family therefore still
contains one dimensionless **shape modulus**, \(r\).

This modulus controls perturbations rather than the structural-clock
normalization:

\[
\boxed{
c_{s,\infty}^2
=
\frac1{1+2r},
}
\]

\[
\boxed{
D_\infty
=
2(1+2r),
}
\]

with the exact invariant

\[
D_\infty c_{s,\infty}^2=2.
\]

The structural lapse, by contrast, is already fixed by \(Z_{\rm BI}\),

\[
{\cal N}_\infty=\Xi_v,
\]

independently of \(r\).

The simplest healthy region is

\[
r>0,
\]

which gives

\[
0<c_s^2<1.
\]

If the mature potential coefficient is also required to remain positive,

\[
U_0>0,
\]

then

\[
r<4.
\]

The reconstructed late action has already entered the negative-\(k_1\) sign
sector.  Requiring the mature coefficient to retain that sign gives

\[
\kappa_1<0
\quad\Longrightarrow\quad
r>1.
\]

Together these conditions define the particularly natural interval

\[
\boxed{
1<r<4.
}
\]

Equivalently,

\[
\boxed{
\frac19<c_{s,\infty}^2<\frac13.
}
\]

The continuity reference obtained by setting the mature sound speed equal to
the accepted present value,

\[
c_{s,0}^2\simeq0.138336,
\]

corresponds to

\[
r\simeq3.11441.
\]

The freely released \(F_\infty\)-matched future action that minimized the
finite-transition residuals instead settled on

\[
r\simeq3.75612,
\]

giving

\[
c_{s,\infty}^2\simeq0.117478.
\]

Both values lie inside the same healthy mature family.

This separates two questions that had previously been mixed together:

1. the **structural normalization problem**, which determines
   \(Z_\star\) and \({\cal N}_\infty\);

2. the **kinetic-shape problem**, which determines \(r\), \(D_\infty\), and
   \(c_{s,\infty}^2\).

The covariant Balanced Identity can resolve the first without resolving the
second.

Therefore even if

\[
{\cal N}_\infty=2.894405018\ldots
\]

is adopted as the structural endpoint, a deeper microscopic Lagrangian is
still needed to predict the mature scalar sound speed rather than selecting it
through a future-tail design criterion.


## 54. Exact relation between the old \(\lambda_\ell=\sqrt2\) endpoint and the new shape modulus

Part VI stated that the final one-scalar potential should approach the
canonical exponential slope

\[
\boxed{
\lambda_\ell=\sqrt2
}
\]

so that the asymptotic SDMC branch remains coasting.

The mature structural-field normal form now allows this earlier statement to
be translated exactly.

Take

\[
G_2
=
\frac1{\sigma^2}
\left(
\kappa_1Z+\kappa_2Z^2-U_0
\right),
\]

with

\[
\kappa_1=2F_\infty(1-r),
\]

\[
\kappa_2=\frac{rF_\infty}{Z_\star},
\]

\[
U_0=F_\infty Z_\star(4-r).
\]

The strictly canonical limit is

\[
\boxed{
r=0.
}
\]

Then

\[
\kappa_2=0,
\qquad
\kappa_1=2F_\infty,
\qquad
U_0=4F_\infty Z_\star,
\]

and the scalar Lagrangian becomes

\[
G_2
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_\star
\right).
\]

Define the canonically normalized field

\[
\boxed{
\phi
=
\sqrt{2F_\infty}\,\ln\sigma.
}
\]

Because

\[
Z=\frac12(\partial\sigma)^2,
\]

the kinetic term transforms as

\[
\frac{2F_\infty Z}{\sigma^2}
=
\frac12(\partial\phi)^2.
\]

The potential becomes

\[
V(\phi)
=
4F_\infty Z_\star
\exp\left(
-\frac{2\phi}{\sqrt{2F_\infty}}
\right).
\]

Using the effective mature Planck scale

\[
M_{\rm eff}^2=F_\infty,
\]

define

\[
\widehat\phi
\equiv
\frac{\phi}{\sqrt{F_\infty}}
=
\sqrt2\,\ln\sigma.
\]

Then

\[
\boxed{
V(\widehat\phi)
\propto
e^{-\sqrt2\,\widehat\phi}.
}
\]

Therefore

\[
\boxed{
r=0
\quad\Longleftrightarrow\quad
\lambda_\ell=\sqrt2
}
\]

for the mature scalar-dominated structural branch.

This is an exact equivalence, not a numerical analogy.

The same limit gives

\[
\boxed{
c_{s,\infty}^2=1,
}
\]

and

\[
\boxed{
D_\infty=2.
}
\]

Thus the earlier canonical late endpoint is one specific member of the modern
inverse-square fixed-point family.

The later reconstructed actions with

\[
r\simeq3.1\text{--}3.8
\]

are genuinely noncanonical mature completions.  They retain the same
background coasting law

\[
w_\sigma=-\frac13,
\qquad
p=1,
\]

but predict a reduced scalar sound speed.

This reveals a second independent legacy question.

The covariant Balanced Identity can select the mature **normalization**

\[
Z_\star=Z_{\rm BI},
\qquad
{\cal N}_\infty=\Xi_v,
\]

while the older Part-VI statement

\[
\lambda_\ell=\sqrt2
\]

selects the mature **shape**

\[
\boxed{
r=0.
}
\]

If both legacy conditions are retained simultaneously, the mature action is no
longer a one-parameter family.  Its quadratic fixed-point coefficients become

\[
\boxed{
\kappa_1=2F_\infty,
}
\]

\[
\boxed{
\kappa_2=0,
}
\]

\[
\boxed{
U_0=4F_\infty Z_{\rm BI}.
}
\]

The resulting endpoint is

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}},
\qquad
c_{s,\infty}^2=1.
}
\]

However, this mathematical identification does not establish that the
accepted late266 Horndeski action can reach the canonical branch smoothly.
The present structural reconstruction has negative \(k_1\) and a substantial
quadratic kinetic sector, so reaching \(r=0\) requires the future action to
remove the \(Z^2\) term and reverse the sign of the linear kinetic coefficient
without crossing a ghost or gradient instability.

A dedicated future-only audit has therefore been added.  It keeps

\[
\chi_\infty=F_\infty,
\qquad
{\cal N}_\infty=\Xi_v,
\]

sets

\[
r=0,
\]

and scans only the unobserved future coefficient-transition rate.  The
question is deliberately narrow:

\[
\boxed{
\text{Can the accepted covariant action reach the old canonical }
\lambda_\ell=\sqrt2
\text{ endpoint while remaining healthy?}
}
\]

If yes, the earlier Part-VI endpoint supplies a natural shape-selection
principle in addition to the Balanced-Identity normalization.

If no, the later noncanonical reconstruction has genuinely superseded the
earlier canonical late-shape assumption, even though both share the same
coasting background.


## 55. The shape modulus separates canonical and kinetic-condensate mature branches

The remaining modulus \(r\) has a direct microscopic interpretation.

For

\[
f(Z)=\kappa_1Z+\kappa_2Z^2-U_0,
\]

the kinetic derivative is

\[
f_{,Z}
=
\kappa_1+2\kappa_2Z.
\]

Using the mature coefficients,

\[
\kappa_1=2F_\infty(1-r),
\qquad
\kappa_2=\frac{rF_\infty}{Z_\star},
\]

gives

\[
\boxed{
f_{,Z}
=
2F_\infty
\left[
1-r+r\frac{Z}{Z_\star}
\right].
}
\]

At the fixed point,

\[
f_{,Z}(Z_\star)=2F_\infty>0
\]

for every value of \(r\), which is why the background coasting solution itself
does not distinguish the shapes.

Away from the fixed point the distinction is sharp.

For

\[
0\le r<1,
\]

the low-kinetic coefficient is positive,

\[
\kappa_1>0,
\]

and

\[
f_{,Z}>0
\]

for every \(Z\ge0\).

The endpoint

\[
r=0
\]

is the strictly canonical exponential branch discussed in Section 54.

For

\[
r>1,
\]

the linear coefficient becomes negative,

\[
\kappa_1<0.
\]

The positive quadratic term then restores a healthy kinetic slope only above

\[
f_{,Z}=0
\quad\Longrightarrow\quad
\boxed{
\frac{Z}{Z_\star}
=
1-\frac1r.
}
\]

Thus the \(r>1\) family has a kinetic-condensate-like structure: a wrong-sign
linear term is stabilized at finite \(Z\) by the \(Z^2\) operator.

The same threshold appeared independently in the nonlinear fixed-point basin
audit as the gradient-stability floor,

\[
\boxed{
\left(
\frac{Z}{Z_\star}
\right)_{\rm grad}
=
1-\frac1r.
}
\]

For the released \(F_\infty\)-matched tail,

\[
r=3.75611884,
\]

so

\[
\boxed{
\left(
\frac{Z}{Z_\star}
\right)_{\rm grad}
\simeq0.73377.
}
\]

The healthy future orbit approaches

\[
Z/Z_\star\to1
\]

from safely above that threshold.

The homogeneous kinetic coefficient has a weaker zero,

\[
{\cal K}
=
f_{,Z}+2Zf_{,ZZ},
\]

which gives

\[
\boxed{
\left(
\frac{Z}{Z_\star}
\right)_{\rm hom}
=
\frac{r-1}{3r}.
}
\]

For the same \(r=3.75612\),

\[
\left(
\frac{Z}{Z_\star}
\right)_{\rm hom}
\simeq0.2446.
\]

Hence gradient stability, rather than the homogeneous kinetic determinant, is
the tighter basin boundary.

The fixed-point energy decomposition is also informative.  The total mature
energy is

\[
\rho_\star=6F_\infty Z_\star.
\]

Its three quadratic-normal-form pieces contribute

\[
\boxed{
\frac{\rho_{\kappa_1}}{\rho_\star}
=
\frac{1-r}{3},
}
\]

\[
\boxed{
\frac{\rho_{\kappa_2}}{\rho_\star}
=
\frac{r}{2},
}
\]

and

\[
\boxed{
\frac{\rho_{U_0}}{\rho_\star}
=
\frac{4-r}{6}.
}
\]

For the canonical branch \(r=0\),

\[
\frac{\rho_{\kappa_1}}{\rho_\star}
=
\frac13,
\qquad
\frac{\rho_{U_0}}{\rho_\star}
=
\frac23,
\qquad
\rho_{\kappa_2}=0,
\]

which is exactly the familiar canonical coasting balance.

For \(r>1\), the linear kinetic contribution is negative while the positive
quadratic contribution overcompensates it.  The total energy and physical
fixed-point kinetic slope remain positive.

Therefore the unresolved shape question is more physical than a choice of
sound speed:

\[
\boxed{
r=0
}
\]

means a canonical exponential mature scalar, whereas

\[
\boxed{
r>1
}
\]

means a finite-kinetic condensate-like mature scalar.

The accepted late266 reconstruction presently has a negative structural
\(k_1\), so its immediate future continuation naturally lies on the
kinetic-condensate side.  The earlier Part-VI \(\lambda_\ell=\sqrt2\) endpoint
instead asks whether the theory eventually crosses back to the canonical
side after the late No-Slip release is complete.

That question is now being tested directly by the legacy-canonical future
audit rather than decided by notation.


## 56. Analytic stability eigenvalues of the canonical mature endpoint

The \(r=0\) endpoint admits a second, independent stability check because it is
exactly equivalent to a canonical exponential scalar in the constant
\(F_\infty\) limit.

Using the effective mature Planck scale and the canonically normalized field,
the potential is

\[
V(\widehat\phi)
=
V_\star e^{-\lambda_\ell\widehat\phi},
\qquad
\lambda_\ell=\sqrt2.
\]

Introduce the standard autonomous variables

\[
x
\equiv
\frac{\dot{\widehat\phi}}{\sqrt6 H},
\qquad
y
\equiv
\frac{\sqrt{V}}{\sqrt3 H},
\]

so that

\[
\Omega_\phi=x^2+y^2.
\]

The scalar-dominated exponential fixed point is

\[
x_\star=\frac{\lambda_\ell}{\sqrt6},
\qquad
y_\star=
\sqrt{1-\frac{\lambda_\ell^2}{6}}.
\]

For

\[
\lambda_\ell^2=2,
\]

this gives

\[
\boxed{
x_\star=\frac1{\sqrt3},
\qquad
y_\star=\sqrt{\frac23}.
}
\]

Hence

\[
w_\phi=x_\star^2-y_\star^2
=
-\frac13,
\]

and therefore

\[
\boxed{
q_\star=0.
}
\]

This reproduces the structural coasting fixed point without using the
inverse-square \(Z\)-space calculation.

For a residual barotropic component with

\[
p_b=(\gamma-1)\rho_b,
\]

the two linear eigenvalues of the scalar-dominated exponential point are

\[
\boxed{
\mu_1=\frac{\lambda_\ell^2-6}{2},
\qquad
\mu_2=\lambda_\ell^2-3\gamma.
}
\]

At the SDMC canonical endpoint,

\[
\boxed{
\mu_1=-2.
}
\]

For pressureless matter,

\[
\gamma=1,
\]

so

\[
\boxed{
\mu_2=-1.
}
\]

For radiation,

\[
\gamma=\frac43,
\]

so

\[
\boxed{
\mu_2=-2.
}
\]

Thus the canonical mature point is linearly attractive against both its
intrinsic scalar mode and the residual matter/radiation directions.

The eigenvalues also explain two asymptotic scalings found independently in
the structural-action calculation.

The dust eigenmode

\[
\mu_2=-1
\]

implies

\[
\boxed{
m\propto a^{-1},
}
\]

which is exactly the residual-matter loading used in the matter-forced
expansion.

The intrinsic scalar eigenmode

\[
\mu_1=-2
\]

gives

\[
\boxed{
\delta_{\rm scalar}\propto a^{-2}.
}
\]

This is precisely the order that remains after the leading
lapse--mapper cancellation,

\[
\frac{{\cal N}K_m}
{{\cal N}_\infty K_{m,\infty}}
=
1+O(a^{-2}).
\]

So the previously empirical hierarchy

\[
O(a^{-1})
\quad\text{matter correction},
\qquad
O(a^{-2})
\quad\text{intrinsic residual}
\]

is the direct eigenmode structure of the legacy canonical
\(\lambda_\ell=\sqrt2\) endpoint.

For \(r=0\),

\[
\kappa=1+2r=1,
\]

and the matter-forced formulas reduce to

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-\frac32 m+O(a^{-2}),
}
\]

\[
\boxed{
p-1
=
-\frac32 m+O(a^{-2}),
}
\]

\[
\boxed{
q
=
-\frac32 m+O(a^{-2}),
}
\]

while

\[
\boxed{
\frac{K_m}{K_{m,\infty}}-1
=
\frac{K_p}{K_{p,\infty}}-1
=
+\frac32 m+O(a^{-2}).
}
\]

The analytic canonical stability calculation therefore agrees with the
covariant structural-action asymptotics at the level of both the fixed point
and the decay exponents.

This substantially strengthens the interpretation of the \(r=0\) branch:
if the accepted late266 action can be joined to it without a transient health
violation, the mature state is not merely an imposed canonical limit; it is an
ordinary late-time attractor of the resulting canonical exponential theory.


## 57. Geometric meaning of the structural lapse and direct mapper reconstruction

The structural lapse has a simple exact geometric interpretation that is useful
for connecting the Lagrangian field to the earlier SDMC radius language.

The structural radius is

\[
R=\ell_P S
=
\ell_P S_0\sigma.
\]

Since

\[
{\cal N}
=
t_P S_0\dot\sigma
\]

and

\[
\frac{\ell_P}{t_P}=c,
\]

we obtain

\[
\boxed{
\dot R=c\,{\cal N}.
}
\]

Thus the structural lapse is literally the rate of change of the structural
radius measured in units of \(c\),

\[
\boxed{
{\cal N}=\frac{\dot R}{c}.
}
\]

This makes the mature endpoints especially transparent.

For the covariant Balanced-Identity branch,

\[
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
=
2.894405018\ldots,
\]

so

\[
\boxed{
\dot R_\infty
=
2.894405018\,c.
}
\]

For the alternative bare-source endpoint,

\[
{\cal N}_\infty=2.86044513,
\]

so

\[
\dot R_\infty
=
2.86044513\,c.
\]

The accepted present structural normalization gives approximately

\[
{\cal N}_0\simeq3.4135,
\]

hence the present structural-radius rate is

\[
\boxed{
\dot R_0\simeq3.4135\,c.
}
\]

This clarifies earlier values near \(3.3\)--\(3.4\): they refer to the
structural-radius speed encoded by the lapse, whereas the mature coasting
endpoint is lower, near \(2.89\).

The bridge exponent satisfies

\[
p
=
\frac{d\ln\sigma}{d\ln a}
=
\frac{\dot R/R}{H}.
\]

Therefore

\[
\boxed{
\frac{HR}{c}
=
\frac{{\cal N}}{p}.
}
\]

At the mature fixed point,

\[
p\to1,
\]

so

\[
\boxed{
\frac{H_\infty R_\infty}{c}
=
{\cal N}_\infty.
}
\]

The matter and baryon mappers are equally direct.

With the present normalization \(a_0=1\),

\[
\boxed{
\frac{K_m(a)}{K_{m0}}
=
\frac{K_p(a)}{K_{p0}}
=
\frac{\sigma(a)}{a}
=
\frac{S(a)}{S_0\,a}.
}
\]

Hence once the structural field is known from the Lagrangian, the mapper does
not require an independent evolution equation.

Equivalently,

\[
\boxed{
K_m(a)
=
K_{m0}\frac{S(a)}{S_0a},
}
\]

and

\[
\boxed{
K_p(a)
=
K_{p0}\frac{S(a)}{S_0a}.
}
\]

The exact differential relation is

\[
\boxed{
\frac{d\ln K_m}{d\ln a}
=
\frac{d\ln K_p}{d\ln a}
=
p-1.
}
\]

This also means that the mapper can be reconstructed directly from the lapse
history.  Since

\[
\dot S=\frac{{\cal N}}{t_P},
\]

we have

\[
S(t)
=
S(t_i)
+
\frac1{t_P}
\int_{t_i}^{t}{\cal N}(t')\,dt',
\]

and therefore

\[
\boxed{
\frac{K_m(t)}{K_{m0}}
=
\frac{
S(t_i)+t_P^{-1}\int_{t_i}^{t}{\cal N}(t')\,dt'
}{
S_0\,a(t)
}.
}
\]

The same formula holds for \(K_p/K_{p0}\).

The current epoch-mapper audit gives the following representative values:

\[
K_{m0}
=
2.1992619194\times10^{20},
\]

\[
K_{p0}
=
1.1743136282\times10^{20},
\]

with

\[
\frac{K_p}{K_m}
=
f_b^{1/3}
=
0.5339580601.
\]

At last scattering,

\[
z_\star=1090,
\qquad
a_\star=9.1659028414\times10^{-4},
\]

the structural reconstruction gives

\[
\frac{S_\star}{S_0}
=
2.7009609742\times10^{-5},
\]

and therefore

\[
\boxed{
\frac{K_m(z_\star)}{K_{m0}}
=
\frac{K_p(z_\star)}{K_{p0}}
=
0.02946748423.
}
\]

Numerically,

\[
\boxed{
K_m(z_\star)
=
6.480671593\times10^{18},
}
\]

and

\[
\boxed{
K_p(z_\star)
=
3.460406832\times10^{18}.
}
\]

This is the distinct last-scattering \(K_p\) value that appeared in the mapper
audit.

It should not be inserted automatically into a sound horizon that has already
been calculated entirely inside the local photon--baryon thermal domain.  The
mapper is a structural/domain conversion factor.  Whether it belongs in a
given observable depends on which domain the corresponding distance or
wavenumber is defined in.

For the released mature future branches the mapper approaches a constant.  In
the \(F_\infty\)-matched noncanonical reference tail,

\[
\boxed{
\frac{K_{m,\infty}}{K_{m0}}
=
\frac{K_{p,\infty}}{K_{p0}}
\simeq0.710728,
}
\]

while the bounded-\(\chi\) reference tail gives approximately

\[
0.712710.
\]

The approach to this constant is precisely the opposite-sign partner of the
lapse correction derived earlier,

\[
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-\frac{3m}{2\kappa}+O(a^{-2}),
\]

\[
\frac{K_m}{K_{m,\infty}}-1
=
\frac{K_p}{K_{p,\infty}}-1
=
+\frac{3m}{2\kappa}+O(a^{-2}).
\]

Consequently the mature product

\[
{\cal N}K_m
\]

and identically

\[
{\cal N}K_p
\]

lose their entire \(O(a^{-1})\) matter correction.

The structural clock, radius speed, and matter/baryon mappers are therefore
three parts of one kinematic identity rather than independent phenomenological
rules.


## 58. Closed mature canonical action selected by both legacy conditions

If the two independent legacy-preserving conditions are imposed together,

\[
\Gamma_\infty=1
\]

from the covariant Balanced Identity and

\[
r=0
\]

from the Part-VI canonical slope \(\lambda_\ell=\sqrt2\), then the mature
normal form becomes completely fixed up to the already measured constant
\(F_\infty\).

The kinetic normalization is

\[
\boxed{
Z_{\rm BI}
=
\frac{4\pi}{3t_P^2S_0^2}.
}
\]

The mature Horndeski functions are then

\[
\boxed{
G_{2,\infty}(\sigma,Z)
=
\frac{1}{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
}
\]

\[
\boxed{
G_{3,\infty}=0,
}
\]

and

\[
\boxed{
G_{4,\infty}
=
\frac{F_\infty}{2}.
}
\]

There is no remaining quadratic \(Z^2\) operator,

\[
\kappa_2=0.
\]

The exact fixed-point solution has

\[
\boxed{
Z=Z_{\rm BI},
}
\]

so

\[
\dot\sigma
=
\sqrt{2Z_{\rm BI}}
=
\frac{1}{t_PS_0}
\sqrt{\frac{8\pi}{3}}.
\]

Consequently,

\[
\boxed{
{\cal N}
=
t_PS_0\dot\sigma
=
\sqrt{\frac{8\pi}{3}}.
}
\]

The mature structural field therefore evolves linearly in chronological time,

\[
\boxed{
\sigma(t)
=
\sigma_c
+
\frac{\sqrt{8\pi/3}}{t_PS_0}
(t-t_c).
}
\]

Because the mapper has frozen,

\[
\frac{\sigma}{a}
\to
M_\infty=\text{constant},
\]

the FRW scale factor is asymptotically proportional to the structural field,

\[
a\propto\sigma.
\]

Hence

\[
\boxed{
H
=
\frac{\dot\sigma}{\sigma},
}
\]

which gives

\[
a(t)\propto t
\]

after an irrelevant origin shift, and therefore

\[
\boxed{
q=0.
}
\]

The scalar energy density is

\[
\rho_\sigma
=
2ZG_{2,Z}-G_2.
\]

For the mature canonical action,

\[
G_{2,Z}
=
\frac{2F_\infty}{\sigma^2},
\]

so at \(Z=Z_{\rm BI}\),

\[
\boxed{
\rho_{\sigma,\infty}
=
\frac{
6F_\infty Z_{\rm BI}
}{
\sigma^2
}.
}
\]

The pressure is simply

\[
p_{\sigma,\infty}
=
G_2
=
-\frac{
2F_\infty Z_{\rm BI}
}{
\sigma^2
},
\]

and therefore

\[
\boxed{
w_{\sigma,\infty}
=
-\frac13.
}
\]

The Friedmann equation closes identically,

\[
3F_\infty H^2
=
\frac{
6F_\infty Z_{\rm BI}
}{
\sigma^2
},
\]

because

\[
H^2
=
\frac{2Z_{\rm BI}}{\sigma^2}.
\]

Thus the coasting solution is an exact solution of the closed mature action,
not merely an imposed asymptotic scaling.

The canonical field representation makes the connection to Part VI explicit.
Define

\[
\boxed{
\widehat\phi
=
\sqrt2\ln\sigma.
}
\]

Then

\[
\sigma^{-2}
=
e^{-\sqrt2\widehat\phi},
\]

and the mature potential is

\[
\boxed{
V_\infty(\widehat\phi)
=
4F_\infty Z_{\rm BI}
e^{-\sqrt2\widehat\phi}.
}
\]

The exponential slope is therefore exactly

\[
\boxed{
\lambda_\ell=\sqrt2.
}
\]

The endpoint also has

\[
\boxed{
D_\infty=2,
\qquad
c_{s,\infty}^2=1,
}
\]

and the scalar/matter eigenmodes derived in Section 56 decay as

\[
a^{-2}
\]

and

\[
a^{-1},
\]

respectively.

This produces a compact candidate for the ultimate mature SDMC action:

\[
\boxed{
\left\{
G_2=
\sigma^{-2}
\left(
2F_\infty Z-4F_\infty Z_{\rm BI}
\right),
\quad
G_3=0,
\quad
G_4=\frac{F_\infty}{2}
\right\}.
}
\]

The remaining difficulty is no longer the endpoint action itself.  It is the
future interpolation problem: whether the accepted present-day Horndeski jet
can reach this canonical action while preserving the native hi_class
background target, No-Slip relation, ghost/gradient health, perturbation
finiteness, and the exact present boundary simultaneously.

That is the purpose of the current native canonical closure audit.


## 59. Exact sound-speed law of the mature quadratic family

The tiny \(c_s^2>1\) excursion found in the most stringent canonical stitch
search can be localized analytically.

For the mature constant-\(F\), \(G_3=0\) inverse-square family,

\[
G_2=\sigma^{-2}f(Z),
\]

with

\[
f(Z)
=
\kappa_1Z+\kappa_2Z^2-U_0,
\]

the scalar sound speed reduces exactly to the k-essence expression

\[
c_s^2
=
\frac{f_{,Z}}
{f_{,Z}+2Zf_{,ZZ}}.
\]

Using

\[
\kappa_1=2F_\infty(1-r),
\qquad
\kappa_2=\frac{rF_\infty}{Z_\star},
\]

and defining

\[
x\equiv\frac{Z}{Z_\star},
\]

gives

\[
f_{,Z}
=
2F_\infty
\left(
1-r+rx
\right),
\]

and

\[
f_{,Z}+2Zf_{,ZZ}
=
2F_\infty
\left(
1-r+3rx
\right).
\]

Therefore

\[
\boxed{
c_s^2(x)
=
\frac{
1-r+rx
}{
1-r+3rx
}.
}
\]

At the fixed point \(x=1\),

\[
\boxed{
c_{s,\star}^2
=
\frac1{1+2r},
}
\]

as obtained previously.

The canonical branch is special.  Setting

\[
r=0
\]

gives

\[
\boxed{
c_s^2=1
}
\]

for **every** value of \(Z\) inside the completed mature canonical normal
form, not merely at the fixed point.

Thus residual matter forcing of \(Z\),

\[
\frac{Z}{Z_\star}-1
=
-3m+O(a^{-2})
\]

for the canonical branch, cannot by itself produce a superluminal sound speed.

This identifies the origin of the tiny numerical maximum found in the strict
future stitch scan,

\[
c_{s,\max}^2
\simeq
1.00000011224.
\]

That maximum occurs around

\[
\ln a\simeq2.83,
\]

before the future coefficient release has become indistinguishable from the
pure \(r=0\) mature action.

Therefore

\[
\boxed{
c_s^2-1
}
\]

at that location is a **transition-sector residual**, involving the remaining
running of \(F\), \(G_3\), and/or the finite-rate \(G_2\) coefficient stitch.
It is not an intrinsic property of the canonical endpoint.

Once the canonical normal form is fully reached,

\[
F'=0,
\qquad
G_3=0,
\qquad
\kappa_2=0,
\]

and luminality is exact:

\[
\boxed{
c_s^2\equiv1.
}
\]

This changes the interpretation of the strict-subluminality search.

The numerical target should not be to push the asymptotic canonical sound
speed below one; doing so would contradict the exact \(r=0\) theory.

The correct requirement is instead

\[
\boxed{
c_s^2\le1
}
\]

through the finite transition, followed by

\[
c_s^2\to1^- \ \text{or}\ 1
\]

as the canonical action is reached.

The remaining \(1.1\times10^{-7}\) excess is consequently a measure of how
well the chosen future interpolation suppresses the transition residue.  It
does not invalidate the canonical endpoint itself.


## 60. Reconciliation of the recurring 2.84, 2.86, 2.89, 3.37, and 3.4135 numbers

Several numerically similar values occur in different parts of the SDMC
manuscript and should not be conflated.

### 60.1 The legacy \(3.37\)

Parts I--III used the working present structural benchmark

\[
\boxed{
\Xi_0
\equiv
\frac{H_0R_0}{c}
\simeq3.37.
}
\]

In the older coasting picture this was also associated numerically with the
global structural-radius rate,

\[
\dot R\simeq3.37\,c.
\]

The later Part-IX chronological closure reconstructed the temporal lapse
independently and obtained

\[
\boxed{
{\cal N}_0=3.41350846
}
\]

for

\[
S_0=2.72\times10^{61}
\]

and the accepted-action age.  The same calculation gives

\[
\boxed{
p_0=1.03421566.
}
\]

The exact identity derived in the present continuation is

\[
\frac{HR}{c}
=
\frac{{\cal N}}{p}.
\]

Therefore the current accepted-action quantities imply

\[
\boxed{
\left(\frac{H_0R_0}{c}\right)_{\rm current}
=
\frac{3.41350846}{1.03421566}
\simeq3.30057704.
}
\]

Thus the modern action separates three statements that were approximately
identified in the early working model:

\[
\boxed{
{\cal N}_0\simeq3.4135,
}
\]

\[
\boxed{
(H_0R_0/c)_{\rm current}\simeq3.3006,
}
\]

and the earlier benchmark

\[
\boxed{
\Xi_0^{\rm legacy}\simeq3.37.
}
\]

The difference is expected because the accepted present branch has

\[
p_0\ne1
\]

and also uses a later refitted \(H_0\).  The manuscript itself already notes
that the \(3.41351\) action-derived lapse differs from the \(3.37\) legacy
working value by about \(1.29\%\).

### 60.2 The structural \(2.894405\)

The value

\[
\boxed{
2.894405018\ldots
=
\sqrt{\frac{8\pi}{3}}
}
\]

is the old ideal vacuum-reference coefficient

\[
\Xi_v
\]

and the Part-III scalar-dominated canonical endpoint.

In the covariant continuation it is recovered as

\[
\boxed{
{\cal N}_\infty=2.894405018\ldots
}
\]

when

\[
\Gamma_\infty
=
\frac{\chi_\infty}{F_\infty}
=1.
\]

This is the endpoint selected by the covariant Balanced-Identity condition and
by continuity with the old \(\lambda_\ell=\sqrt2\) canonical endpoint.

### 60.3 The structural \(2.860445\)

The distinct lower value

\[
\boxed{
{\cal N}_\infty=2.86044513
}
\]

belongs to the alternative closure

\[
\chi_\infty=1.
\]

Since

\[
F_\infty>1,
\]

this gives

\[
\Gamma_\infty=\frac1{F_\infty}<1
\]

and lowers the mature lapse by about \(1.17\%\) relative to the
Balanced-Identity value.

Thus the genuine structural endpoint interval found in the modern action is

\[
\boxed{
2.86044513
\le
{\cal N}_\infty
\le
2.89440502.
}
\]

### 60.4 The manuscript's \(2.84892452\) is not a lapse

The remembered value near \(2.84\) is also present in Part IX, but in a
different context.

The geometry-refined likelihood milestone reports

\[
\boxed{
\Delta\chi^2_{\rm fair,Planck+DESI}
=
+2.84892452.
}
\]

That number is a fair Planck+DESI best-fit difference against the optimized
\(\Lambda\)CDM control.  It is **not** a structural lapse, sound speed, mapper,
or asymptotic coefficient.

An adjacent earlier likelihood point is

\[
\boxed{
\Delta\chi^2_{\rm fair,Planck+DESI,Q-best}
=
+2.86998458,
}
\]

which can easily be mistaken for the lower structural endpoint because both
begin with \(2.86\).  They are unrelated quantities.

The likelihood continuation subsequently crossed zero and therefore neither
\(2.8489\) nor \(2.86998\) should be imported into the Lagrangian endpoint
derivation.

### 60.5 Compact dictionary

The recurring numbers can therefore be kept separate as

\[
\boxed{
3.37
}
\quad
\text{legacy present structural benchmark},
\]

\[
\boxed{
3.41350846
}
\quad
\text{Part-IX action-derived present temporal lapse},
\]

\[
\boxed{
3.30057704
}
\quad
\text{current }H_0R_0/c={\cal N}_0/p_0,
\]

\[
\boxed{
2.894405018
}
\quad
\text{Balanced-Identity / canonical mature lapse},
\]

\[
\boxed{
2.86044513
}
\quad
\text{bare-source mature lapse},
\]

\[
\boxed{
2.84892452
}
\quad
\text{historical fair Planck+DESI }\Delta\chi^2\text{ value}.
\]

This resolves the apparent numerical overlap without discarding any of the
earlier derivations.


## 61. Long-horizon native closure of the canonical endpoint

The final precision test extends the native legacy-canonical branch from

\[
\ln a=5
\]

to

\[
\boxed{
\ln a=10,
\qquad
a=2.20264658\times10^4.
}
\]

The mature target is

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
=
2.894405018\ldots,
}
\]

\[
\boxed{
p_\infty=1,
\qquad
q_\infty=0,
\qquad
D_\infty=2,
\qquad
c_{s,\infty}^2=1.
}
\]

The native hi_class background at the \(\ln a=10\) endpoint gives

\[
\boxed{
{\cal N}
=
2.8943619267,
}
\]

\[
\boxed{
p
=
0.9999850300,
}
\]

\[
\boxed{
q
=
-1.45726\times10^{-5},
}
\]

\[
\boxed{
D
=
1.9999401171,
}
\]

and

\[
\boxed{
c_s^2
=
0.999999999853.
}
\]

The fractional lapse error is therefore

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}-1
=
-1.48879\times10^{-5}.
}
\]

The No-Slip combination at the final point is

\[
\boxed{
\alpha_B+2\alpha_M
\simeq
-1.20\times10^{-23}.
}
\]

Across the entire native future interval the health diagnostics remain

\[
\boxed{
D_{\min}
=
0.479233>0,
}
\]

\[
\boxed{
c_{s,\min}^2
=
0.250051>0.
}
\]

The maximum native sound speed is

\[
\boxed{
c_{s,\max}^2
=
1.000000000161.
}
\]

The excess above unity is therefore only

\[
\boxed{
1.61\times10^{-10},
}
\]

which is at the level expected from the finite table/interpolation
representation of an endpoint whose analytic \(r=0\) normal form has

\[
c_s^2\equiv1.
\]

### 61.1 Agreement with the exact homogeneous action

The independently integrated C3 No-Slip-refined action with the future rates

\[
\boxed{
\mu_{k_1}=50,
\qquad
\mu_{k_2}=5,
\qquad
\mu_V=120
}
\]

is a strict feasible canonical solution.

Its \(\ln a=10\) values are

\[
{\cal N}_{\rm action}
=
2.8943616600,
\]

\[
p_{\rm action}
=
0.9999850208,
\]

\[
q_{\rm action}
=
-1.49807\times10^{-5}.
\]

The native/action differences are only

\[
\boxed{
\frac{
{\cal N}_{\rm native}
-
{\cal N}_{\rm action}
}{
{\cal N}_{\rm action}
}
=
9.21\times10^{-8},
}
\]

\[
\boxed{
p_{\rm native}-p_{\rm action}
=
9.22\times10^{-9},
}
\]

and

\[
\boxed{
q_{\rm native}-q_{\rm action}
=
4.08\times10^{-7}.
}
\]

Thus the native implementation is reproducing the independently integrated
homogeneous Lagrangian trajectory, rather than merely approaching the same
qualitative endpoint.

The same refined action has

\[
\boxed{
c_{s,\max}^2
=
1.0000000000000004,
}
\]

which differs from exact unity only at double-precision floating-point level.
Its minimum sound speed is

\[
c_{s,\min}^2
=
0.0477261,
\]

and

\[
D_{\min}
=
0.386661>0.
\]

This resolves the earlier \(10^{-7}\)-level sound-speed excess found for the
nearby \((46,7,100)\) exploratory stitch.  That small excess belonged to a
different finite transition profile, not to the canonical endpoint itself.
The \((50,5,120)\) C3 No-Slip-refined continuation reaches the same \(r=0\)
endpoint without a resolved action-level superluminal excursion.

### 61.2 Native perturbations to \(\ln a=10\)

The same canonical native branch was propagated to \(\ln a=10\) for

\[
k
=
10^{-4},
\quad
10^{-3},
\quad
10^{-2}
\ {\rm Mpc}^{-1}.
\]

All stored perturbation columns remain finite.

A stronger late-time metric gate was applied from approximately

\[
\ln a=5
\]

to

\[
\ln a=10.
\]

For the longest wavelength,

\[
k=10^{-4}\ {\rm Mpc}^{-1},
\]

the final metric potentials are only about

\[
\boxed{
0.707
}
\]

of their \(\ln a\simeq5.1\) amplitudes.

For

\[
k=10^{-3}\ {\rm Mpc}^{-1},
\]

the ratio is approximately

\[
\boxed{
5.93\times10^{-3},
}
\]

and for

\[
k=10^{-2}\ {\rm Mpc}^{-1},
\]

approximately

\[
\boxed{
6.84\times10^{-3}.
}
\]

None of the three metric modes develops a late-time growth envelope.

The raw synchronous structural-scalar variable also decreases over the same
late interval.  Its final-to-start absolute-amplitude ratios are approximately

\[
\boxed{
0.774,
\quad
0.0148,
\quad
0.0134
}
\]

for the three increasing wavenumbers respectively.  These raw scalar values
remain gauge-dependent and are therefore reported as diagnostics rather than
used as a standalone physical stability definition.

The strengthened long-horizon workflow returns

\[
\boxed{
\texttt{NATIVE\_CANONICAL\_LONG\_HORIZON\_GATE\_PASS}.
}
\]

### 61.3 Closure status

At this stage the same mature solution has been reached through four
independent routes:

1. the covariant Balanced Identity fixes

\[
\Gamma_\infty=1
\]

and hence

\[
{\cal N}_\infty=\sqrt{8\pi/3};
\]

2. the legacy Part-VI condition

\[
\lambda_\ell=\sqrt2
\]

fixes the mature shape modulus

\[
r=0;
\]

3. the exact homogeneous Horndeski equations evolve the accepted present
action toward that fixed point while remaining healthy;

4. native hi_class background and scalar perturbation evolution reproduce the
same trajectory through

\[
\ln a=10.
\]

The mature canonical candidate can therefore be written as

\[
\boxed{
G_2
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
\qquad
G_3=0,
\qquad
G_4=\frac{F_\infty}{2},
}
\]

with

\[
\boxed{
Z_{\rm BI}
=
\frac{4\pi}{3t_P^2S_0^2}.
}
\]

Within the present numerical precision tests, the remaining task is no longer
to identify the mature endpoint.  It is to determine whether a deeper
microscopic or variational principle can derive the Balanced-Identity
normalization and the \(\lambda_\ell=\sqrt2\) shape condition from one common
fundamental statement, rather than retaining them as two independently
motivated structural boundary conditions.

## 62. No-go theorem for a background-only derivation of the canonical shape

The long-horizon closure makes it possible to ask a sharper question:

\[
\text{Can one single background structural identity derive both }
\Gamma_\infty=1
\text{ and }
r=0?
\]

For the mature scale-covariant normal form

\[
G_2=\frac{f(Z)}{\sigma^2},
\qquad
G_3\to0,
\qquad
F\to F_\infty,
\]

the answer is no if the additional identity depends only on the homogeneous
background energy, pressure, and geometry.

At the coasting fixed point the exact background equations give

\[
\boxed{
f(Z_\star)+Z_\star f_{,Z}(Z_\star)=0,
}
\]

and

\[
\boxed{
f_{,Z}(Z_\star)=2F_\infty.
}
\]

Hence

\[
\boxed{
f(Z_\star)=-2F_\infty Z_\star.
}
\]

The structural energy is

\[
\rho_{\sigma,\star}
=
2Z_\star f_{,Z}(Z_\star)-f(Z_\star)
=
6F_\infty Z_\star,
\]

while

\[
p_{\sigma,\star}=f(Z_\star)=-2F_\infty Z_\star,
\]

so

\[
w_{\sigma,\star}=-\frac13.
\]

The covariant Balanced Identity fixes the remaining **normalization**

\[
\boxed{
Z_\star=Z_{\rm BI},
}
\]

or equivalently

\[
\boxed{
\Gamma_\infty=1.
}
\]

But none of these equations contains

\[
f_{,ZZ}(Z_\star).
\]

The mature shape modulus is precisely second-derivative data.  For the
quadratic normal form,

\[
f(Z)=\kappa_1Z+\kappa_2Z^2-U_0,
\]

one has

\[
f_{,ZZ}=2\kappa_2,
\]

and therefore

\[
\boxed{
r
=
\frac{Z_\star f_{,ZZ}(Z_\star)}
{2F_\infty}.
}
\]

Thus two theories can have the same

\[
Z_\star,
\qquad
f(Z_\star),
\qquad
f_{,Z}(Z_\star),
\qquad
\rho_\star,
\qquad
p_\star,
\qquad
H_\star,
\qquad
{\cal N}_\star,
\]

while possessing different \(r\).

Their homogeneous mature background is identical, but their scalar response is
different.

This can be seen directly from

\[
\boxed{
c_{s,\star}^2
=
\frac{
f_{,Z}
}{
f_{,Z}+2Z_\star f_{,ZZ}
}
=
\frac1{1+2r}.
}
\]

Therefore \(r\) is an **off-shell/perturbative shape datum**, not a background
normalization datum.

This gives a useful no-go statement:

\[
\boxed{
\text{No identity involving only the mature homogeneous }
(\rho,p,H,R,S)
\text{ state can uniquely select }r.
}
\]

The Unified Balanced Identity can fix the gravitational normalization

\[
\Gamma_\infty=1,
\]

and inverse-square structural scaling with a frozen mapper fixes

\[
w_\infty=-\frac13,
\]

but a continuum of kinetic functions share those same background properties.

### 62.1 The minimal extra condition that does select \(r=0\)

The shape is fixed if one supplements the background state identity with an
asymptotic response condition.

The simplest such condition is exact mature scalar luminality,

\[
\boxed{
c_{s,\infty}^2=1.
}
\]

Using

\[
c_{s,\infty}^2=\frac1{1+2r}
\]

gives immediately

\[
\boxed{
r=0.
}
\]

Equivalently,

\[
\boxed{
f_{,ZZ}(Z_{\rm BI})=0.
}
\]

For the quadratic family this means

\[
\boxed{
\kappa_2=0.
}
\]

The two mature closure conditions can therefore be written in a particularly
compact way:

\[
\boxed{
{\cal I}_{\rm BI}
=
\frac{G_{\rm eff}\rho_XR^2}{c^2}
=1,
}
\]

and

\[
\boxed{
c_s^2=1.
}
\]

The first fixes the on-shell normalization,

\[
Z_\star=Z_{\rm BI},
\]

while the second fixes the local kinetic curvature,

\[
f_{,ZZ}(Z_\star)=0.
\]

Together they uniquely recover

\[
\boxed{
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z-4F_\infty Z_{\rm BI}
\right).
}
\]

### 62.2 A possible unified interpretation, but not yet a derivation

These two equations may be viewed as two aspects of a single proposed
**mature structural-vacuum recovery principle**:

1. the active source has the same gravitational strength-radius relation as the
original Unified Balanced Identity;

2. small scalar disturbances recover the undeformed local causal cone of the
canonical mature vacuum.

Under that interpretation,

\[
{\cal I}_{\rm BI}=1
\]

fixes the background state and

\[
c_s^2=1
\]

fixes the perturbative response around that state.

This is conceptually economical, but it must not be confused with a
first-principles derivation.  The current manuscript establishes the
consequences of these conditions and verifies that the accepted late action can
reach them smoothly.  It does not yet derive both from one deeper symmetry or
microscopic action.

The result is nevertheless restrictive: any deeper SDMC principle that claims
to produce the mature canonical endpoint must constrain **both**

\[
f(Z_\star)
\ \text{and}\ 
f_{,Z}(Z_\star)
\]

through the background equations and

\[
f_{,ZZ}(Z_\star)
\]

through the perturbative sector.

A purely background variational condition is insufficient.

## 63. Resonant second-order approach to the canonical fixed point

The first-order matter-forced result,

\[
{\cal N}/{\cal N}_\infty-1
=
p-1
=
q
=
-\frac32 m+O(a^{-2}),
\]

can be sharpened analytically for the \(r=0\) canonical endpoint.

Use the canonical exponential variables

\[
x
=
\frac{\dot{\widehat\phi}}{\sqrt6 H},
\qquad
y
=
\frac{\sqrt V}{\sqrt3 H},
\]

with

\[
\lambda_\ell=\sqrt2
\]

and pressureless matter.

The autonomous equations are

\[
x'
=
-3x+\sqrt3\,y^2
+\frac32x(1+x^2-y^2),
\]

\[
y'
=
-\sqrt3\,xy
+\frac32y(1+x^2-y^2),
\]

where a prime denotes \(d/d\ln a\).

The mature fixed point is

\[
x_\star=\frac1{\sqrt3},
\qquad
y_\star=\sqrt{\frac23}.
\]

Its eigenvalues are

\[
-1,\qquad -2.
\]

The \(-1\) mode is the residual matter mode and the \(-2\) mode is the
intrinsic scalar mode.

Let

\[
u\equiv a^{-1}=e^{-N}.
\]

Because twice the matter eigenvalue equals the scalar eigenvalue,

\[
2(-1)=-2,
\]

the nonlinear matter self-coupling is resonant with the intrinsic scalar mode.
Consequently the second-order asymptotic expansion necessarily contains
\(u^2\ln u\) terms.

Write

\[
x
=
\frac1{\sqrt3}
+
A u
+
\left(
B+\sqrt3 A^2\ln u
\right)u^2
+O(u^3\ln^2u),
\]

\[
y
=
\sqrt{\frac23}
+
\left(
D-\frac{\sqrt6}{2}A^2\ln u
\right)u^2
+O(u^3\ln^2u).
\]

The remaining second-order coefficients obey

\[
\boxed{
B+\sqrt2D
=
-\frac{5\sqrt3}{2}A^2.
}
\]

One free combination remains because the genuine \(-2\) scalar eigenmode can
be added independently.

Define a positive leading matter amplitude \(m_1\) by

\[
A=-\frac{\sqrt3}{2}m_1.
\]

Then the matter fraction becomes

\[
\boxed{
\Omega_m
=
m_1u
+
3m_1^2u^2
+
O(u^3\ln u).
}
\]

This expression is independent of the free intrinsic scalar-mode amplitude at
this order.

Using

\[
p=\sqrt3\,x,
\]

the bridge exponent is

\[
\boxed{
p
=
1
-\frac32m_1u
+
u^2
\left[
\frac94m_1^2\ln u
-\frac{45}{8}m_1^2
-\sqrt6\,D
\right]
+\cdots .
}
\]

The deceleration parameter is

\[
\boxed{
q
=
-\frac32m_1u
+
u^2
\left[
\frac92m_1^2\ln u
-\frac92m_1^2
-2\sqrt6\,D
\right]
+\cdots .
}
\]

For the canonical structural lapse,

\[
\frac{{\cal N}}{{\cal N}_\infty}
=
\frac{\sqrt2\,x}{y},
\]

which gives

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-\frac32m_1u
+
u^2
\left[
\frac{27}{8}m_1^2\ln u
-\frac{45}{8}m_1^2
-\frac{3\sqrt6}{2}D
\right]
+\cdots .
}
\]

The matter/baryon mapper follows from

\[
\frac{d\ln M}{dN}=p-1,
\]

and therefore

\[
\boxed{
\frac{M}{M_\infty}
=
1
+\frac32m_1u
+
u^2
\left[
-\frac98m_1^2\ln u
+\frac92m_1^2
+\frac{\sqrt6}{2}D
\right]
+\cdots ,
}
\]

where

\[
M
=
\frac{K_m}{K_{m0}}
=
\frac{K_p}{K_{p0}}.
\]

The entire \(O(a^{-1})\) correction cancels in the lapse-mapper product as
already found.  The first surviving term is now explicit,

\[
\boxed{
\frac{{\cal N}M}
{{\cal N}_\infty M_\infty}
=
1
+
u^2
\left[
\frac94m_1^2\ln u
-\frac{27}{8}m_1^2
-\sqrt6D
\right]
+\cdots .
}
\]

Similarly,

\[
\boxed{
\frac{{\cal N}/p}
{({\cal N}/p)_\infty}
=
1
+
u^2
\left[
\frac98m_1^2\ln u
-\frac{\sqrt6}{2}D
\right]
+\cdots .
}
\]

Since

\[
\ln u=-\ln a,
\]

the generic second-order correction is not merely \(a^{-2}\), but

\[
\boxed{
a^{-2}\ln a
}
\]

plus the independent \(a^{-2}\) scalar eigenmode.

This resonance explains why the first-order lapse-mapper cancellation is so
efficient while a smaller structured residual survives.  It also supplies a
more precise asymptotic template for the new native \(\ln a=10\) data.

The next audit fits the native/exact canonical trajectory to this resonant
template and tests whether the inferred matter-mode amplitude and intrinsic
scalar-mode coefficient become constant in the mature regime.


## 64. Mapper-memory closure fixes the leading canonical matter mode

The resonant fit of Section 63 leaves one important question:

\[
m_1
\]

appears there as the amplitude of the asymptotic matter eigenmode,

\[
\Omega_m
=
m_1 a^{-1}
+
3m_1^2a^{-2}
+\cdots.
\]

Is this coefficient an additional free late-time constant, or is it already
determined by the accepted present cosmology and the structural mapper?

For the canonical mature branch the answer is that it is **not free**.

Define the normalized matter/baryon mapper

\[
M(a)
\equiv
\frac{\sigma}{a}
=
\frac{K_m(a)}{K_{m0}}
=
\frac{K_p(a)}{K_{p0}}.
\]

Since

\[
\frac{d\ln M}{d\ln a}=p-1,
\]

and

\[
M(1)=1,
\]

the frozen future mapper is exactly

\[
\boxed{
M_\infty
=
\exp\left[
\int_0^\infty
\left(p(N)-1\right)dN
\right].
}
\]

Thus \(M_\infty\) is a memory of the complete finite transition from the
accepted present action into the mature coasting action.  It is not fixed by
the local mature normal form alone.

For the refined canonical continuation,

\[
\boxed{
M_\infty
=
0.75880292265.
}
\]

Equivalently,

\[
\boxed{
\ln M_\infty
=
-0.27601318926.
}
\]

This gives the absolute future mapper values

\[
\boxed{
K_{m,\infty}
=
1.66880637\times10^{20},
}
\]

and

\[
\boxed{
K_{p,\infty}
=
8.91072613\times10^{19}.
}
\]

These values belong to the refined canonical future history.  They should not
be confused with the earlier noncanonical reference-tail values near
\(0.71\,K_{m0}\) and \(0.71\,K_{p0}\), because the absolute mapper endpoint
retains information about the finite transition path even when the final
coasting fixed point is the same type of structural state.

### 64.1 Algebraic prediction of the resonant coefficient

On the mature canonical branch,

\[
\sigma
\sim
M_\infty a,
\]

and

\[
\dot\sigma
\to
v_\infty
=
\sqrt{2Z_{\rm BI}}.
\]

Since

\[
H
=
\frac{\dot\sigma}{\sigma}
\]

at the fixed point,

\[
H
\sim
\frac{
v_\infty
}{
M_\infty a
}.
\]

Ordinary pressureless matter remains conserved,

\[
\rho_m
=
\rho_{m0}a^{-3}.
\]

Therefore

\[
\Omega_m
=
\frac{
\rho_m
}{
3F_\infty H^2
}
\]

has the asymptotic form

\[
\Omega_m
=
\frac{
\rho_{m0}M_\infty^2
}{
3F_\infty v_\infty^2
}
a^{-1}
+\cdots.
\]

Hence

\[
\boxed{
m_1
=
\frac{
\rho_{m0}M_\infty^2
}{
3F_\infty v_\infty^2
}.
}
\]

Using

\[
v_\infty
=
v_0
\frac{{\cal N}_\infty}{{\cal N}_0},
\]

and

\[
p_0
=
\frac{v_0}{H_0},
\]

this becomes the dimensionless present-normalization formula

\[
\boxed{
m_1
=
\frac{
\Omega_{m0}M_\infty^2
}{
F_\infty p_0^2
}
\left(
\frac{{\cal N}_0}{{\cal N}_\infty}
\right)^2.
}
\]

No asymptotic fit parameter appears on the right-hand side.

For the accepted-action quantities,

\[
\Omega_{m0}
=
0.30075556084,
\]

\[
F_\infty
=
1.02388542769,
\]

\[
p_0
=
1.03421565856,
\]

\[
{\cal N}_0
=
3.41350846465,
\]

\[
{\cal N}_\infty
=
2.89440501823,
\]

and

\[
M_\infty
=
0.75880292265,
\]

the algebraic prediction is

\[
\boxed{
m_1^{\rm pred}
=
0.21992846314.
}
\]

The independent resonant fit of Section 63 gives

\[
\boxed{
m_1^{\rm fit}
=
0.21992943871.
}
\]

The absolute difference is

\[
9.76\times10^{-7},
\]

and the fractional difference is only

\[
\boxed{
4.44\times10^{-6}.
}
\]

Thus the asymptotic matter eigenmode amplitude is closed by the present matter
normalization plus the integrated structural mapper.

### 64.2 The leading approach coefficient is therefore fixed

The universal first-order canonical formulas become

\[
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-\frac32m_1a^{-1}
+\cdots,
\]

\[
p
=
1
-\frac32m_1a^{-1}
+\cdots,
\]

\[
q
=
-\frac32m_1a^{-1}
+\cdots,
\]

and

\[
\frac{M}{M_\infty}
=
1
+\frac32m_1a^{-1}
+\cdots.
\]

For the accepted refined canonical trajectory,

\[
\boxed{
\frac32m_1
=
0.32989269471.
}
\]

Hence

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-
0.32989269471\,a^{-1}
+
O(a^{-2}\ln a).
}
\]

At

\[
\ln a=10,
\]

this first-order prediction alone gives

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}
=
0.999985022895.
}
\]

The exact homogeneous canonical trajectory gives

\[
0.999985019997,
\]

so the first-order expression is already correct to about

\[
3\times10^{-9}
\]

in absolute normalized lapse at that epoch.  The remaining difference is
accounted for by the resonant \(a^{-2}\ln a\) and intrinsic \(a^{-2}\) terms
derived in Section 63.

### 64.3 Structural interpretation

The frozen mapper also determines the asymptotic relation between the
structural field and the FRW scale factor,

\[
\boxed{
\sigma
=
M_\infty a
+
\frac32m_1M_\infty
+
O(a^{-1}\ln a).
}
\]

For the refined canonical branch,

\[
\boxed{
\frac32m_1M_\infty
\simeq
0.250324.
}
\]

Thus the late structural scale is not merely proportional to \(a\); it carries
a finite additive memory of the matter-loaded transition.

Equivalently,

\[
aH
\to
\frac{
v_\infty
}{
M_\infty
}
=
\text{constant}.
\]

The mature canonical universe therefore has a constant comoving Hubble radius,

\[
\boxed{
\frac{c}{aH}
=
\text{constant},
}
\]

as expected for exact coasting.

Since

\[
R
=
\ell_PS_0\sigma,
\]

one also has

\[
\frac{R}{a}
\to
\ell_PS_0M_\infty
=
\text{constant},
\]

and therefore

\[
\boxed{
\frac{R}{c/H}
=
\frac{HR}{c}
\to
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}.
}
\]

The mature lapse is therefore simultaneously

1. the structural-clock rate \(d t_s/dt\);
2. the structural-radius speed in units of \(c\);
3. the asymptotic ratio of the structural radius to the Hubble radius.

The new mapper-memory closure shows that the **value of the fixed point** and
the **rate at which the universe approaches it** are not independent pieces of
the theory.  Once the accepted present normalization and the finite future
mapper history are specified, the leading matter-loaded approach coefficient
is already determined.


## 65. A minimal-derivative principle selects the canonical shape

Section 62 proved that the Balanced Identity and the mature homogeneous
background cannot by themselves determine the kinetic-shape modulus \(r\),
because \(r\) depends on

\[
f_{,ZZ}(Z_\star).
\]

The canonical condition was therefore written there as an additional response
condition,

\[
c_{s,\infty}^2=1.
\]

There is, however, a more structural way to state the same selection.

Assume that after the finite late activation has completely switched off, the
mature structural vacuum is described by the **leading two-derivative
scale-covariant scalar action**.

The inverse-square structural scaling already requires

\[
G_{2,\infty}
=
\frac{f(Z)}{\sigma^2}.
\]

Here

\[
Z=\frac12(\partial\sigma)^2
\]

contains two first derivatives of the structural field.

A term linear in \(Z\) is therefore part of the leading two-derivative scalar
theory.  By contrast,

\[
Z^2,\ Z^3,\ldots
\]

contain higher powers of derivatives and belong to higher-order operators in
the derivative expansion, even though the special Horndeski/k-essence
structure keeps the field equations second order.

If the mature vacuum retains only the leading derivative order, then

\[
\boxed{
f(Z)=A Z-B.
}
\]

Equivalently,

\[
\boxed{
f_{,ZZ}=0.
}
\]

For the quadratic family this is precisely

\[
\boxed{
r=0.
}
\]

Thus canonical shape selection can be interpreted as **mature EFT
minimality**, rather than as an independently inserted numerical sound-speed
condition.

### 65.1 The coefficients are then fixed uniquely

The mature coasting equations are

\[
f(Z_\star)+Z_\star f_{,Z}(Z_\star)=0,
\]

and

\[
f_{,Z}(Z_\star)=2F_\infty.
\]

For

\[
f(Z)=AZ-B,
\]

the second equation immediately gives

\[
\boxed{
A=2F_\infty.
}
\]

The first then gives

\[
AZ_\star-B+AZ_\star=0,
\]

hence

\[
\boxed{
B=2AZ_\star.
}
\]

The covariant Balanced Identity fixes

\[
Z_\star=Z_{\rm BI},
\]

so

\[
\boxed{
B=4F_\infty Z_{\rm BI}.
}
\]

Therefore the three assumptions

\[
\boxed{
\text{inverse-square scale covariance},
}
\]

\[
\boxed{
\text{leading two-derivative mature dynamics},
}
\]

and

\[
\boxed{
{\cal I}_{\rm BI}=1
}
\]

produce uniquely

\[
\boxed{
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right).
}
\]

Together with the already established mature limits

\[
G_{3,\infty}=0,
\qquad
G_{4,\infty}=\frac{F_\infty}{2},
\]

this is exactly the canonical action reached by the successful native
long-horizon continuation.

### 65.2 Luminality becomes a consequence

Because

\[
f_{,ZZ}=0,
\]

the mature k-essence sound speed is

\[
c_s^2
=
\frac{f_{,Z}}
{f_{,Z}+2Zf_{,ZZ}}
=
1.
\]

Thus

\[
\boxed{
c_{s,\infty}^2=1
}
\]

need not be treated as an independent mature boundary number if
two-derivative minimality is adopted.  It follows automatically.

Likewise the quadratic shape modulus obeys

\[
r
=
\frac{
Z_\star f_{,ZZ}(Z_\star)
}{
2F_\infty
}
=
0.
\]

The Part-VI exponential slope then follows from the canonical field
redefinition,

\[
\widehat\phi=\sqrt2\ln\sigma,
\]

giving

\[
V(\widehat\phi)
\propto
e^{-\sqrt2\widehat\phi},
\]

and therefore

\[
\boxed{
\lambda_\ell=\sqrt2.
}
\]

So the chain becomes

\[
\boxed{
\text{two-derivative mature EFT}
\Longrightarrow
r=0
\Longrightarrow
c_s^2=1
\Longrightarrow
\lambda_\ell=\sqrt2
}
\]

once inverse-square scale covariance and the coasting fixed-point equations are
imposed.

### 65.3 Scientific status of this step

This is stronger than simply choosing \(r=0\), but it is still a conditional
derivation.

The current results establish that:

1. the accepted late Horndeski action can evolve smoothly into this
two-derivative canonical endpoint;
2. the endpoint is a stable coasting attractor;
3. native background and scalar perturbations remain healthy through
\(\ln a=10\);
4. the endpoint reproduces the legacy \(\lambda_\ell=\sqrt2\) result.

What has **not** yet been proved is that the microscopic SDMC theory must
eliminate all higher-derivative structural operators in the mature vacuum.

The new statement is therefore a candidate fundamental principle:

\[
\boxed{
\textit{The fully relaxed structural vacuum is the minimal
two-derivative scale-covariant completion of the Unified Balanced Identity.}
}
\]

If this principle is accepted, the mature Lagrangian is no longer selected by
two unrelated conditions.  The Balanced Identity fixes its normalization,
while derivative minimality fixes its kinetic shape, and the coasting field
equations determine the remaining coefficient.

A deeper microscopic derivation would have to explain why the finite
late-time Horndeski operators are generated during the activation era but
decouple as the structural vacuum relaxes.


## 66. Decoupling hierarchy of the finite Horndeski operators

The successful canonical future stitch also shows explicitly how the
non-minimal operators disappear.

Write

\[
x\equiv\ln\sigma.
\]

The future tail is constructed in the rescaled \(G_2\) coefficients

\[
\bar k_1\equiv\sigma^2k_1,
\qquad
\bar k_2\equiv\sigma^2k_2,
\qquad
\bar V\equiv\sigma^2V.
\]

For the refined canonical continuation

\[
(\mu_{k_1},\mu_{k_2},\mu_V)
=
(50,5,120),
\]

the asymptotic tails have the form

\[
\bar k_1
=
2F_\infty
+
e^{-50x}P_{k_1}(x),
\]

\[
\bar k_2
=
e^{-5x}P_{k_2}(x),
\]

and

\[
\bar V
=
4F_\infty Z_{\rm BI}
+
e^{-120x}P_V(x),
\]

where the \(P(x)\) are finite matching polynomials required to preserve the
present action jet.

The non-minimal gravitational sector relaxes as

\[
F
=
F_\infty
+
e^{-\mu_Fx}P_F(x),
\]

with

\[
\boxed{
\mu_F
=
5.57322517\ldots,
}
\]

while the linear braiding coefficient obeys

\[
g
=
e^{-\mu_gx}P_g(x),
\]

with

\[
\boxed{
\mu_g
=
\mu_F+1
=
6.57322517\ldots .
}
\]

Returning to the physical \(G_2\) coefficients gives

\[
\boxed{
k_1
=
\frac{2F_\infty}{\sigma^2}
+
O\!\left(
\sigma^{-52}\,\mathrm{poly}(\ln\sigma)
\right),
}
\]

\[
\boxed{
k_2
=
O\!\left(
\sigma^{-7}\,\mathrm{poly}(\ln\sigma)
\right),
}
\]

and

\[
\boxed{
V
=
\frac{4F_\infty Z_{\rm BI}}{\sigma^2}
+
O\!\left(
\sigma^{-122}\,\mathrm{poly}(\ln\sigma)
\right).
}
\]

Therefore the relative correction from the higher-kinetic operator

\[
k_2Z^2
\]

to the leading mature \(G_2\) density decays as

\[
\boxed{
\sigma^{-5}\,\mathrm{poly}(\ln\sigma).
}
\]

The \(F\) correction decays approximately as

\[
\sigma^{-5.573},
\]

and the braiding coefficient as

\[
\sigma^{-6.573}.
\]

### 66.1 Matter becomes the slow mode

On the canonical attractor,

\[
\sigma
\sim
M_\infty a.
\]

Hence all of these finite-action transition remnants disappear at least as
fast as approximately

\[
a^{-5}
\]

relative to the mature \(G_2\) sector.

Ordinary conserved matter, however, produces the eigenmode

\[
\Omega_m
\sim
m_1a^{-1}.
\]

Its nonlinear self-coupling produces the slower resonant correction

\[
a^{-2}\ln a.
\]

Thus, once the canonical transition has completed,

\[
\boxed{
a^{-1}
\gg
a^{-2}\ln a
\gg
a^{-5}
}
\]

as \(a\to\infty\).

The leading deviation from exact coasting is therefore **not** caused by a
slowly dying Horndeski operator.  It is caused by the residual conserved
matter component.

This explains why the Section-63 resonant canonical formulas reproduce the
long-horizon numerical trajectory with extremely high accuracy even though the
full future action began from a noncanonical, braided scalar-tensor state.

### 66.2 EFT interpretation

The successful future action therefore realizes the minimal-derivative
principle of Section 65 dynamically rather than by an abrupt projection.

The finite late transition contains

\[
F_{,\sigma}\ne0,
\qquad
g\ne0,
\qquad
k_2\ne0,
\]

but the corresponding operators become irrelevant along the mature flow,

\[
F\to F_\infty,
\qquad
g\to0,
\qquad
k_2\to0.
\]

The surviving action is

\[
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
\]

\[
G_{3,\infty}=0,
\qquad
G_{4,\infty}=\frac{F_\infty}{2}.
\]

In this sense the canonical mature theory acts as an infrared fixed action for
the successful future completion.

The numerical decay exponents used in the current stitch are not claimed to
be fundamental renormalization-group critical exponents.  They are one
explicit healthy interpolation compatible with the accepted present jet.

What is robust is the hierarchy required for the interpretation:

\[
\boxed{
\text{finite activation operators must decay faster than the physical
matter eigenmode if the canonical mature attractor is to control the
late-time expansion.}
}
\]

The current continuation satisfies that requirement by a wide margin.


## 67. Robust-luminality theorem for the mature structural vacuum

Section 62 showed that exact luminality at a **single** fixed point,

\[
c_s^2(Z_\star)=1,
\]

only constrains the local curvature

\[
f_{,ZZ}(Z_\star)=0.
\]

That still leaves open the logical possibility that higher-order kinetic
operators conspire to have zero curvature at one special value of \(Z\).

The matter-loaded canonical attractor supplies a stronger condition.

Even after the finite Horndeski activation has disappeared, residual conserved
matter perturbs the mature solution away from the exact fixed point.  Therefore
the late trajectory samples a continuous interval of kinetic values,

\[
Z=Z_\star+\delta Z(a),
\qquad
\delta Z\neq0
\]

for every finite \(a\), with

\[
\delta Z\to0
\]

only asymptotically.

For the mature constant-\(F\), \(G_3=0\) structural action

\[
G_2=\frac{f(Z)}{\sigma^2},
\]

the exact scalar sound speed is

\[
\boxed{
c_s^2(Z)
=
\frac{f_{,Z}}
{f_{,Z}+2Zf_{,ZZ}}.
}
\]

Suppose the fully relaxed structural vacuum is required to recover the
undeformed scalar causal cone not only at the exact fixed point but throughout
a neighborhood sampled by the residual-matter trajectory,

\[
\boxed{
c_s^2(Z)=1
\quad
\text{for all }
Z\in I
}
\]

for some open interval \(I\) containing \(Z_\star\).

Assuming

\[
Z>0
\]

and a nondegenerate kinetic response

\[
f_{,Z}\neq0,
\]

the sound-speed equation gives

\[
f_{,Z}
=
f_{,Z}+2Zf_{,ZZ},
\]

hence

\[
\boxed{
f_{,ZZ}(Z)=0
\qquad
\forall Z\in I.
}
\]

Therefore

\[
\boxed{
f(Z)=AZ-B
}
\]

throughout that interval.

If \(f\) is analytic on the mature branch, this local result extends to the
connected mature domain.

Thus **robust mature luminality** uniquely selects the canonical kinetic shape
without separately assuming a polynomial truncation or two-derivative
minimality.

### 67.1 Why a single-point condition is weaker

For a general analytic expansion

\[
f(Z)
=
\sum_{n=0}^{\infty}c_n Z^n,
\]

the single-point condition

\[
f_{,ZZ}(Z_\star)=0
\]

is only

\[
\sum_{n\ge2}
n(n-1)c_n Z_\star^{\,n-2}
=
0.
\]

Different higher-order coefficients can cancel in this one equation.

So

\[
c_s^2(Z_\star)=1
\]

alone does **not** prove

\[
c_{n\ge2}=0.
\]

By contrast, if

\[
f_{,ZZ}(Z)=0
\]

on an open interval, analyticity requires every higher derivative contribution
to vanish there, giving

\[
\boxed{
c_{n\ge2}=0
}
\]

for the mature analytic branch.

This strengthens the shape-selection statement substantially.

### 67.2 Connection to the matter-forced orbit

For the canonical attractor,

\[
\Omega_m
=
m_1a^{-1}
+
3m_1^2a^{-2}
+\cdots,
\]

and the scalar velocity is correspondingly displaced from its exact fixed-point
value.

Therefore the residual-matter universe naturally performs the required
neighborhood test: it approaches \(Z_\star\) continuously rather than sitting
at one isolated value.

The successful native trajectory simultaneously shows

\[
c_s^2\to1,
\]

while the exact mature canonical normal form gives

\[
c_s^2\equiv1
\]

for every \(Z\) after the finite activation operators have decoupled.

Hence the canonical completion passes the stronger criterion automatically.

### 67.3 Scale covariance alone does not suppress higher powers of \(Z\)

This also exposes a useful limitation.

Consider the general scale-covariant mature action

\[
G_2
=
\frac1{\sigma^2}
\sum_{n=0}^{\infty}c_n Z^n.
\]

At the coasting fixed point,

\[
Z\to Z_\star=\text{constant}.
\]

Every term therefore scales as

\[
\frac{c_n Z_\star^n}{\sigma^2}.
\]

Consequently all powers \(Z^n\) have the **same cosmological**
\(\sigma^{-2}\) scaling on the fixed-point background.

Therefore

\[
\boxed{
\text{inverse-square scale covariance by itself does not make }
Z^2,Z^3,\ldots
\text{ irrelevant}.
}
\]

The disappearance of the higher-kinetic sector in the successful future
completion is a genuine coefficient flow,

\[
c_{n\ge2}\to0,
\]

not a trivial consequence of the background expansion.

This is the sharp form of the remaining microscopic problem.

### 67.4 Strongest current mature-vacuum closure statement

The mature endpoint can now be characterized by two physically distinct but
tightly connected requirements:

\[
\boxed{
{\cal I}_{\rm BI}
=
\frac{G_{\rm eff}\rho_XR^2}{c^2}
=1
}
\]

for the background normalization, and

\[
\boxed{
c_s^2(Z)=1
\quad
\text{through a neighborhood of }
Z_{\rm BI}
}
\]

for the local response.

The first gives

\[
Z_\star=Z_{\rm BI},
\]

while the second gives

\[
f(Z)=AZ-B.
\]

The coasting equations then fix

\[
A=2F_\infty,
\qquad
B=4F_\infty Z_{\rm BI}.
\]

Hence

\[
\boxed{
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
}
\]

with

\[
G_{3,\infty}=0,
\qquad
G_{4,\infty}=\frac{F_\infty}{2}.
\]

This is stronger than imposing \(c_s^2=1\) at one isolated endpoint and gives a
precise target for a future microscopic derivation:

\[
\boxed{
\text{Why should the relaxed SDMC vacuum recover a robust undeformed scalar
causal cone over nearby matter-loaded states?}
}
\]

Answering that question would explain dynamically why the higher-order
kinetic operators flow away rather than merely documenting that a healthy
future interpolation can make them do so.


## 68. Causal-cone shape modulus and the meaning of the \(3.1144\) value

The mature kinetic-shape modulus has a direct causal interpretation.

For the constant-\(F\), unbraided mature action

\[
G_2=\frac{f(Z)}{\sigma^2},
\]

the exact scalar sound speed is

\[
c_s^2
=
\frac{f_{,Z}}
{f_{,Z}+2Zf_{,ZZ}}.
\]

Define the dimensionless causal-shape variable

\[
\boxed{
r_{\rm cone}
\equiv
\frac12
\left(
\frac1{c_s^2}-1
\right).
}
\]

Then

\[
\boxed{
r_{\rm cone}
=
\frac{Zf_{,ZZ}}{f_{,Z}}.
}
\]

For the quadratic mature family,

\[
f(Z)
=
\kappa_1Z+\kappa_2Z^2-U_0,
\]

this becomes, at the fixed point,

\[
\boxed{
r_{\rm cone}=r.
}
\]

So the shape parameter that previously appeared algebraically in

\[
\kappa_2
=
\frac{rF_\infty}{Z_\star}
\]

is exactly the amount by which the scalar causal cone is narrowed relative to
the luminal canonical cone.

The relation may be inverted as

\[
\boxed{
c_{s,\infty}^2
=
\frac1{1+2r}.
}
\]

This immediately reproduces the two noncanonical mature reference tails.

For the \(F_\infty\)-matched tail,

\[
r=3.75611884,
\]

so

\[
c_{s,\infty}^2
=
0.1174779\ldots .
\]

For the bounded-\(\chi\) reference tail,

\[
r=3.43684444,
\]

giving

\[
c_{s,\infty}^2
=
0.127005\ldots .
\]

And for the canonical branch,

\[
\boxed{
r=0
\quad\Longleftrightarrow\quad
c_s^2=1.
}
\]

### 68.1 The present \(3.114415\) number

The accepted present action has

\[
\boxed{
c_{s,0}^2
=
0.13833496979.
}
\]

If one translates that present causal-cone width into the equivalent
k-essence curvature using the definition above, one finds

\[
\boxed{
r_{{\rm cone},0}^{\rm equiv}
=
\frac12
\left(
\frac1{0.13833496979}-1
\right)
=
3.11441507.
}
\]

This explains the recurring value near

\[
\boxed{3.1144}
\]

seen during the earlier action investigation.

It is important, however, not to overinterpret it.

At the accepted present epoch,

\[
\alpha_M\neq0,
\qquad
\alpha_B\neq0,
\]

and the scalar is still described by the full Horndeski system rather than the
constant-\(F\), \(G_3=0\) mature k-essence limit.

Therefore

\[
\boxed{
3.11441507
}
\]

is a **causal-cone-equivalent shape diagnostic** at the present epoch, not the
literal mature quadratic coefficient \(r\).

Only after

\[
F_{,\sigma}\to0,
\qquad
G_3\to0
\]

does

\[
r_{\rm cone}
\]

become identical to the action's kinetic-shape modulus.

### 68.2 A direct observable measure of shape relaxation

The relation

\[
1+2r_{\rm cone}
=
\frac1{c_s^2}
\]

gives

\[
\boxed{
\frac{d\ln(1+2r_{\rm cone})}{dN}
=
-
\frac{d\ln c_s^2}{dN}.
}
\]

Thus the entire mature shape relaxation can be tracked directly through the
scalar sound-speed history.

The canonical future completion has

\[
r_{\rm cone}\to0,
\]

while the old noncanonical tails approach finite positive constants.

This gives a useful conceptual separation:

\[
\boxed{
\Gamma
=
\frac{\chi}{F}
}
\]

measures **background normalization relaxation**, while

\[
\boxed{
r_{\rm cone}
}
\]

measures **kinetic causal-shape relaxation**.

The mature canonical SDMC endpoint is therefore the simultaneous limit

\[
\boxed{
\Gamma\to1,
\qquad
r_{\rm cone}\to0.
}
\]

The first condition restores the covariant Unified Balanced Identity,

\[
G_{\rm eff}\rho_XR^2/c^2=1,
\]

while the second restores the undeformed scalar causal cone.

### 68.3 Two-dimensional closure plane

It is therefore useful to regard the late structural theory as evolving in a
two-dimensional closure plane,

\[
\boxed{
(\Gamma,r_{\rm cone}).
}
\]

The important mature states occupy different points:

\[
\text{bounded-source noncanonical}
\quad\to\quad
\left(
\frac1{F_\infty},
\,3.43684
\right),
\]

\[
\text{Balanced-Identity noncanonical}
\quad\to\quad
\left(
1,
\,3.75612
\right),
\]

whereas the legacy-preserving canonical closure is

\[
\boxed{
(\Gamma_\infty,r_{\rm cone,\infty})
=
(1,0).
}
\]

This representation makes clear why background normalization and kinetic shape
were initially independent.

The Unified Balanced Identity fixes the horizontal coordinate,

\[
\Gamma_\infty=1,
\]

but does not determine the causal curvature coordinate.

Robust luminality fixes the vertical coordinate,

\[
r_{\rm cone,\infty}=0.
\]

Only together do they select the closed canonical mature action.

The outstanding microscopic problem can therefore be stated even more
sharply:

\[
\boxed{
\text{What underlying structural dynamics drives }
(\Gamma,r_{\rm cone})
\to
(1,0)?
}
\]

Any future microscopic completion should explain both flows, rather than only
reproducing the final background expansion.


## 69. Two-beta form of the remaining microscopic closure

The late structural state is now described by two independent coordinates,

\[
\Gamma=\frac{\chi}{F},
\qquad
r_{\rm cone}
=
\frac12\left(\frac1{c_s^2}-1\right).
\]

The closed canonical endpoint is

\[
\boxed{(\Gamma,r_{\rm cone})=(1,0)}.
\]

For the normalization sector,

\[
\boxed{
\beta_\Gamma
\equiv
\frac{d\ln\Gamma}{d\ln\sigma}
=
2-m_X-\frac{\alpha_M}{p}.
}
\]

Since \(d\ln\sigma=p\,dN\),

\[
\boxed{
\frac{d\Gamma}{dN}
=
p\Gamma\beta_\Gamma.
}
\]

For the causal-shape sector,

\[
1+2r_{\rm cone}=\frac1{c_s^2},
\]

so

\[
\boxed{
\frac{dr_{\rm cone}}{dN}
=
-\frac{1+2r_{\rm cone}}{2}
\frac{d\ln c_s^2}{dN}.
}
\]

Equivalently,

\[
\boxed{
\beta_r
\equiv
\frac{dr_{\rm cone}}{d\ln\sigma}
=
-\frac{1+2r_{\rm cone}}{2p}
\frac{d\ln c_s^2}{dN}.
}
\]

Thus the unresolved microscopic problem can be written as a two-dimensional
flow toward \((1,0)\).

Let

\[
\epsilon_\Gamma=1-\Gamma,
\qquad
\epsilon_r=r_{\rm cone}.
\]

Near the mature state, a general linear relaxation law has the form

\[
\frac{d}{dN}
\begin{pmatrix}
\epsilon_\Gamma\\
\epsilon_r
\end{pmatrix}
=
-
\begin{pmatrix}
\lambda_{\Gamma\Gamma}&\lambda_{\Gamma r}\\
\lambda_{r\Gamma}&\lambda_{rr}
\end{pmatrix}
\begin{pmatrix}
\epsilon_\Gamma\\
\epsilon_r
\end{pmatrix}
+
O(\epsilon^2).
\]

Local attraction requires the relaxation matrix to have eigenvalues with
positive real parts.  For a real two-dimensional system this is equivalent to

\[
\boxed{{\rm tr}\,{\bf M}>0},
\qquad
\boxed{\det{\bf M}>0}.
\]

If the flow is approximately diagonal,

\[
\epsilon_\Gamma'\simeq-\lambda_\Gamma\epsilon_\Gamma,
\qquad
\epsilon_r'\simeq-\lambda_r\epsilon_r,
\]

with positive rates, then

\[
{\cal L}
=
\frac12
\left(
\epsilon_\Gamma^2+\eta\epsilon_r^2
\right),
\qquad \eta>0,
\]

decreases locally as

\[
\frac{d{\cal L}}{dN}
=
-\lambda_\Gamma\epsilon_\Gamma^2
-\eta\lambda_r\epsilon_r^2
+O(\epsilon^3).
\]

This gives a precise target for a microscopic completion: derive the beta
functions that make the Balanced-Identity normalization and canonical causal
shape a joint attractor.

The successful future continuation already proves that at least one healthy
trajectory exists with

\[
\Gamma\to1,
\qquad
r_{\rm cone}\to0,
\]

while retaining positive \(F\), positive kinetic determinant, positive scalar
sound speed, No-Slip closure, and finite native perturbations.

What remains unknown is whether the two beta functions arise from one common
structural order parameter or from two coupled relaxation sectors.  This is
now the cleanest formulation of the remaining microscopic question.


## 70. The finite canonical handoff is genuinely two-dimensional

The two-beta formulation raises a natural question: can the successful
canonical trajectory actually be reduced to one global order parameter?

The completed closure-plane audit answers this in the negative for the full
future handoff.

Define

\[
u\equiv\ln\Gamma
\]

and

\[
v\equiv
\ln(1+2r_{\rm cone})
=
-\ln c_s^2.
\]

A particularly simple candidate common relaxation measure is

\[
\boxed{
{\cal L}_\lambda
=
\frac12
\left(
u^2+\lambda v^2
\right),
\qquad
\lambda>0.
}
\]

Its derivative is

\[
\frac{d{\cal L}_\lambda}{dN}
=
u\,u'
+
\lambda v\,v'.
\]

If one constant positive \(\lambda\) made this quantity non-positive over the
entire transition, the demonstrated future path would admit a simple global
quadratic Lyapunov distance to the canonical point.

The numerical audit finds that no such constant exists over the complete
interval beginning at the present boundary.

The pointwise inequalities require simultaneously

\[
\boxed{
\lambda
\gtrsim0.9452
}
\]

from one part of the trajectory, while another part requires

\[
\boxed{
\lambda
\lesssim3.07\times10^{-4}.
}
\]

These conditions are incompatible.

Therefore

\[
\boxed{
\text{the complete finite handoff is not described by a single
constant-metric quadratic relaxation distance in }
(\ln\Gamma,-\ln c_s^2).
}
\]

This is not a failure of the canonical endpoint.  It means that the rapid
matching layer contains genuinely multi-component dynamics.

### 70.1 Neither closure coordinate is globally monotone

The normalization coordinate itself is not monotonic through the complete
handoff.

Starting from

\[
\Gamma_0=0.91951178,
\]

the refined canonical path first rises to approximately

\[
\Gamma\simeq0.98455
\]

near

\[
\ln a\simeq0.07,
\]

then falls to approximately

\[
\Gamma\simeq0.92912
\]

near

\[
\ln a\simeq0.69,
\]

before turning again and finally approaching

\[
\Gamma_\infty=1.
\]

Thus the derivative of \(\Gamma\) changes sign twice.

The causal-shape coordinate is even more dramatic in the immediate matching
layer.  Although its accepted present value is

\[
r_{{\rm cone},0}\simeq3.11439,
\]

the rapid future coefficient release produces a short positive excursion to

\[
\boxed{
r_{\rm cone}\simeq9.97646
}
\]

near

\[
\ln a\simeq0.005,
\]

after which it rapidly collapses toward zero.

The scalar sound speed remains positive during this excursion.  The feature is
therefore a transient narrowing of the scalar cone, not a ghost or gradient
instability.

### 70.2 \(\Gamma\) cannot be the hidden common order parameter

The same value of \(\Gamma\) can occur at very different kinetic shapes.

For example, the audit finds two points with nearly equal normalization,

\[
\Gamma\simeq0.929,
\]

but with

\[
r_{\rm cone}\simeq9.98
\]

in the immediate handoff and

\[
r_{\rm cone}\simeq1.8\times10^{-2}
\]

later in the relaxation.

Hence

\[
\boxed{
r_{\rm cone}\neq r_{\rm cone}(\Gamma)
}
\]

as a single-valued global relation through the complete transition.

Therefore the normalization variable \(\Gamma\) by itself cannot serve as the
microscopic order parameter controlling both background and kinetic closure.

The reverse reduction also fails because \(r_{\rm cone}\) is not globally
monotonic.

### 70.3 Physical interpretation

The result divides the future evolution into two conceptually distinct stages.

First comes a **matching layer**, during which the accepted present Horndeski
jet is rearranged:

\[
F_{,\sigma},
\quad
G_3,
\quad
k_1,
\quad
k_2,
\quad
V
\]

all relax at different rates while continuity through the present boundary is
maintained.

During this stage, normalization and causal shape can move in different
directions.

After that layer, the higher Horndeski operators become small and the
trajectory enters the **mature relaxation regime**, where both closure
coordinates approach

\[
(\Gamma,r_{\rm cone})=(1,0).
\]

Thus a microscopic theory need not possess a one-dimensional order parameter
valid from the present epoch onward.

A more realistic possibility is

\[
\boxed{
\text{multi-field/multi-operator matching}
\quad\longrightarrow\quad
\text{lower-dimensional infrared relaxation}.
}
\]

This is compatible with the decoupling hierarchy of Section 66: several
finite-activation operators are important during the matching layer, but they
become irrelevant much faster than the residual matter eigenmode.

The next question is therefore more precise than asking for a global common
order parameter:

\[
\boxed{
\text{Does a genuine Lyapunov description emerge after the finite
matching layer has ended?}
}
\]

That question can be tested directly by moving the lower boundary of the
closure-plane audit forward in \(\ln a\) and searching for the first interval
on which a single positive quadratic metric becomes monotonic.


## 71. Collapse onto a one-dimensional canonical infrared manifold

The failure of a global one-parameter description in Section 70 applies to the
finite matching layer.  Once the higher Horndeski operators have decoupled,
the situation simplifies dramatically.

For the mature canonical action,

\[
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_\star
\right),
\]

with

\[
Z_\star=Z_{\rm BI},
\]

the positive structural source is exactly

\[
\rho_X
=
2ZG_{2,Z}-G_2
=
\frac{
2F_\infty
\left(
Z+2Z_\star
\right)
}{
\sigma^2
}.
\]

At the Balanced-Identity fixed point,

\[
Z=Z_\star,
\]

so

\[
\rho_{X,\star}
=
\frac{
6F_\infty Z_\star
}{
\sigma^2
}.
\]

The corresponding bare inverse-square structural source is therefore

\[
\boxed{
\rho_v
=
\frac{
6Z_\star
}{
\sigma^2
}.
}
\]

Consequently, throughout the canonical mature action,

\[
\Gamma
=
\frac{\rho_X/F_\infty}{\rho_v}
\]

obeys

\[
\boxed{
\Gamma_{\rm can}
=
\frac13
\left(
2+\frac{Z}{Z_\star}
\right).
}
\]

But the structural lapse satisfies

\[
{\cal N}
=
t_PS_0\sqrt{2Z},
\]

so

\[
\frac{Z}{Z_\star}
=
\left(
\frac{{\cal N}}{{\cal N}_\infty}
\right)^2.
\]

Hence the canonical infrared manifold is

\[
\boxed{
\Gamma_{\rm can}
=
\frac13
\left[
2+
\left(
\frac{{\cal N}}{{\cal N}_\infty}
\right)^2
\right].
}
\]

This is an exact action-level relation once the pure canonical mature normal
form has been reached.

At the same time,

\[
\boxed{
r_{\rm cone}=0
}
\]

identically.

Thus the two-dimensional closure plane collapses onto the one-dimensional
canonical manifold

\[
\boxed{
\left(
\Gamma,r_{\rm cone}
\right)
=
\left(
\frac{
2+({\cal N}/{\cal N}_\infty)^2
}{3},
\,0
\right).
}
\]

The finite handoff is therefore genuinely two-dimensional, but the infrared
relaxation is one-dimensional.

### 71.1 Matter controls the remaining coordinate

Using the resonant canonical expansion

\[
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-\frac32m_1a^{-1}
+
O(a^{-2}\ln a),
\]

the exact manifold relation gives

\[
\boxed{
\Gamma_{\rm can}
=
1
-
m_1a^{-1}
+
O(a^{-2}\ln a).
}
\]

At second order,

\[
\boxed{
\Gamma_{\rm can}
=
1
-m_1u
+
u^2
\left[
\frac94m_1^2\ln u
-
3m_1^2
-
\sqrt6\,s_2
\right]
+\cdots,
}
\]

where

\[
u=a^{-1}.
\]

Thus after the kinetic shape has relaxed,

\[
r_{\rm cone}=0,
\]

the remaining normalization displacement

\[
1-\Gamma
\]

is simply another representation of the residual matter eigenmode.

The late canonical theory therefore needs only one dynamical relaxation
coordinate.

### 71.2 The reconstructed trajectory reveals a constant calibration factor

The numerical structural source is normalized to the independently
reconstructed accepted present value

\[
\chi_0=0.941089569\ldots .
\]

When the refined canonical trajectory is compared with the ideal manifold, the
late samples satisfy

\[
\Gamma_{\rm action}
=
C_{\rm cal}
\,
\frac{
2+({\cal N}/{\cal N}_\infty)^2
}{3}
\]

with an essentially constant factor

\[
\boxed{
C_{\rm cal}
\simeq
0.9999704552.
}
\]

For example, the inferred factor is approximately

\[
0.9999704558
\]

at

\[
\ln a=5,
\]

\[
0.9999704552
\]

at

\[
\ln a=8,
\]

and

\[
0.9999704552
\]

at

\[
\ln a=10.
\]

Its stability over the mature interval shows that this is not another
late-time dynamical mode.

It is a fixed normalization offset between the reconstructed
present-calibrated \(\chi\) convention and the ideal canonical
Balanced-Identity normalization.

The size of the offset is

\[
\boxed{
C_{\rm cal}-1
\simeq
-2.95\times10^{-5}.
}
\]

The ideal continuum closure corresponds to

\[
C_{\rm cal}=1.
\]

Whether the remaining \(2.95\times10^{-5}\) should be absorbed into a refined
present normalization or retained as reconstruction uncertainty is a
calibration question, not a change of the mature equations of motion.

### 71.3 Three-stage structure of the future solution

The complete future development can therefore be organized into three stages.

First is the **multi-operator matching layer**,

\[
(\Gamma,r_{\rm cone})
\quad\text{genuinely two-dimensional},
\]

during which the present Horndeski jet is rearranged.

Second is the **canonicalization layer**,

\[
r_{\rm cone}\to0,
\qquad
F_{,\sigma}\to0,
\qquad
G_3\to0,
\qquad
k_2\to0.
\]

Third is the **matter-loaded canonical infrared flow**,

\[
\boxed{
r_{\rm cone}=0,
\qquad
\Gamma\to1,
}
\]

with

\[
1-\Gamma
\sim
m_1a^{-1}.
\]

This resolves the apparent tension between Sections 69 and 70.

A single global order parameter does not describe the complete handoff, but a
one-dimensional order parameter emerges naturally after canonicalization.

The remaining microscopic problem is therefore not required to make the
entire present-to-future trajectory one-dimensional.

It must instead explain why the multi-operator matching flow is attracted onto
the canonical infrared manifold.


## 72. Emergent Lyapunov flow and normal attraction of the canonical infrared manifold

The refined closure-plane audit can now be sharpened beyond the result of
Section 70.

The complete interval beginning exactly at the present boundary does not admit
one constant-metric quadratic Lyapunov function because the immediate matching
jet produces a very short causal-shape excursion.  However, this obstruction
is confined to the boundary layer itself.

Use the natural logarithmic closure coordinates

\[
u\equiv\ln\Gamma,
\qquad
v\equiv-\ln c_s^2
=
\ln(1+2r_{\rm cone}),
\]

and the equal-weight distance

\[
\boxed{
{\cal L}
=
\frac12
\left(
u^2+v^2
\right).
}
\]

The new moving-window audit finds that after

\[
\boxed{
\ln a\simeq0.005,
}
\]

the quantity \({\cal L}\) decreases monotonically throughout the tested
interval to

\[
\ln a=9.8.
\]

The corresponding scale factor is only

\[
\boxed{
a\simeq1.0050125,
}
\]

approximately a \(0.50\%\) increase beyond its present normalization.

For the present accepted \(H_0\), this corresponds roughly to

\[
\boxed{
\Delta t\sim7\times10^7\ {\rm yr},
}
\]

with the precise value depending on the short future variation of \(H\).

Thus the failure of a common Lyapunov distance is localized to an extremely
thin matching layer.  Once that layer is crossed, the demonstrated canonical
future path already behaves as a genuine relaxation flow even though
\(\Gamma\) itself still executes its broader non-monotonic excursion.

This distinction is important:

\[
\boxed{
\text{individual coordinates need not be monotonic for the total
distance to the attractor to decrease monotonically.}
}
\]

### 72.1 Tangent and transverse coordinates

After calibration of the tiny fixed normalization offset found in Section 71,
define the tangent displacement

\[
\boxed{
\epsilon_\parallel
\equiv
\left|
1-\frac{\Gamma}{C_{\rm cal}}
\right|,
}
\]

and the transverse kinetic-shape displacement

\[
\boxed{
\epsilon_\perp
\equiv
|r_{\rm cone}|.
}
\]

The first measures motion **along** the matter-loaded canonical infrared
manifold.

The second measures departure **away from** the canonical luminal surface.

The fitted tangent decay exponent approaches unity:

\[
\epsilon_\parallel\propto a^{-\nu_\parallel},
\]

with

\[
\nu_\parallel
=
0.9812
\quad
(\ln a=2\text{--}3),
\]

\[
\nu_\parallel
=
1.0097
\quad
(\ln a=3\text{--}4),
\]

and

\[
\nu_\parallel
=
1.0094
\quad
(\ln a=4\text{--}5).
\]

This reproduces the independently derived matter eigenvalue

\[
\boxed{
\nu_\parallel\to1.
}
\]

The transverse causal-shape mode decays much faster:

\[
\epsilon_\perp\propto a^{-\nu_\perp},
\]

with

\[
\nu_\perp
=
3.989
\quad
(\ln a=2\text{--}3),
\]

\[
\nu_\perp
=
4.342
\quad
(\ln a=3\text{--}4),
\]

and

\[
\nu_\perp
=
4.518
\quad
(\ln a=4\text{--}5).
\]

The local logarithmic rate reaches approximately

\[
\boxed{
\nu_\perp\simeq4.66
}
\]

by

\[
\ln a=6,
\]

while

\[
\nu_\parallel\simeq1.00.
\]

The numerical trend is therefore toward the expected transition-tail
hierarchy

\[
\boxed{
\nu_\perp\to5,
\qquad
\nu_\parallel\to1,
}
\]

consistent with the relative \(k_2\) suppression

\[
k_2Z^2/G_{2,\rm lead}\sim a^{-5}
\]

derived in Section 66.

### 72.2 Asymptotic normal attraction

The ratio of local contraction rates grows from approximately

\[
4.0
\]

near

\[
\ln a\simeq2
\]

to

\[
\boxed{
4.64
}
\]

by

\[
\ln a=6.
\]

Therefore departures transverse to the canonical manifold disappear several
times faster than the residual matter motion tangent to it.

This supplies numerical evidence that the canonical infrared manifold is
**asymptotically normally attractive**:

\[
\boxed{
|\nu_\perp|>|\nu_\parallel|.
}
\]

In dynamical-systems language, the fast variables

\[
F_{,\sigma},
\quad
G_3,
\quad
k_2,
\quad
r_{\rm cone}
\]

collapse rapidly onto a slow manifold, while the remaining matter-loaded
normalization coordinate evolves with the much smaller eigenvalue \(1\).

Thus the future hierarchy is not merely

\[
\text{two-dimensional}\to\text{one-dimensional}
\]

kinematically.  It is dynamically hierarchical:

\[
\boxed{
\text{fast canonicalization}
\quad\longrightarrow\quad
\text{slow matter relaxation}.
}
\]

This is precisely the structure expected if the mature canonical theory is an
infrared attractor rather than an externally imposed endpoint.

## 73. Exact canonical-manifold beta law and infrared Lyapunov exponent

The one-dimensional infrared manifold also permits an analytic flow equation.

Use the standard canonical exponential variables

\[
x
=
\frac{\dot{\widehat\phi}}{\sqrt6H},
\qquad
y
=
\frac{\sqrt V}{\sqrt3H},
\]

with

\[
\lambda_\ell=\sqrt2.
\]

For pressureless matter,

\[
x'
=
-3x+\sqrt3\,y^2
+\frac32x(1+x^2-y^2),
\]

\[
y'
=
-\sqrt3\,xy
+\frac32y(1+x^2-y^2).
\]

Define

\[
n
\equiv
\frac{{\cal N}}{{\cal N}_\infty}.
\]

The canonical structural identities give

\[
\boxed{
n=\frac{\sqrt2\,x}{y}
}
\]

and

\[
\boxed{
\Gamma
=
\frac{2+n^2}{3}.
}
\]

Taking a logarithmic derivative of \(n\),

\[
\frac{n'}{n}
=
\frac{x'}x-\frac{y'}y,
\]

and inserting the autonomous equations yields

\[
\frac{n'}n
=
-3
+
\sqrt3
\left(
x+\frac{y^2}{x}
\right).
\]

Since

\[
p=\sqrt3\,x
\]

and

\[
\frac{y^2}{x}
=
\frac{2x}{n^2},
\]

one obtains the exact dust-era canonical relation

\[
\boxed{
\frac{n'}n
=
-3
+
p
\left(
1+\frac2{n^2}
\right).
}
\]

Using

\[
n^2=3\Gamma-2,
\]

the normalization coordinate therefore obeys

\[
\boxed{
\Gamma'
=
2
\left[
2+(p-3)\Gamma
\right].
}
\]

This equation is exact on the pure canonical dust manifold.  Radiation adds
only a rapidly vanishing late-time correction.

### 73.1 Universal infrared beta function

The resonant matter expansion gives

\[
\Gamma
=
1-m_1a^{-1}
+
O(a^{-2}\ln a),
\]

and

\[
p
=
1-\frac32m_1a^{-1}
+
O(a^{-2}\ln a).
\]

Let

\[
\epsilon
\equiv
1-\Gamma.
\]

Then the exact manifold equation reduces to

\[
\boxed{
\epsilon'
=
-\epsilon
+
O(\epsilon^2\ln\epsilon).
}
\]

Hence the tangent infrared eigenvalue is exactly

\[
\boxed{
\lambda_{\rm IR}=1.
}
\]

The normalization beta function becomes

\[
\boxed{
\beta_\Gamma
=
\frac{d\ln\Gamma}{d\ln\sigma}
=
1-\Gamma
+
O\!\left(
(1-\Gamma)^2\ln(1-\Gamma)
\right).
}
\]

Thus the canonical infrared flow possesses a universal leading beta law
independent of the fitted interpolation rates.

Those rates control how rapidly the theory reaches the manifold; they do not
control the final matter relaxation once the manifold has been reached.

### 73.2 Analytic infrared Lyapunov function

A minimal one-dimensional Lyapunov function is

\[
\boxed{
{\cal V}_{\rm IR}
=
\frac12
(1-\Gamma)^2.
}
\]

Using

\[
\epsilon'=-\epsilon+\cdots,
\]

one obtains

\[
\frac{d{\cal V}_{\rm IR}}{dN}
=
\epsilon\epsilon'
=
-\epsilon^2
+
O(\epsilon^3\ln\epsilon).
\]

Therefore

\[
\boxed{
\frac{d{\cal V}_{\rm IR}}{dN}
=
-2{\cal V}_{\rm IR}
+
O({\cal V}_{\rm IR}^{3/2}\ln{\cal V}_{\rm IR}).
}
\]

To leading order,

\[
\boxed{
{\cal V}_{\rm IR}\propto a^{-2}.
}
\]

This \(a^{-2}\) decay of the Lyapunov distance is the square of the physical
matter eigenmode

\[
1-\Gamma\propto a^{-1}.
\]

The result supplies the analytic counterpart of the numerical closure-plane
audit:

\[
\boxed{
\text{the finite handoff is multi-operator, but the mature canonical
infrared flow has a simple one-dimensional Lyapunov law.}
}
\]

The remaining microscopic task is therefore reduced still further.  A deeper
SDMC completion must explain the fast transverse attraction onto the canonical
manifold.  Once that attraction has occurred, the subsequent normalization
flow is already fixed by the canonical action and ordinary matter
conservation.


## 74. Planck-braiding tail-root bifurcation and interpolation non-uniqueness

The C3 matching problem for the non-minimal gravitational sector has now been
solved analytically enough to expose a genuine branch structure.

Let

\[
x\equiv\ln\sigma,
\]

and denote by

\[
F_n
=
\left.
\frac{d^nF}{dx^n}
\right|_{x=0},
\qquad
g_n
=
\left.
\frac{d^ng}{dx^n}
\right|_{x=0}
\]

the accepted present C3 jets.

The future Planck-mass tail is written as

\[
F-F_\infty
=
e^{-\mu_Fx}P_F(x),
\]

while braiding decays as

\[
g
=
e^{-\mu_gx}P_g(x).
\]

Because

\[
F_{,\sigma}
=
\sigma^{-1}F_{,x},
\]

asymptotic No-Slip requires

\[
\boxed{
\mu_g=\mu_F+1.
}
\]

Matching the complete C3 present jet gives the cubic tail coefficients

\[
A_{F3}(\mu)
=
F_3
+
3\mu F_2
+
3\mu^2F_1
+
\mu^3(F_0-F_\infty),
\]

and

\[
A_{g3}(\mu+1)
=
g_3
+
3(\mu+1)g_2
+
3(\mu+1)^2g_1
+
(\mu+1)^3g_0.
\]

The leading asymptotic No-Slip relation

\[
g
=
-\frac{F_{,\sigma}}{2Z_\star}
\]

then reduces to the quartic condition

\[
\boxed{
\mu_F A_{F3}(\mu_F)
=
2Z_\star
A_{g3}(\mu_F+1).
}
\]

For the accepted present C3 jet and the canonical

\[
Z_\star=Z_{\rm BI},
\]

this quartic has two positive real roots,

\[
\boxed{
\mu_F^{(-)}
=
2.86106741391,
}
\]

and

\[
\boxed{
\mu_F^{(+)}
=
5.57321899721.
}
\]

The corresponding braiding exponents are

\[
\boxed{
\mu_g^{(-)}
=
3.86106741391,
}
\]

and

\[
\boxed{
\mu_g^{(+)}
=
6.57321899721.
}
\]

The remaining two quartic roots form a complex-conjugate pair and do not give
monotonic real exponential tails of the same form.

### 74.1 Both positive roots reach the same mature canonical action

The two positive roots were propagated independently through the exact
homogeneous action equations while keeping

\[
r=0,
\]

\[
\mu_{k_1}=50,
\qquad
\mu_{k_2}=5,
\qquad
\mu_V=120,
\]

and the same accepted present C3 jet.

At

\[
\ln a=10,
\]

the slow root gives approximately

\[
{\cal N}
=
2.8943616574,
\]

\[
p
=
0.9999850199,
\]

\[
q
=
-1.49816\times10^{-5},
\]

while the fast root gives

\[
{\cal N}
=
2.8943616597,
\]

\[
p
=
0.9999850207,
\]

\[
q
=
-1.49808\times10^{-5}.
\]

Both have

\[
D\simeq1.99994008,
\]

\[
c_s^2\simeq1,
\]

and

\[
\Gamma\simeq0.99996047.
\]

Thus the difference between the two roots is not a difference in mature
physics.

They are distinct finite interpolation branches connecting the same accepted
present jet to the same canonical infrared endpoint.

### 74.2 Native hi_class does not eliminate the slow branch

The slower root has now also been installed directly into the native hi_class
future table and tested independently.

Its native background reaches

\[
\ln a=5
\]

with

\[
\boxed{
{\cal N}=2.88788923481,
}
\]

\[
\boxed{
p=0.99775389658,
}
\]

\[
\boxed{
q=-0.00225372304,
}
\]

\[
\boxed{
D=1.9910269831,
}
\]

and

\[
\boxed{
c_s^2=0.999999521917.
}
\]

The corresponding fast branch gives

\[
{\cal N}=2.88789408597,
\]

\[
p=0.99775516727,
\]

\[
q=-0.00225182497,
\]

\[
D=1.99103075464,
\]

and

\[
c_s^2=0.999999997658.
\]

The differences are tiny.

The slow-root native future also satisfies the same background health gates,

\[
F>0,
\qquad
D>0,
\qquad
c_s^2>0,
\]

and the scalar perturbation propagation for

\[
k=
10^{-4},
\quad
10^{-3},
\quad
10^{-2}\ {\rm Mpc}^{-1}
\]

completes successfully with finite stored variables.

The native workflow returns

\[
\boxed{
\texttt{NATIVE\_FUTURE\_GATE\_PASS}
}
\]

with the slow branch included.

Therefore native perturbation propagation does **not** select the fast
\(\mu_F\simeq5.5732\) root over the slow
\(\mu_F\simeq2.8611\) root.

### 74.3 What is unique and what is not

The present investigation therefore separates two kinds of uniqueness.

The mature canonical endpoint is strongly fixed:

\[
\boxed{
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
}
\]

\[
\boxed{
G_{3,\infty}=0,
\qquad
G_{4,\infty}=\frac{F_\infty}{2}.
}
\]

Its structural lapse is

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}.
}
\]

Its kinetic shape is

\[
\boxed{
r_{\rm cone,\infty}=0.
}
\]

Its tangent infrared eigenvalue is

\[
\boxed{
\lambda_{\rm IR}=1.
}
\]

By contrast, the finite Planck/braiding handoff is not unique within the
current C3 exponential-polynomial ansatz.

At least two positive real branches satisfy

1. the same accepted present C3 jet;
2. the same asymptotic No-Slip relation;
3. the same canonical endpoint;
4. action-level background and health conditions;
5. native hi_class background evolution;
6. native scalar perturbation propagation.

Thus

\[
\boxed{
\text{unique endpoint}
\;\not\Rightarrow\;
\text{unique interpolation}.
}
\]

### 74.4 Consequence for the microscopic problem

The specific value

\[
\mu_F\simeq5.5732
\]

must therefore not be promoted to a fundamental SDMC constant on the basis of
the present evidence.

It is one healthy root of the C3 matching problem.

Likewise,

\[
\mu_F\simeq2.8611
\]

is another healthy root.

A deeper microscopic theory may ultimately select one branch, but if so it
must supply information not contained in

\[
\text{present C3 matching}
+
\text{asymptotic No-Slip}
+
\text{canonical endpoint}
+
\text{current background/perturbation health}.
\]

Possible additional selectors include a microscopic action principle, entropy
or dissipation law, analyticity requirement stronger than C3 matching,
UV completion, or a dynamical activation equation.

Until such a selector is derived, the scientifically correct statement is

\[
\boxed{
\text{the mature canonical action is much more constrained than the
finite path by which the accepted late-time Horndeski theory reaches it.}
}
\]

This is useful rather than problematic: it identifies exactly which parts of
the current future construction are physical predictions and which remain
interpolation freedom.


## 75. Canonical universality class versus finite transition memory

The two positive Planck/braiding roots found in Section 74 clarify which
quantities belong to the infrared theory itself and which quantities retain a
memory of the finite route into that theory.

The canonical endpoint is common to both branches:

\[
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}},
\]

\[
p_\infty=1,
\qquad
q_\infty=0,
\]

\[
D_\infty=2,
\qquad
c_{s,\infty}^2=1,
\]

\[
r_{\rm cone,\infty}=0.
\]

The tangent infrared eigenvalue is also common,

\[
\boxed{
\lambda_{\rm IR}=1.
}
\]

These are therefore **universal infrared data** of the mature canonical
action.

By contrast, the frozen mapper

\[
M_\infty
=
\lim_{a\to\infty}\frac{\sigma}{a}
\]

retains information about the finite transition.

Section 64 showed that the amplitude of the physical matter eigenmode is

\[
\boxed{
m_1
=
\frac{
\Omega_{m0}M_\infty^2
}{
F_\infty p_0^2
}
\left(
\frac{{\cal N}_0}{{\cal N}_\infty}
\right)^2.
}
\]

Therefore two healthy future histories that share the same accepted present
normalization and the same mature endpoint can still have slightly different
late matter-mode amplitudes if their frozen mapper memories differ.

For two branches \(A\) and \(B\),

\[
\boxed{
\frac{m_1^{(B)}}{m_1^{(A)}}
=
\left(
\frac{M_\infty^{(B)}}{M_\infty^{(A)}}
\right)^2.
}
\]

For a small mapper difference,

\[
\boxed{
\frac{\delta m_1}{m_1}
=
2\frac{\delta M_\infty}{M_\infty}
+
O\!\left[
\left(\frac{\delta M_\infty}{M_\infty}\right)^2
\right].
}
\]

The action-level branch comparison already shows this effect at
\(\ln a=10\).  The mapper values there are approximately

\[
M_{10}^{(-)}
=
0.7588342511
\]

for the slow Planck root and

\[
M_{10}^{(+)}
=
0.7588174881
\]

for the fast root.

Their fractional difference is only

\[
\boxed{
\frac{
M_{10}^{(+)}
}{
M_{10}^{(-)}
}
-1
\simeq
-2.21\times10^{-5}.
}
\]

The corresponding expected fractional shift in the asymptotic matter-mode
amplitude is therefore of order

\[
\boxed{
\frac{
m_1^{(+)}
}{
m_1^{(-)}
}
-1
\simeq
-4.42\times10^{-5},
}
\]

up to the tiny remaining difference between \(M_{10}\) and the exact frozen
\(M_\infty\).

This is small enough that the two branches are essentially
indistinguishable at the level of the mature background endpoint, while still
being mathematically distinct histories.

### 75.1 Universal exponents, non-universal amplitudes

The canonical matter-loaded infrared flow has the generic form

\[
1-\Gamma
=
m_1a^{-1}
+\cdots,
\]

\[
1-\frac{{\cal N}}{{\cal N}_\infty}
=
\frac32m_1a^{-1}
+\cdots.
\]

The exponent

\[
a^{-1}
\]

is universal because it is fixed by the canonical dust eigenvalue.

The coefficient

\[
m_1
\]

is not fully universal because it contains \(M_\infty\), which remembers the
finite transition.

Thus the clean dynamical distinction is

\[
\boxed{
\text{infrared exponent}
=
\text{universal},
}
\]

while

\[
\boxed{
\text{infrared amplitude}
=
\text{transition-memory dependent}.
}
\]

This is the same pattern familiar in many attractor systems: different
trajectories can fall into the same fixed point with the same critical
exponents while carrying different amplitudes along the allowed decaying
modes.

### 75.2 The intrinsic scalar amplitude can also retain path memory

The resonant second-order expansion contains an independent scalar-mode
coefficient \(s_2\),

\[
a^{-2}
\]

together with the matter-generated resonant contribution

\[
a^{-2}\ln a.
\]

The eigenvalue \(-2\) is fixed by the mature canonical action, but its free
homogeneous amplitude is set by how the finite transition enters the canonical
basin.

Therefore the expected hierarchy is

\[
\boxed{
\lambda_1=1,
\qquad
\lambda_2=2
}
\]

as universal infrared eigenvalues, while

\[
\boxed{
m_1,
\qquad
s_2
}
\]

can retain small branch-dependent memories.

The dedicated tail-root memory audit now measures these quantities
independently for both positive roots.

### 75.3 What a microscopic selector would have to determine

A microscopic theory that selects one Planck/braiding root over the other
would therefore not be changing the mature canonical Lagrangian.

It would be selecting the **basin-entry history**, and hence the amplitudes
with which the universal infrared modes are excited.

That is a much narrower task than deriving the endpoint itself.

The current picture is therefore

\[
\boxed{
\text{many admissible finite handoffs}
\longrightarrow
\text{one canonical infrared universality class}.
}
\]

A future microscopic principle may choose a preferred handoff, but the
canonical endpoint and its physical infrared exponents do not depend on that
choice.


## 76. Long-horizon native validation selects the fast root operationally

The two positive C3 Planck/braiding roots remain almost indistinguishable in
the exact homogeneous action equations, and both survive the native
\(\ln a=5\) background and perturbation tests.

The stricter native extension to

\[
\ln a=10
\]

now separates them operationally.

### 76.1 Fast-root native closure remains successful

For the fast branch,

\[
\mu_F^{(+)}
=
5.57321899721,
\qquad
\mu_g^{(+)}
=
6.57321899721,
\]

the native background reaches

\[
\ln a=10
\]

with

\[
\boxed{
{\cal N}
=
2.89436190703,
}
\]

\[
\boxed{
p
=
0.99998501061,
}
\]

\[
\boxed{
q
=
-1.46175\times10^{-5},
}
\]

\[
\boxed{
D
=
1.99994003958,
}
\]

and

\[
\boxed{
c_s^2
=
0.9999999999994.
}
\]

The three native scalar perturbation modes

\[
k=
10^{-4},
\quad
10^{-3},
\quad
10^{-2}\ {\rm Mpc}^{-1}
\]

also remain finite through the same horizon and satisfy the long-horizon metric
potential gates.

Thus the fast branch again returns

\[
\boxed{
\texttt{NATIVE\_CANONICAL\_LONG\_HORIZON\_GATE\_PASS}.
}
\]

### 76.2 Slow-root native background loses viability before \(\ln a=10\)

For the slower branch,

\[
\mu_F^{(-)}
=
2.86106741391,
\qquad
\mu_g^{(-)}
=
3.86106741391,
\]

the same native installation is healthy through

\[
\ln a=5,
\]

but the attempted extension to

\[
\ln a=10
\]

terminates inside the native background solver with

\[
\boxed{
\rho_{\rm crit}
=
-1.466883\times10^{-14},
}
\]

where hi_class requires

\[
\rho_{\rm crit}>0.
\]

The corresponding native long-horizon background file is therefore not
produced, and the slow-root long-horizon gate returns failure.

This is the first test in the investigation that distinguishes the two
positive Planck/braiding roots.

### 76.3 This is an operational, not yet fundamental, branch selection

The result must be interpreted carefully.

The exact action-level slow-root trajectory remains healthy through

\[
\ln a=10,
\]

with

\[
{\cal N}
=
2.89436165740,
\]

\[
p
=
0.99998501991,
\]

\[
q
=
-1.49816\times10^{-5},
\]

\[
D
\simeq
1.99994008,
\]

and

\[
c_s^2
\simeq
1.
\]

It also passes native background and perturbation propagation through

\[
\ln a=5.
\]

Therefore the present evidence does **not** prove that the slower root is
mathematically inconsistent as a covariant action.

Instead it proves the narrower statement

\[
\boxed{
\text{the current native tabulated hi\_class realization of the slow root
does not remain background-viable to }\ln a=10.
}
\]

The failure could represent

1. a genuine incompatibility that appears only in the complete native
background system;
2. accumulated mismatch between the exact structural action and its tabulated
logstruct realization;
3. a conditioning/resolution issue associated with the much longer
Planck/braiding transient.

These possibilities require a dedicated convergence and failure-onset audit.

### 76.4 Production-branch status

Until that audit demonstrates otherwise, the branch hierarchy is now

\[
\boxed{
\mu_F^{(+)}
=
5.57321899721
}
\]

as the **fully validated production branch**, because it passes

- exact homogeneous action evolution;
- background health;
- No-Slip closure;
- action-level perturbation-envelope checks;
- native background propagation;
- native scalar perturbations;
- native \(\ln a=10\) long-horizon closure.

The slower

\[
\boxed{
\mu_F^{(-)}
=
2.86106741391
}
\]

remains a mathematically healthy action-level alternative and a successful
short-horizon native branch, but it is not yet long-horizon-native validated.

Thus the mature canonical endpoint remains branch independent, while the
current implementation evidence now gives a practical preference for the
faster Planck/braiding decay.

### 76.5 Refined universality statement

The branch-memory audit simultaneously confirms that the two roots excite
nearly the same universal infrared matter mode.

For the slow branch,

\[
M_\infty
=
0.75882288887,
\]

\[
m_1
=
0.21994106005,
\]

while for the fast branch,

\[
M_\infty
=
0.75880612254,
\]

\[
m_1
=
0.21993129363.
\]

Hence

\[
\boxed{
\frac{\Delta M_\infty}{M_\infty}
=
-2.21\times10^{-5},
}
\]

and

\[
\boxed{
\frac{\Delta m_1}{m_1}
=
-4.44\times10^{-5}.
}
\]

Both branches nevertheless recover the same tangent exponent,

\[
\boxed{
\nu_\parallel\to1.
}
\]

Their transverse causal-shape relaxation is different.

The fast branch approaches approximately

\[
\nu_\perp
\simeq
4.62
\]

over

\[
5\le\ln a\le6,
\]

continuing toward the \(k_2\)-controlled value near \(5\).

The slow branch instead gives approximately

\[
\boxed{
\nu_\perp
\simeq
2.23
}
\]

over the same interval.

Thus the physical matter eigenvalue is genuinely universal, whereas the rate
at which the finite noncanonical sector collapses onto the canonical manifold
retains substantial interpolation-branch dependence.

This sharpens the distinction:

\[
\boxed{
\text{canonical IR dynamics}
=
\text{universal},
}
\]

but

\[
\boxed{
\text{approach to the canonical IR manifold}
=
\text{branch dependent}.
}
\]

The next audit should locate the slow-root native failure boundary and test its
convergence with future-table resolution before deciding whether the
long-horizon failure is physical or numerical.


## 77. Consolidated closed equation set for the final discussion

This section collects the equations that now define the mature SDMC
structural-clock Lagrangian, its canonical fixed point, its matter-loaded
approach, and the finite branch structure.  The purpose is to separate
quantities that are fully derived from those that remain implementation or
microscopic-selection questions.

### 77.1 Structural variables and exact kinematics

Normalize the structural scale by

\[
\boxed{
\sigma
\equiv
\frac{S}{S_0}
=
e^\psi .
}
\]

The structural radius is

\[
\boxed{
R
=
\ell_P S
=
\ell_P S_0\sigma .
}
\]

For a homogeneous structural field define

\[
\boxed{
Z
\equiv
-\frac12
g^{\mu\nu}
\partial_\mu\sigma
\partial_\nu\sigma
=
\frac12\dot\sigma^2 .
}
\]

The structural lapse is

\[
\boxed{
{\cal N}
=
t_PS_0\dot\sigma
=
\frac{\dot R}{c}.
}
\]

The bridge exponent is

\[
\boxed{
p
\equiv
\frac{d\ln\sigma}{d\ln a}
=
\frac{\dot\sigma}{H\sigma}.
}
\]

Therefore

\[
\boxed{
\frac{HR}{c}
=
\frac{{\cal N}}{p}.
}
\]

Define the normalized matter/baryon mapper

\[
\boxed{
M(a)
\equiv
\frac{\sigma}{a}
=
\frac{K_m(a)}{K_{m0}}
=
\frac{K_p(a)}{K_{p0}}.
}
\]

It obeys the exact differential identity

\[
\boxed{
\frac{d\ln M}{d\ln a}
=
p-1.
}
\]

Hence

\[
\boxed{
M_\infty
=
\exp\left[
\int_0^\infty
(p-1)\,dN
\right].
}
\]

### 77.2 Horndeski structural action used in the reconstruction

The accepted finite late-time action is represented in the linear-\(G_3\)
Horndeski sector as

\[
\boxed{
S
=
\int d^4x\sqrt{-g}
\left[
G_2(\sigma,Z)
-
G_3(\sigma,Z)\Box\sigma
+
G_4(\sigma)R
\right]
+
S_m .
}
\]

The reconstructed coefficient form is

\[
\boxed{
G_2
=
k_1(\sigma)Z
+
k_2(\sigma)Z^2
-
V(\sigma),
}
\]

\[
\boxed{
G_3
=
g(\sigma)Z,
}
\]

\[
\boxed{
G_4
=
\frac{F(\sigma)}{2}.
}
\]

The No-Slip trajectory condition is

\[
\boxed{
\alpha_B
=
-2\alpha_M.
}
\]

The scalar health conditions are

\[
\boxed{
F>0,
\qquad
D>0,
\qquad
c_s^2>0.
}
\]

### 77.3 Covariant Balanced Identity

Define the mature normalization coordinate

\[
\boxed{
\Gamma
\equiv
\frac{\chi}{F}.
}
\]

The covariant Unified Balanced Identity is

\[
\boxed{
{\cal I}_{\rm BI}
\equiv
\frac{
G_{\rm eff}\rho_XR^2
}{
c^2
}
=
1.
}
\]

For the mature structural solution this gives

\[
\boxed{
\Gamma_\infty=1.
}
\]

The corresponding structural lapse is

\[
\boxed{
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}
=
2.89440501823\ldots .
}
\]

The fixed kinetic scale is

\[
\boxed{
Z_{\rm BI}
=
\frac{
4\pi
}{
3t_P^2S_0^2
}.
}
\]

Indeed,

\[
{\cal N}_\infty^2
=
2t_P^2S_0^2Z_{\rm BI}
=
\frac{8\pi}{3}.
\]

### 77.4 General scale-covariant mature family

The mature constant-\(F\), unbraided action has the scale-covariant form

\[
\boxed{
G_{2,\infty}
=
\frac{f(Z)}{\sigma^2},
}
\]

\[
\boxed{
G_{3,\infty}=0,
\qquad
G_{4,\infty}
=
\frac{F_\infty}{2}.
}
\]

The coasting fixed-point equations are

\[
\boxed{
f(Z_\star)
+
Z_\star f_{,Z}(Z_\star)
=
0,
}
\]

and

\[
\boxed{
f_{,Z}(Z_\star)
=
2F_\infty.
}
\]

For the quadratic family,

\[
f(Z)
=
\kappa_1 Z
+
\kappa_2 Z^2
-
U_0,
\]

write

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

The exact mature scalar sound speed is

\[
\boxed{
c_s^2(Z)
=
\frac{
f_{,Z}
}{
f_{,Z}
+
2Zf_{,ZZ}
}.
}
\]

At the fixed point,

\[
\boxed{
c_{s,\star}^2
=
\frac{1}{1+2r}.
}
\]

### 77.5 Causal-shape modulus and canonical selection

Define

\[
\boxed{
r_{\rm cone}
\equiv
\frac12
\left(
\frac1{c_s^2}-1
\right).
}
\]

In the mature constant-\(F\), unbraided theory,

\[
\boxed{
r_{\rm cone}
=
\frac{
Zf_{,ZZ}
}{
f_{,Z}
}.
}
\]

For the quadratic family at the fixed point,

\[
\boxed{
r_{\rm cone}=r.
}
\]

Robust luminality over a neighborhood of the matter-loaded fixed point,

\[
\boxed{
c_s^2(Z)=1
}
\]

for an open interval of \(Z\), implies

\[
\boxed{
f_{,ZZ}=0,
}
\]

and therefore

\[
\boxed{
f(Z)=AZ-B.
}
\]

The coasting equations then give

\[
\boxed{
A=2F_\infty,
}
\]

\[
\boxed{
B=4F_\infty Z_{\rm BI}.
}
\]

Thus the mature canonical action is

\[
\boxed{
G_{2,\infty}
=
\frac1{\sigma^2}
\left(
2F_\infty Z
-
4F_\infty Z_{\rm BI}
\right),
}
\]

\[
\boxed{
G_{3,\infty}=0,
}
\]

\[
\boxed{
G_{4,\infty}
=
\frac{F_\infty}{2}.
}
\]

Equivalently,

\[
\boxed{
r=0,
\qquad
r_{\rm cone}=0,
\qquad
c_s^2=1.
}
\]

### 77.6 Canonical-field representation

Define the canonically shaped logarithmic structural field

\[
\boxed{
\widehat\phi
=
\sqrt2\ln\sigma.
}
\]

Then

\[
\sigma^{-2}
=
e^{-\sqrt2\widehat\phi},
\]

and the mature potential is

\[
\boxed{
V_\infty(\widehat\phi)
=
4F_\infty Z_{\rm BI}
e^{-\sqrt2\widehat\phi}.
}
\]

Hence the exponential slope is

\[
\boxed{
\lambda_\ell
=
\sqrt2.
}
\]

### 77.7 Exact mature fixed point

At the canonical fixed point,

\[
\boxed{
Z=Z_{\rm BI}.
}
\]

Therefore

\[
\boxed{
\dot\sigma
=
\sqrt{2Z_{\rm BI}}
=
\frac1{t_PS_0}
\sqrt{\frac{8\pi}{3}}.
}
\]

The structural lapse becomes

\[
\boxed{
{\cal N}
=
{\cal N}_\infty
=
\sqrt{\frac{8\pi}{3}}.
}
\]

The structural field evolves linearly,

\[
\boxed{
\sigma(t)
=
\sigma_c
+
\frac{
\sqrt{8\pi/3}
}{
t_PS_0
}
(t-t_c).
}
\]

Since

\[
M=\sigma/a\to M_\infty,
\]

one has

\[
\boxed{
a\propto\sigma\propto t,
}
\]

\[
\boxed{
H
=
\frac{\dot\sigma}{\sigma},
}
\]

and

\[
\boxed{
q=0.
}
\]

The mature scalar density and pressure are

\[
\boxed{
\rho_{\sigma,\infty}
=
\frac{
6F_\infty Z_{\rm BI}
}{
\sigma^2
},
}
\]

\[
\boxed{
p_{\sigma,\infty}
=
-\frac{
2F_\infty Z_{\rm BI}
}{
\sigma^2
}.
}
\]

Hence

\[
\boxed{
w_{\sigma,\infty}
=
-\frac13.
}
\]

The canonical health limits are

\[
\boxed{
D_\infty=2,
\qquad
c_{s,\infty}^2=1.
}
\]

### 77.8 Autonomous fixed-point system

Using

\[
\boxed{
x
=
\frac{
\dot{\widehat\phi}
}{
\sqrt6H
},
}
\]

and

\[
\boxed{
y
=
\frac{
\sqrt V
}{
\sqrt3H
},
}
\]

the dust-era system for

\[
\lambda_\ell=\sqrt2
\]

is

\[
\boxed{
x'
=
-3x+\sqrt3\,y^2
+
\frac32x
\left(
1+x^2-y^2
\right),
}
\]

\[
\boxed{
y'
=
-\sqrt3\,xy
+
\frac32y
\left(
1+x^2-y^2
\right).
}
\]

The scalar-dominated fixed point is

\[
\boxed{
x_\star
=
\frac1{\sqrt3},
}
\]

\[
\boxed{
y_\star
=
\sqrt{\frac23}.
}
\]

Its eigenvalues are

\[
\boxed{
\mu_1=-1,
\qquad
\mu_2=-2.
}
\]

The \(-1\) mode is the physical residual-matter mode and the \(-2\) mode is the
intrinsic scalar mode.

### 77.9 Resonant matter-loaded asymptotics

Let

\[
\boxed{
u=a^{-1}.
}
\]

Because

\[
2(-1)=-2,
\]

the matter mode is resonant with the intrinsic scalar mode and second-order
terms contain

\[
u^2\ln u.
\]

Define the leading matter amplitude \(m_1\) and the independent intrinsic
scalar amplitude \(s_2\).

Then

\[
\boxed{
\Omega_m
=
m_1u
+
3m_1^2u^2
+
O(u^3\ln u).
}
\]

The bridge exponent is

\[
\boxed{
p
=
1
-
\frac32m_1u
+
u^2
\left[
\frac94m_1^2\ln u
-
\frac{45}{8}m_1^2
-
\sqrt6\,s_2
\right]
+\cdots .
}
\]

The deceleration parameter is

\[
\boxed{
q
=
-\frac32m_1u
+
u^2
\left[
\frac92m_1^2\ln u
-
\frac92m_1^2
-
2\sqrt6\,s_2
\right]
+\cdots .
}
\]

The structural lapse is

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-
\frac32m_1u
+
u^2
\left[
\frac{27}{8}m_1^2\ln u
-
\frac{45}{8}m_1^2
-
\frac{3\sqrt6}{2}s_2
\right]
+\cdots .
}
\]

The mapper approaches

\[
\boxed{
\frac{M}{M_\infty}
=
1
+
\frac32m_1u
+
u^2
\left[
-\frac98m_1^2\ln u
+
\frac92m_1^2
+
\frac{\sqrt6}{2}s_2
\right]
+\cdots .
}
\]

Therefore the complete \(O(a^{-1})\) correction cancels in

\[
{\cal N}M,
\]

leaving

\[
\boxed{
\frac{
{\cal N}M
}{
{\cal N}_\infty M_\infty
}
=
1
+
u^2
\left[
\frac94m_1^2\ln u
-
\frac{27}{8}m_1^2
-
\sqrt6\,s_2
\right]
+\cdots .
}
\]

### 77.10 Mapper-memory closure of the matter eigenmode

The asymptotic matter amplitude is not free.

It is fixed by

\[
\boxed{
m_1
=
\frac{
\Omega_{m0}M_\infty^2
}{
F_\infty p_0^2
}
\left(
\frac{{\cal N}_0}{{\cal N}_\infty}
\right)^2 .
}
\]

For the current fast production branch,

\[
\boxed{
M_\infty
=
0.75880612254,
}
\]

and

\[
\boxed{
m_1
=
0.21993129363.
}
\]

Hence the leading canonical approach is

\[
\boxed{
\frac{{\cal N}}{{\cal N}_\infty}
=
1
-
0.32989694044\,a^{-1}
+
O(a^{-2}\ln a),
}
\]

with the exact coefficient

\[
\frac32m_1
=
0.32989694044\ldots .
\]

### 77.11 Exact one-dimensional canonical infrared manifold

For the pure canonical mature action,

\[
\rho_X
=
\frac{
2F_\infty
\left(
Z+2Z_{\rm BI}
\right)
}{
\sigma^2
}.
\]

The canonical closure coordinate therefore obeys

\[
\boxed{
\Gamma_{\rm can}
=
\frac13
\left(
2+\frac{Z}{Z_{\rm BI}}
\right).
}
\]

Since

\[
\frac{Z}{Z_{\rm BI}}
=
\left(
\frac{{\cal N}}{{\cal N}_\infty}
\right)^2,
\]

the exact manifold is

\[
\boxed{
\Gamma_{\rm can}
=
\frac13
\left[
2+
\left(
\frac{{\cal N}}{{\cal N}_\infty}
\right)^2
\right],
}
\]

with

\[
\boxed{
r_{\rm cone}=0.
}
\]

Thus

\[
\boxed{
(\Gamma,r_{\rm cone})
=
\left(
\frac{
2+({\cal N}/{\cal N}_\infty)^2
}{3},
0
\right)
}
\]

after canonicalization.

The matter-loaded expansion is

\[
\boxed{
\Gamma
=
1
-
m_1a^{-1}
+
O(a^{-2}\ln a).
}
\]

### 77.12 Exact canonical beta law

Define

\[
\boxed{
n
=
\frac{{\cal N}}{{\cal N}_\infty}.
}
\]

The canonical variables obey

\[
\boxed{
n
=
\frac{\sqrt2\,x}{y}.
}
\]

The exact dust-era relation is

\[
\boxed{
\frac{n'}{n}
=
-3
+
p
\left(
1+\frac2{n^2}
\right).
}
\]

Using

\[
n^2=3\Gamma-2,
\]

one obtains

\[
\boxed{
\Gamma'
=
2
\left[
2+(p-3)\Gamma
\right].
}
\]

For

\[
\epsilon
=
1-\Gamma,
\]

the infrared limit becomes

\[
\boxed{
\epsilon'
=
-\epsilon
+
O(\epsilon^2\ln\epsilon).
}
\]

Therefore the universal tangent eigenvalue is

\[
\boxed{
\lambda_{\rm IR}=1.
}
\]

The corresponding one-dimensional Lyapunov function is

\[
\boxed{
{\cal V}_{\rm IR}
=
\frac12(1-\Gamma)^2,
}
\]

with

\[
\boxed{
\frac{d{\cal V}_{\rm IR}}{dN}
=
-2{\cal V}_{\rm IR}
+
O(
{\cal V}_{\rm IR}^{3/2}
\ln{\cal V}_{\rm IR}
).
}
\]

### 77.13 Full two-coordinate closure before canonicalization

Before the higher Horndeski operators have decoupled, use

\[
\boxed{
(\Gamma,r_{\rm cone}).
}
\]

The normalization beta function is

\[
\boxed{
\beta_\Gamma
\equiv
\frac{d\ln\Gamma}{d\ln\sigma}
=
2-m_X-\frac{\alpha_M}{p}.
}
\]

Hence

\[
\boxed{
\frac{d\Gamma}{dN}
=
p\Gamma\beta_\Gamma.
}
\]

The causal-shape flow is

\[
\boxed{
\frac{dr_{\rm cone}}{dN}
=
-\frac{
1+2r_{\rm cone}
}{2}
\frac{d\ln c_s^2}{dN}.
}
\]

Equivalently,

\[
\boxed{
\beta_r
=
-\frac{
1+2r_{\rm cone}
}{
2p
}
\frac{d\ln c_s^2}{dN}.
}
\]

The canonical endpoint is

\[
\boxed{
(\Gamma,r_{\rm cone})
\to
(1,0).
}
\]

The complete finite handoff is genuinely two-dimensional, but after the
short matching/canonicalization layer the flow collapses onto the
one-dimensional manifold of Section 77.11.

### 77.14 Planck-braiding branch equation

The C3 future Planck tail is

\[
F-F_\infty
=
e^{-\mu_Fx}P_F(x),
\]

with

\[
x=\ln\sigma,
\]

while

\[
g
=
e^{-\mu_gx}P_g(x).
\]

Asymptotic No-Slip requires

\[
\boxed{
\mu_g=\mu_F+1.
}
\]

The C3 matching equation is

\[
\boxed{
\mu_F A_{F3}(\mu_F)
=
2Z_{\rm BI}
A_{g3}(\mu_F+1),
}
\]

where

\[
A_{F3}(\mu)
=
F_3
+
3\mu F_2
+
3\mu^2F_1
+
\mu^3(F_0-F_\infty),
\]

and

\[
A_{g3}(\mu+1)
=
g_3
+
3(\mu+1)g_2
+
3(\mu+1)^2g_1
+
(\mu+1)^3g_0.
\]

The two positive real roots are

\[
\boxed{
\mu_F^{(-)}
=
2.86106741391,
}
\]

and

\[
\boxed{
\mu_F^{(+)}
=
5.57321899721.
}
\]

Their braiding exponents are

\[
\boxed{
\mu_g^{(-)}
=
3.86106741391,
}
\]

and

\[
\boxed{
\mu_g^{(+)}
=
6.57321899721.
}
\]

Both reach the same action-level canonical endpoint.

The current fully native-validated production branch is

\[
\boxed{
\mu_F^{(+)}
=
5.57321899721.
}
\]

The slow branch remains action-level healthy but is subject to the dedicated
native convergence audit described after Section 76.

### 77.15 Production numerical closure

For the fast production branch, native hi_class reaches

\[
\ln a=10
\]

with approximately

\[
\boxed{
{\cal N}
=
2.89436191,
}
\]

\[
\boxed{
p
=
0.99998501,
}
\]

\[
\boxed{
q
=
-1.46\times10^{-5},
}
\]

\[
\boxed{
D
=
1.99994004,
}
\]

and

\[
\boxed{
c_s^2
=
0.999999999999.
}
\]

Thus

\[
\boxed{
{\cal N}
\to
\sqrt{\frac{8\pi}{3}},
\qquad
p\to1,
\qquad
q\to0,
\qquad
D\to2,
\qquad
c_s^2\to1.
}
\]

These equations constitute the current closed mathematical core of the
Lagrangian development.

What remains outside this closed core is not the mature action itself, but

1. the microscopic origin of the transverse attraction onto the canonical
   manifold;
2. the possible microscopic selection of one finite C3 handoff branch;
3. the numerical convergence classification of the slow-root native
   long-horizon failure.

Those points should be treated separately in the final discussion rather than
mixed into the derived mature equations.
