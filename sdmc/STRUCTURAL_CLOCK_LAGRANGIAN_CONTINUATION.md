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
