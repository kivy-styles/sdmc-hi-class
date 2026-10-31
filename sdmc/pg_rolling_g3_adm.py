"""Exact unitary-gauge ADM map for the accepted rolling G2/G3/F scalar action.

Covariant sector:
 S_J = int sqrt(-g)[G2(sigma,Z)-G3(sigma,Z) Box(sigma)+F(sigma)R/2],
 G2=k1 Z+k2 Z²-V, G3=g Z, Z=-g^ab sig_a sig_b/2.

On monotonic sigma(t_c), with ADM lapse N and Z=v²/(2N²),
up to covariant boundary terms, the *exact unitary-gauge ADM*
density is N sqrt(h)[A2(N)+A3(N)K+F/2(R3+KijKij-K²)]:
 A2 = k1 Z + (k2-2g_sigma/3)Z²-V
 A3 = -F_sigma v/N + g v³/(3N³).

This implements (i) the lapse/shift QUADRATIC kinetic Hessian in flat
FLRW unitary gauge, (ii) the kinetic coefficient Qs analytically and
numerically, (iii) the F/G3 Horndeski alpha EFT map and the exact
identity Qs=2 F D/(2-alphaB)²,
(iv) an isolated scalar-only gradient checkpoint and a conditional
No-Slip background gradient gate.

It is not the full coupled supporting-space ADM action with dynamical
L_ij, curvature W, and geometric averaged stress. No CMB fit or CDM
replacement follows from the kinetic identities alone.
"""
from dataclasses import dataclass
import json
import numpy as np

@dataclass(frozen=True)
class RollingAction:
    F:float
    F_sigma:float
    H:float
    sigma_dot:float
    k1:float
    k2:float
    g:float
    g_sigma:float
    V:float=0.
    a:float=1.

    def __post_init__(self):
        vals=np.asarray(list(vars(self).values()),dtype=float)
        if not np.isfinite(vals).all():raise ValueError("finite action coefficients")
        if self.F<=0 or self.a<=0 or self.H==0:
            raise ValueError("F,a>0 and H nonzero required")
        if self.sigma_dot<=0:
            raise ValueError("unitary ADM orientation assumes sigma_dot>0")

    @property
    def Z0(self):return .5*self.sigma_dot**2

def coefficients_at_lapse(action,N):
    """The exact A2 and A3 for spatially varying positive ADM lapse."""
    N=np.asarray(N,dtype=float)
    if np.any(N<=0):raise ValueError("ADM lapse N must remain >0")
    X=action.Z0/N**2
    A2=action.k1*X+(action.k2-2*action.g_sigma/3)*X**2-action.V
    A3=(-action.F_sigma*action.sigma_dot/N+
        action.g*action.sigma_dot**3/(3*N**3))
    return A2,A3

def lapse_derivatives(action):
    """A2, A3 and first two lapse derivatives evaluated at N=1."""
    X=action.Z0
    p=action.k2-2*action.g_sigma/3
    v=action.sigma_dot
    A20=action.k1*X+p*X**2-action.V
    A21=-2*action.k1*X-4*p*X**2
    A22=6*action.k1*X+20*p*X**2
    A30=-action.F_sigma*v+action.g*v**3/3
    A31=action.F_sigma*v-action.g*v**3
    A32=-2*action.F_sigma*v+4*action.g*v**3
    return dict(A20=A20,A21=A21,A22=A22,A30=A30,
                A31=A31,A32=A32)

def alphas(action):
    """Bellini-Sawicki EFT parameters for the accepted G2/G3/F class."""
    F=action.F;H=action.H;v=action.sigma_dot;X=action.Z0
    am=v*action.F_sigma/(H*F)
    ab=v*(2*X*action.g-action.F_sigma)/(H*F)
    ak=(2*X/(H**2*F))*(action.k1+6*action.k2*X-
                           4*action.g_sigma*X)+12*v*X*action.g/(H*F)
    D=ak+1.5*ab**2
    return dict(alpha_M=am,alpha_B=ab,alpha_K=ak,D=D,
                no_slip_residual=ab+2*am)

def adm_kinetic_hessian(action):
    """Exact quadratic ADM Hessian on (alpha, B, zeta_dot) in a cos(kx) mode.

    N=1+alpha cos(kx), shift potential beta with
    B=(k²/a²)beta_amplitude, spatial metric
    h_ij=a² exp[2 zeta cos(kx)] delta_ij.

    For zeta amplitude set to zero, the spatially averaged quadratic
    ADM Lagrangian density divided by a³ is

      L2 = U alpha² + 3 Theta alpha zeta_dot
            + Theta alpha B - F B zeta_dot -3F zeta_dot²/2.

    Theta=F H + A3_N/2,
    U = A2_N/2 + A2_NN/4 +3H A3_NN/4 -3F H²/2.

    Hessian is with respect to AMPLITUDES for the above cosine mode.
    This includes the coefficient 1/2 from averaging cos² over space.
    """
    d=lapse_derivatives(action)
    Theta=action.F*action.H+.5*d["A31"]
    U=.5*d["A21"]+.25*d["A22"]+.75*action.H*d["A32"]-1.5*action.F*action.H**2
    H=np.array([[2*U,Theta,3*Theta],
                [Theta,0.,-action.F],
                [3*Theta,-action.F,-3*action.F]],dtype=float)
    return dict(Hessian=H,Theta=Theta,U=U)

def kinetic_from_constraints(action):
    """Eliminate both lapse and shift in the scalar-only unitary ADM kinetic block.

    Shift gives alpha=F zeta_dot/Theta; lapse then sets B.
    Scalar Q_s follows from averaged L2=Q_s zeta_dot²/2.
    """
    b=adm_kinetic_hessian(action)
    Theta=b["Theta"];U=b["U"]
    if abs(Theta)<1e-12:
        raise ValueError("Theta=0: degenerate lapse/shift constraint")
    alpha_per_zdot=action.F/Theta
    B_per_zdot=-(2*U*alpha_per_zdot+3*Theta)/Theta
    qs_adm=3*action.F+2*U*action.F**2/Theta**2
    ef=alphas(action)
    if abs(2-ef["alpha_B"])<1e-12:
        raise ValueError("2-alpha_B=0: strong braiding degeneracy")
    qs_eft=2*action.F*ef["D"]/(2-ef["alpha_B"])**2
    return dict(Qs_ADM=qs_adm,Qs_EFT=qs_eft,
          identity_residual=qs_adm-qs_eft,
          alpha_per_zdot=alpha_per_zdot,
          B_per_zdot=B_per_zdot,
          Theta=Theta,U=U,
          scalar_no_ghost=qs_adm>0,
          alphas=ef)

def exact_adm_density_cosine(action,k,x,alpha=0.,B=0.,
                             zeta=0.,zeta_dot=0.):
    """Evaluate exact unitary ADM action at a 1D periodic scalar cosine.

    h_ij=a² exp(2zeta(x)) delta_ij,
    N=1+alpha cos(kx), N_x=partial_x beta.
    B=k² beta_amplitude/a², k>0.
    Average produces physical-volume action density per background a³.

    This provides an independent finite-difference check of Hessian.
    """
    if k<=0:raise ValueError("k>0 required")
    x=np.asarray(x,dtype=float)
    cs=np.cos(k*x);sn=np.sin(k*x)
    zz=zeta*cs;ztime=zeta_dot*cs
    zx=-k*zeta*sn
    zxx=-k*k*zeta*cs
    N=1+alpha*cs
    if np.any(N<=0):raise ValueError("nonpositive lapse")
    bx=-action.a**2*B*sn/k
    bxx=-action.a**2*B*cs
    e=np.exp(-2*zz)/action.a**2
    Kx=(action.H+ztime-e*(bxx-zx*bx))/N
    Ky=(action.H+ztime-e*zx*bx)/N
    K=Kx+2*Ky
    extrinsic=-4*Kx*Ky-2*Ky**2
    R3=e*(-4*zxx-2*zx**2)
    A2,A3=coefficients_at_lapse(action,N)
    lag=N*np.exp(3*zz)*(A2+A3*K+.5*action.F*(R3+extrinsic))
    return float(np.mean(lag))

def numeric_kinetic_hessian(action,k=2.,n=512,eps=1e-4):
    """Finite-difference the EXACT ADM lapse/shift/curvature action.

    Not inferred from the EFT alpha mapping or quadratic expansion.
    """
    xx=np.linspace(0,2*np.pi/k,n,endpoint=False)
    def L(q):
        return exact_adm_density_cosine(action,k,xx,
                     alpha=q[0],B=q[1],zeta_dot=q[2])
    out=np.empty((3,3));zero=np.zeros(3)
    for i in range(3):
        for j in range(3):
            if i==j:
                up=zero.copy();down=zero.copy()
                up[i]=eps;down[i]=-eps
                out[i,i]=(L(up)-2*L(zero)+L(down))/eps**2
            else:
                a=zero.copy();b=zero.copy();c=zero.copy();d=zero.copy()
                a[i]+=eps;a[j]+=eps;b[i]+=eps;b[j]-=eps
                c[i]-=eps;c[j]+=eps;d[i]-=eps;d[j]-=eps
                out[i,j]=(L(a)-L(b)-L(c)+L(d))/(4*eps*eps)
    return out

def isolated_scalar_gradient(action,Hdot,alpha_B_dot=0.):
    """Scalar-only gradient checkpoint from complete Horndeski alpha branch.

    F_s=(1/a)d/dt(a F²/Theta)-F, Theta=F H(1-alphaB/2).
    This equals gradient coefficient for the closed vacuum-scalar
    Horndeski curvature quadratic action. In a baryon+photon+support
    system, matter/structural kinetic mixing must also be included.
    """
    eft=alphas(action)
    ab=eft["alpha_B"]
    p=1-.5*ab
    if abs(p)<1e-12:raise ValueError("Theta=0")
    Fdot=action.F_sigma*action.sigma_dot
    theta=action.F*action.H*p
    theta_dot=((Fdot*action.H+action.F*Hdot)*p
                 -.5*action.F*action.H*alpha_B_dot)
    Fs=(action.F**2/theta)*(action.H+
                     2*Fdot/action.F-theta_dot/theta)-action.F
    qs=kinetic_from_constraints(action)["Qs_ADM"]
    return dict(Fs=Fs,Qs=qs,cs2=Fs/qs if qs else None,
                gradient_positive=Fs>0,kinetic_positive=qs>0,
                caveat="isolated scalar-only, not full baryon/pg spectrum")

def no_slip_gradient_identity(action,Hdot,alpha_M_dot):
    """No-Slip reduction of the isolated scalar gradient:
    Fs = F/(1+alpha_M)*[-Hdot/H²-alpha_M_dot/(H*(1+alpha_M))].
    Requires g Z=-F_sigma/2, alpha_B=-2alpha_M.
    """
    ef=alphas(action)
    am=ef["alpha_M"]
    if abs(ef["no_slip_residual"])>1e-9:
        raise ValueError("input is not on the No-Slip branch")
    if abs(1+am)<1e-12:raise ValueError("No-Slip Theta=0")
    p=1+am
    Fs=action.F/p*(-Hdot/action.H**2-alpha_M_dot/(action.H*p))
    alpha_B_dot=-2*alpha_M_dot
    result=isolated_scalar_gradient(action,Hdot,alpha_B_dot)
    return dict(Fs_NoSlip=Fs,Fs_direct=result["Fs"],
                residual=Fs-result["Fs"],
                cs2_NoSlip=Fs/result["Qs"] if result["Qs"] else None,
                positive_gradient=Fs>0,
                kinetic_positive=result["kinetic_positive"])

def diagnostics():
    action=RollingAction(F=1.05,F_sigma=.12,H=.8,sigma_dot=.4,
          k1=1.7,k2=.3,g=-.75,g_sigma=.15)
    h_analytic=adm_kinetic_hessian(action)["Hessian"]
    h_numeric=numeric_kinetic_hessian(action)
    kinetic=kinetic_from_constraints(action)
    no_slip=RollingAction(F=1.05,F_sigma=.12,H=.8,sigma_dot=.4,
          k1=1.7,k2=.3,g=-.75,g_sigma=.15)
    grad=no_slip_gradient_identity(no_slip,Hdot=-.1,alpha_M_dot=.02)
    return dict(scope="exact accepted rolling G3 ADM kinetic/isolated gradient; no support coupling",
          A2_A3_derivatives=lapse_derivatives(action),
          analytic_Hessian=h_analytic.tolist(),
          numeric_Hessian=h_numeric.tolist(),
          max_Hessian_difference=float(np.max(np.abs(h_numeric-h_analytic))),
          kinetic={k:v for k,v in kinetic.items() if k!="alphas"},
          alphas=kinetic["alphas"],isolated_NoSlip_gradient=grad,
          missing=["full supporting-space ADM Hessian and clock mixing",
                   "physical averaged pg stress","matter perturbations",
                   "background equations, complete cs2, no-CDM CMB/lensing"])
if __name__=="__main__":
    print(json.dumps(diagnostics(),indent=2))
