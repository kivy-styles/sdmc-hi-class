#!/usr/bin/env python3
"""Autonomous curvature/activation handoff stress test.

Not the reconstructed late300 action. Every run has a *fixed*, analytic
Jordan-frame canonical scalar-tensor action:
    G2=X-V(phi), G3=0, G4=F(phi)/2, G5=0,
    F=exp(beta*phi),
    V=Vs exp(-lambda_e*phi)+Vir exp(-sqrt(2)*phi).
rho_m and rho_r conserve separately. No H(z), F(N), z_t or activation
window is imposed. All action amplitudes and initial data are disclosed;
the model has NOT passed Planck, DESI, stability or No-Slip tests.

Uses 8 pi G_ref=1 and a reference H_ref=1; rho=3 Omega_ref H_ref^2.
"""
import json
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LAMBDA_E = 18.40625       # borrowed from benchmark, not first-principles derived
LAMBDA_IR = math.sqrt(2)   # manuscript canonical mature-branch slope
OMEGA_M_REF = 0.2985059603044
OMEGA_R_REF = 4.17998772e-5 / (0.6971482083085**2)
V_IR = 3.0*0.7014         # chosen reference potential scale, NOT derived
N_START = -12.0
N_END = 1.5

def evolve(beta, steep_ratio, phi_ini=-0.5, y_ini=0.0):
    V_STEEP = V_IR*steep_ratio
    def point(n, phi, y):
        F = math.exp(beta*phi)
        Fp = beta*F
        Fpp = beta*beta*F
        steep = V_STEEP*math.exp(-LAMBDA_E*phi)
        mature = V_IR*math.exp(-LAMBDA_IR*phi)
        V = steep + mature
        Vp = -LAMBDA_E*steep - LAMBDA_IR*mature
        rm = 3*OMEGA_M_REF*math.exp(-3*n)
        rr = 3*OMEGA_R_REF*math.exp(-4*n)
        denom = 3*F+3*Fp*y-y*y/2
        if denom <= 0:
            raise ValueError("Unhealthy Hamiltonian denominator")
        H2 = (rm+rr+V)/denom
        # Scalar: y_N+(3+u)y+V_phi/H^2-3F_phi(2+u)=0
        rhs_s = -3*y - Vp/H2 + 6*Fp
        # Raychaudhuri: -2 F u = (rm+4rr/3)/H2+y^2+F_NN+(u-1)F_N
        rhs_r = -(rm+4*rr/3)/H2 - (1+Fpp)*y*y + Fp*y
        u = (rhs_r-Fp*rhs_s)/(2*F+3*Fp*Fp)
        y_N = rhs_s-(y-3*Fp)*u
        rhoX = 3*F*H2-rm-rr
        return dict(phi=phi, y=y, y_N=y_N, F=F,
                    alphaM=beta*y, H2=H2, u=u, q=-1-u,
                    QJ=rhoX/rm, rhoX=rhoX, steep_share=steep/V)
    def dy(n, state):
        d = point(n, float(state[0]), float(state[1]))
        return [state[1], d["y_N"]]
    sol = solve_ivp(dy,(N_START,N_END),[phi_ini,y_ini],
                    dense_output=True,method="DOP853",rtol=2e-10,
                    atol=3e-12,max_step=0.035)
    if not sol.success: raise RuntimeError(sol.message)
    def at_n(n):
        phi, y = map(float, sol.sol(n))
        return point(n,phi,y)
    grid = np.linspace(-6,0,1001)
    transitions = {}
    for key, target in (("q_equals_zero",("q",0)),
                        ("QJ_equals_one",("QJ",1)),
                        ("steep_mature_equal",("steep_share",0.5))):
        attr, val = target
        vv = [at_n(float(n))[attr]-val for n in grid]
        changes = [i for i in range(len(grid)-1) if vv[i]*vv[i+1]<0]
        roots = [brentq(lambda n:at_n(n)[attr]-val,grid[i],grid[i+1])
                 for i in changes]
        transitions[key] = [math.expm1(-r) for r in roots]
    def sample(z):
        d = at_n(-math.log1p(z))
        return {x:float(d[x]) for x in
                ("phi","y","F","alphaM","H2","q","QJ","steep_share")}
    bianchi_max = 0.0
    for n in np.linspace(-8.95,-.05,56):
        eps=1.e-4
        h1=at_n(n+eps)["H2"]; h0=at_n(n-eps)["H2"]
        numerical=(math.log(h1)-math.log(h0))/(2*eps)
        bianchi_max=max(bianchi_max,abs(numerical-2*at_n(n)["u"]))
    return dict(action_coefficients=dict(beta=beta,steep_ratio=steep_ratio,
                lambda_early=LAMBDA_E,lambda_mature=LAMBDA_IR,
                V_mature=V_IR,V_early=V_STEEP),
                initial=dict(N=N_START,phi=phi_ini,phi_N=y_ini),
                radiation=sample(1090),
                F_window=sample(3.927876388467848),
                present=sample(0),
                transitions_redshift=transitions,
                max_abs_dlnH2_dN_minus_2H_N_over_H=bianchi_max)

def main():
    runs=[evolve(b,r) for b in
          (0.0,0.02048146490100771,-0.02048146490100771)
          for r in (0.01,1.0)]
    maxcheck=max(x["max_abs_dlnH2_dN_minus_2H_N_over_H"]
                 for x in runs)
    assert maxcheck < 1e-5
    output=dict(status="autonomous restricted-model numerical countercheck",
        passed_background_constraint_derivative_check=True,
        max_Bianchi_difference=maxcheck,
        caveat="Not accepted late300, not a unique microphysical SDMC action. "
        "G3=0 and time-varying F generally do not enforce the accepted "
        "No-Slip relation alpha_B=-2 alpha_M. phi was not proven identical "
        "to the SDMC structural clock ln(S/S0). Absolute potential scale "
        "and early initial state remain externally chosen.",runs=runs)
    print("AUTONOMOUS_ACTIVATION_TEST",json.dumps(output,sort_keys=True))
if __name__=="__main__":
    main()
