#!/usr/bin/env python3
"""
Structural Planck recovery around the SN-friendly late107 late background.

Low-z geometry is frozen at the SN-friendly point:
  A_late=0.018589897081255913
  B_late=0.013815921754576268
  H0=70.09917331933976

Primordial amplitude uses the tiny Q improvement found in the preceding
ordinary pass:
  Q=2.943463801247999

Reopen only perturbation/kinetic structure:
  A_F, z_c, width, D_floor

Because the explicit sdmc_full background is fixed, these directions are
intended to recover CMB likelihood without undoing the SN distance gain.
"""
from pathlib import Path
import json,math,re,subprocess,textwrap
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc
from cobaya.likelihoods.planck_2018_highl_plik.TTTEEE_lite_native import TTTEEE_lite_native
from cobaya.likelihoods.planck_2018_lowl.TT import TT
from cobaya.likelihoods.planck_2018_lowl.EE import EE
from cobaya.likelihoods.planck_2018_lensing import native as LensingNative

OUT=Path("output/snclosure107_structural"); OUT.mkdir(parents=True,exist_ok=True)
TCMB=2.7255; CAL_SIGMA=.0025; OR=4.17998772e-5

# Frozen SN-friendly background + ordinary sector
H0=70.09917331933976
OB=0.022083219194622913; OC=0.12299536722293603
NS=0.9625227132590487; TAU=0.055202901571989066
Q=2.943463801247999; AS=math.exp(Q+2*TAU)/1e10
AL=0.0013601897284388; BL=0.0147392870020121
LAM=18.40625; ZT=16.173189924377947; DN=.5; TAUA=.25; TAUB=1.5
D0=0.34231919445927034

# edge024 structural center
AF0=0.024221895329654217; ZC0=4.145109392143786
W0=0.3370092310383916; DF0=0.04453643597662449

# Fair-ledger constants used only for promotion ranking.
EDGE_P_LITE=1023.801347
EDGE_PD_FAIR=-0.123919
EDGE_BAO=13.433001
NEW_BAO=13.433001
SN_PP=1.631268
SN_U3=1.395445
SN_D5=3.275618

LOW=np.array([0.022,3.55,0.28,0.044])
HIGH=np.array([0.044,4.40,0.43,0.070])

high=TTTEEE_lite_native(packages_path="planck_packages")
lowT=TT(packages_path="planck_packages"); lowE=EE(packages_path="planck_packages")
lens=LensingNative(packages_path="planck_packages")

def table(path):
    lines=Path(path).read_text().splitlines()
    hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
    ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
    for i,m in enumerate(ms):
        e=ms[i+1].start() if i+1<len(ms) else len(hdr)
        names.append(hdr[m.end():e].strip())
    return pd.DataFrame(np.loadtxt(path),columns=names)

def load_cls(path):
    a=np.loadtxt(path); e=a[:,0].astype(int); n=int(e.max())+1; cv=(TCMB*1e6)**2
    tt=np.zeros(n); ee=np.zeros(n); te=np.zeros(n); pp=np.zeros(n)
    tt[e]=a[:,1]*cv; ee[e]=a[:,2]*cv; te[e]=a[:,3]*cv
    L=e.astype(float); pp[e]=a[:,5]*L*(L+1)
    return np.column_stack([e,tt[e],te[e],ee[e]]),{"tt":tt,"te":te,"ee":ee,"pp":pp}

def pscore(path):
    ha,d=load_cls(path)
    def f(A):
        x=float(high.chi_squared(ha,A_planck=A))
        x+=float(-2*lowT.log_likelihood(d["tt"],calib=A))
        x+=float(-2*lowE.log_likelihood(d["ee"],calib=A))
        lp={lens.calibration_param:A} if getattr(lens,"calibration_param",None) else {}
        x+=float(-2*lens.log_likelihood(d,**lp))
        x+=((A-1)/CAL_SIGMA)**2
        return x
    op=minimize_scalar(f,bounds=(.97,1.03),method="bounded",options={"xatol":1e-9})
    return float(op.fun),float(op.x)

def ini(root,AF,ZC,W,DF):
    h=H0/100.; ox=1.-(OB+OC+OR)/(h*h)
    return textwrap.dedent(f"""\
H0={H0:.17g}
omega_b={OB:.17g}
omega_cdm={OC:.17g}
N_ncdm=0
N_ur=3.046
T_cmb=2.7255
YHe=0.2453
A_s={AS:.17e}
n_s={NS:.17g}
tau_reio={TAU:.17g}
Omega_Lambda=0
Omega_fld=3.1443554e-8
fluid_equation_of_state=SDMC_TRACKER
cs2_fld=0.003
use_ppf=no
Omega_smg=-1
gravity_model=sdmc_v3_independent_kinetic
parameters_smg={AF:.17g},{ZC:.17g},{W:.17g},{D0:.17g},1.0,{DF:.17g}
expansion_model=sdmc_full
expansion_smg={ox:.17g},{LAM},{ZT},{DN},{AL},{TAUA},{BL},{TAUB}
pert_initial_conditions_smg=zero
method_qs_smg=fully_dynamic
a_ini_over_a_today_default=1.e-8
a_ini_test_qs_smg=1.e-8
pert_ic_ini_z_ref_smg=1.e7
a_min_stability_test_smg=1.e-8
output_background_smg=3
modes=s
output=tCl,pCl,lCl
lensing=yes
l_max_scalars=3000
write background=yes
root={root}
format=class
input_verbose=0
background_verbose=0
thermodynamics_verbose=0
perturbations_verbose=0
spectra_verbose=0
lensing_verbose=0
output_verbose=0
""")

def cand(tag,AF,ZC,W,DF):
    root=str(OUT/(tag+"_")); ip=OUT/(tag+".ini"); ip.write_text(ini(root,AF,ZC,W,DF))
    cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=300)
    r={"id":tag,"AF":float(AF),"ZC":float(ZC),"width":float(W),"D_floor":float(DF),"status":"FAIL"}
    bgp=Path(root+"00_background.dat"); clp=Path(root+"00_cl_lensed.dat")
    if cp.returncode==0 and bgp.exists() and clp.exists():
        bg=table(bgp)
        r.update(min_D=float(bg["kin (D)"].min()),
                 min_cs2=float(bg["c_s^2"].min()),
                 max_cs2=float(bg["c_s^2"].max()),
                 max_abs_noslip=float(np.max(np.abs(bg["braiding_smg"]+2*bg["M2_running_smg"]))))
        stable=r["min_D"]>0 and r["min_cs2"]>0 and r["max_cs2"]<=1
        r["stable_subluminal"]=bool(stable)
        if stable:
            pc,Ap=pscore(clp); r.update(status="OK",chi2_planck=pc,A_planck=Ap)
            r["delta_planck_vs_edge024"]=pc-EDGE_P_LITE
            r["pd_proxy"]=EDGE_PD_FAIR+r["delta_planck_vs_edge024"]+.25*(NEW_BAO-EDGE_BAO)
            r["joint_pp_proxy"]=r["pd_proxy"]+SN_PP
            r["joint_u3_proxy"]=r["pd_proxy"]+SN_U3
            r["joint_d5_proxy"]=r["pd_proxy"]+SN_D5
            js=[r["joint_pp_proxy"],r["joint_u3_proxy"],r["joint_d5_proxy"]]
            r["second_joint_proxy"]=sorted(js)[1]
            r["n_sn_closed_proxy"]=sum(x<0 for x in js)
            r["goal_score"]=pc
    if r["status"]!="OK": r["error"]=cp.stdout[-500:].replace("\n"," | ")
    print("SN107S_POINT",json.dumps(r,sort_keys=True),flush=True)
    return r

pts=[("center",AF0,ZC0,W0,DF0)]
# coordinate anchors
pts += [
 ("afm",AF0-.004,ZC0,W0,DF0),("afp",AF0+.004,ZC0,W0,DF0),
 ("zcm",AF0,ZC0-.20,W0,DF0),("zcp",AF0,ZC0+.20,W0,DF0),
 ("wm",AF0,ZC0,W0-.035,DF0),("wp",AF0,ZC0,W0+.035,DF0),
 ("dfm",AF0,ZC0,W0,DF0-.004),("dfp",AF0,ZC0,W0,DF0+.004)]
sam=qmc.Sobol(d=4,scramble=True,seed=22341)
for i,p in enumerate(qmc.scale(sam.random_base2(m=7),LOW,HIGH)):
    pts.append((f"sobol{i:03d}",*map(float,p)))
rows=[cand(*p) for p in pts]
df=pd.DataFrame(rows); df.to_csv(OUT/"snclosure107_structural.csv",index=False)
ok=df[df.status=="OK"].sort_values(["chi2_planck","min_cs2"],ascending=[True,False])
best=ok.head(12).to_dict("records")
summary={"background":{"H0":H0,"A":AL,"B":BL},"Q":Q,"n_ok":int(len(ok)),
         "best":best,"target":"minimize Planck-lite at fixed late107 SN-friendly background; exact raw DESI/full-Plik promotion required"}
(OUT/"snclosure107_structural_summary.json").write_text(json.dumps(summary,indent=2))
print("SN107S_BEST",json.dumps(summary,sort_keys=True),flush=True)
