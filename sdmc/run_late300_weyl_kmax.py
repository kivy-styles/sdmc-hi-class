#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,csv,numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import jv
S=importlib.util.spec_from_file_location("b","sdmc/run_late300_desy3_weyl_tomography.py")
b=importlib.util.module_from_spec(S); S.loader.exec_module(b)
DATA=Path("des_y3_data/likelihood/des-y3/2pt_NG_final_2ptunblind_02_24_21_wnz_covupdate.v2.fits")
CUTS=Path("des_y3_data/examples/des-y3-scale-cuts.ini")
hL=0.6971482083084993; OmL=(0.022083219194622913+0.12299536722293603)/hL**2
hC=0.6856859; OmC=(0.02240637+0.11824151)/hC**2
late=b.load_model("late300","output/weyl_spectra",hL,OmL)
lcdm=b.load_model("local021","output/weyl_spectra",hC,OmC)
des=b.load_des(DATA,CUTS,False)

def theory(m,kmaxh):
 zmax=min(m["zq"].max(),m["zp"].max(),float(np.max(des["zsrc"])))
 z,H,chi,nz,g,nchi=b.prepare_sources(des,m["bg"],zmax)
 zm=np.broadcast_to(z[None,:],(len(b.ELL),len(z)))
 km=(b.ELL[:,None]+0.5)/np.maximum(chi[None,:],1e-12)
 Q=b.eval_cube(m["Iq"],zm,km); P=b.eval_cube(m["Ip"],zm,km)
 ok=np.isfinite(Q)&np.isfinite(P)&(km<=kmaxh*m["h"])
 Q=np.where(ok,Q,0.0); P=np.where(ok,P,0.0)
 pref=np.exp(m["Ip"](np.column_stack([z,np.full_like(z,np.log(0.05))])))
 p0=float(np.exp(m["Ip"]([[0.0,np.log(0.05)]])[0]))
 F=-b.C1RHO*m["Om"]/np.maximum(np.sqrt(pref/p0),1e-4)
 X=np.sqrt(np.maximum(Q*P,0.0)); pair={}
 for i in range(4):
  for j in range(i,4):
   gg=np.trapezoid(g[i][None,:]*g[j][None,:]*Q,x=chi,axis=1)
   gi=np.trapezoid(((g[i]*nchi[j]+g[j]*nchi[i])*F/np.maximum(chi,1e-12))[None,:]*X,x=chi,axis=1)
   ii=np.trapezoid((nchi[i]*nchi[j]*F**2/np.maximum(chi,1e-12)**2)[None,:]*P,x=chi,axis=1)
   pair[(i+1,j+1)]=(gg,gi,ii)
 out=[[],[],[]]
 for r in des["rows"]:
  cls=pair[(min(r["b1"],r["b2"]),max(r["b1"],r["b2"]))]
  th=np.deg2rad(r["ang"]/60.); wt=b.ELL/(2*np.pi)*jv(0 if r["typ"]=="xip" else 4,b.ELL*th)
  for q,cl in enumerate(cls): out[q].append(np.trapezoid(wt*cl,x=b.ELL))
 return tuple(np.asarray(x) for x in out),float(np.mean(ok))

def prof(T):
 L=np.linalg.cholesky(des["cov"]); d=des["data"]
 def f(A):
  y=np.linalg.solve(L,d-(T[0]+A*T[1]+A*A*T[2])); return float(y@y)
 r=minimize_scalar(f,bounds=(-5,5),method="bounded",options={"xatol":1e-4})
 return float(r.fun),float(r.x)

rows=[]; out={}
for km in [0.08,0.10,0.15,0.20,0.30,0.50,1.0]:
 rec={}
 for m in [late,lcdm]:
  T,frac=theory(m,km); chi,A=prof(T)
  rec[m["name"]]={"chi2":chi,"AIA":A,"grid_fraction":frac}
 rec["delta_chi2"]=rec["late300"]["chi2"]-rec["local021"]["chi2"]; out[str(km)]=rec
 rows.append([km,rec["late300"]["chi2"],rec["local021"]["chi2"],rec["delta_chi2"],rec["late300"]["AIA"],rec["local021"]["AIA"]])
 print("WEYL_KMAX",km,json.dumps(rec,sort_keys=True))
R=Path("output/weyl_kmax"); R.mkdir(parents=True,exist_ok=True)
with (R/"late300_weyl_kmax.csv").open("w",newline="") as f:
 w=csv.writer(f); w.writerow(["kmax_h_Mpc","chi2_late300","chi2_local021","delta_chi2","AIA_late300","AIA_local021"]); w.writerows(rows)
(R/"late300_weyl_kmax.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
