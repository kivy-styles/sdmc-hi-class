#!/usr/bin/env python3
from pathlib import Path
import argparse, json, math, re, subprocess, sys, textwrap
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import minimize

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--screen",required=True)
    ap.add_argument("--data-dir",default="desi_likelihood")
    ap.add_argument("--desi-repo",default="desi-kp-cosmological-likelihoods")
    ap.add_argument("--out",default="output/stage2")
    args=ap.parse_args()

    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    screen=pd.read_csv(args.screen)
    sys.path.insert(0,str(Path(args.desi_repo)/"dr1/cobaya"))
    import lsstypes as types
    from cosmoprimo import PowerSpectrumInterpolator1D,PowerSpectrumBAOFilter,Cosmology
    from cosmoprimo.fiducial import DESI
    from velocileptors.EPT.ept_fullresum_varyDz_nu_fftw import REPT
    from desi_fs_bao_all import list_zrange,dataset_fn,get_tracer_label,get_physical_stochastic_settings

    blocks=[]
    records=[]
    for tracer,iz,zrange in list_zrange:
        if "lya" in tracer.lower():
            continue
        fn=dataset_fn(args.data_dir,tracer,zrange,observable_name="spectrum-poles-rotated")
        d=types.read(fn); sp=d.observable.get("spectrum")
        th=d.window.theory.get("spectrum")
        kin=np.asarray(th.get(0).coords("k"),float)
        win=d.window.at.observable.get("spectrum").at.theory.get("spectrum").value()
        namespace=f"{get_tracer_label(tracer)}_z{iz}"
        zeff=float(sp.attrs["zeff"])
        blocks.append(dict(tracer=tracer,namespace=namespace,zeff=zeff,
                           data=sp.value(),precision=np.linalg.inv(d.covariance.value()),
                           window=win,kin=kin,
                           shotnoise=float(np.mean(sp.get(0).values("shotnoise")))))
        records.append(dict(tracer=tracer,iz=iz,zrange=list(zrange),namespace=namespace,zeff=zeff))
    (out/"desi_zeff.json").write_text(json.dumps(records,indent=2))
    zvals=[0.0]+sorted({b["zeff"] for b in blocks})
    zpk=",".join(f"{z:.12g}" for z in zvals)
    print("STAGE2_DESI_Z",zvals,flush=True)

    def sdmc_ini(row,root):
        As=math.exp(3.076)/1e10
        return textwrap.dedent(f"""\
        H0 = 70.8514
        omega_b = 0.02239952
        omega_cdm = 0.12444227328918850
        N_ncdm = 0
        N_ur = 3.046
        T_cmb = 2.7255
        YHe = 0.2453
        A_s = {As:.16e}
        n_s = 0.964
        tau_reio = 0.0544
        Omega_Lambda = 0
        Omega_fld = 3.1443554e-8
        fluid_equation_of_state = SDMC_TRACKER
        cs2_fld = 0.003
        use_ppf = no
        Omega_smg = -1
        gravity_model = sdmc_v3_independent_kinetic
        parameters_smg = {float(row.AF):.12g}, {float(row.zc):.12g}, {float(row.width):.12g}, 0.40, 1.0, 0.0001
        expansion_model = sdmc_full
        expansion_smg = 0.7073985893,20.0,21.10,0.5,0.01105624999,0.25,0.01951933685,1.5
        pert_initial_conditions_smg = zero
        method_qs_smg = fully_dynamic
        a_ini_over_a_today_default = 1.e-8
        a_ini_test_qs_smg = 1.e-8
        pert_ic_ini_z_ref_smg = 1.e7
        a_min_stability_test_smg = 1.e-8
        gauge = synchronous
        modes = s
        output = mPk,dTk,vTk
        extra_metric_transfer_functions = yes
        matter_source_in_current_gauge = no
        P_k_max_h/Mpc = 2.0
        z_pk = {zpk}
        output_background_smg = 3
        write background = yes
        write thermodynamics = yes
        root = {root}
        format = class
        input_verbose = 0
        background_verbose = 0
        thermodynamics_verbose = 0
        perturbations_verbose = 0
        spectra_verbose = 0
        output_verbose = 0
        """)

    def lcdm_ini(root):
        return textwrap.dedent(f"""\
        H0 = 67.36
        omega_b = 0.02237
        omega_cdm = 0.1200
        N_ncdm = 0
        N_ur = 3.046
        T_cmb = 2.7255
        YHe = 0.2453
        A_s = 2.10e-9
        n_s = 0.9649
        tau_reio = 0.0544
        gauge = synchronous
        modes = s
        output = mPk,dTk,vTk
        extra_metric_transfer_functions = yes
        matter_source_in_current_gauge = no
        P_k_max_h/Mpc = 2.0
        z_pk = {zpk}
        write background = yes
        write thermodynamics = yes
        root = {root}
        format = class
        input_verbose = 0
        background_verbose = 0
        thermodynamics_verbose = 0
        perturbations_verbose = 0
        spectra_verbose = 0
        output_verbose = 0
        """)

    def run_ini(tag,text):
        ip=out/f"{tag}.ini"; ip.write_text(text)
        cp=subprocess.run(["./class",str(ip)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                          text=True,timeout=360)
        (out/f"{tag}.log").write_text(cp.stdout)
        print("STAGE2_TRANSFER",tag,"returncode",cp.returncode,flush=True)
        if cp.returncode:
            raise RuntimeError(f"{tag} failed\n{cp.stdout[-2500:]}")

    for row in screen.itertuples(index=False):
        run_ini(str(row.id),sdmc_ini(row,str(out/(str(row.id)+"_"))))
    run_ini("lcdm",lcdm_ini(str(out/"lcdm_")))

    C=299792.458
    z_to_i={round(z,8):i+1 for i,z in enumerate(zvals)}

    def tab(path):
        lines=Path(path).read_text().splitlines()
        hdr=[l for l in lines if l.startswith("#") and re.search(r"1\s*:",l)][-1].lstrip("#").strip()
        ms=list(re.finditer(r"(\d+)\s*:\s*",hdr)); names=[]
        for i,m in enumerate(ms):
            e=ms[i+1].start() if i+1<len(ms) else len(hdr)
            names.append(hdr[m.end():e].strip())
        df=pd.DataFrame(np.loadtxt(path),columns=names)
        return df.sort_values("z") if "z" in df.columns else df

    def rd_from(prefix,bg):
        th=tab(out/f"{prefix}_00_thermodynamics.dat").sort_values("z")
        zz=th.z.to_numpy(); kb=th.kappa_b.to_numpy(); cross=[]
        for i in range(len(zz)-1):
            if (kb[i]-1)*(kb[i+1]-1)<=0 and kb[i]!=kb[i+1]:
                q=(1-kb[i])/(kb[i+1]-kb[i]); cross.append(zz[i]+q*(zz[i+1]-zz[i]))
        if not cross:
            raise RuntimeError("no drag crossing for "+prefix)
        zd=float(cross[0])
        return float(np.interp(zd,bg.z,bg["comov.snd.hrz."]))

    def model(prefix,H0,ob,oc,ns):
        h=H0/100.; fb=ob/(ob+oc); fc=1-fb
        bg=tab(out/f"{prefix}_00_background.dat").sort_values("z")
        rd=rd_from(prefix,bg); series={}
        for z in zvals:
            ii=z_to_i[round(z,8)]
            pk=np.loadtxt(out/f"{prefix}_00_z{ii}_pk.dat")
            tk=tab(out/f"{prefix}_00_z{ii}_tk.dat")
            kh=tk["k (h/Mpc)"].to_numpy()
            dm=tk["d_m"].to_numpy()
            tb=tk["t_b"].to_numpy()
            hp=tk["h_prime"].to_numpy(); ep=tk["eta_prime"].to_numpy()
            shift=.5*(hp+6*ep)
            a=1/(1+z); H=float(np.interp(z,bg.z,bg["H [1/Mpc]"])); Hc=a*H
            vb=-(tb+shift)/Hc
            vc=-shift/Hc
            vcb=fb*vb+fc*vc
            ratio=np.divide(vcb,dm,out=np.zeros_like(vcb),where=np.abs(dm)>1e-300)
            kpk=pk[:,0]; Pdd=pk[:,1]; rv=np.interp(kpk,kh,ratio)
            series[z]=(kpk,Pdd,Pdd*rv**2)
        return dict(prefix=prefix,H0=H0,h=h,ob=ob,oc=oc,ns=ns,bg=bg,rd=rd,series=series)

    kobs=blocks[0]["kin"]
    fid=DESI(engine="camb")
    starts={"BGS_z0":(1.11348,0.660148,-0.223088),
            "LRG_z0":(1.1351,-0.139101,-0.910397),
            "LRG_z1":(1.22446,-0.671727,-0.22624),
            "LRG_z2":(1.04564,-0.493139,-0.24021),
            "ELG_z1":(0.127936,-0.578801,-1.45356),
            "QSO_z0":(0.780983,0.353843,0.135725)}
    seeds_extra=[np.array([1.,0.,0.]),np.array([.55,4.,-4.]),
                 np.array([1.8,-4.,4.]),np.array([2.6,8.,-8.]),
                 np.array([1.4,-10.,-10.])]
    all_marg=["alpha0p","alpha2p","alpha4p","alpha6p","sn0p","sn2p","sn4p"]
    scales=np.array([12.5]*4+[2.]+[5.]*2)
    marg=["alpha0p","alpha2p","sn0p","sn2p"]
    gi=np.array([all_marg.index(x) for x in marg])
    prior_hess=-np.diag(scales[gi]**-2)

    def prepare(m):
        kin=np.geomspace(min(5e-4,kobs[0]/2),max(1.0,kobs[-1]*2),500)
        Pdd=np.stack([np.interp(kin,*m["series"][z][:2]) for z in zvals],axis=-1)
        Ptt=np.stack([np.interp(kin,m["series"][z][0],m["series"][z][2]) for z in zvals],axis=-1)
        pki=PowerSpectrumInterpolator1D(kin,Pdd)
        cos=Cosmology(n_s=m["ns"],Omega_b=m["ob"]/m["h"]**2,
                      Omega_cdm=m["oc"]/m["h"]**2,Omega_ncdm=0.,H0=m["H0"])
        cos.rs_drag=m["rd"]*m["h"]
        filt=PowerSpectrumBAOFilter(pki,engine="peakaverage",cosmo=cos,cosmo_fid=fid)
        filt(pki,cosmo=cos); Pnw=filt.smooth_pk_interpolator()(kin)
        sig8=np.asarray(pki.sigma8())
        fs8=np.asarray(PowerSpectrumInterpolator1D(kin,Ptt).sigma8())
        pt=REPT(kin,Pdd[:,0],pnw=Pnw[:,0],kmin=kobs[0],kmax=kobs[-1],nk=200,
                rbao=110,sbao=None,beyond_gauss=True,one_loop=True,shear=True,
                cutoff=20,jn=5,N=4000,threads=2,extrap_min=-4,extrap_max=3,import_wisdom=False)
        extk=np.append(pt.kv,1.)
        def rg(arr):
            return 10**interp1d(np.log10(kin),np.log10(np.maximum(arr,1e-300)),
                                kind="cubic",fill_value="extrapolate",axis=0,
                                assume_sorted=True)(np.log10(extk))
        pcb,pnw,ptt=rg(Pdd),rg(Pnw),rg(Ptt)
        ans={}
        for b in blocks:
            z=b["zeff"]; iz=zvals.index(z)
            Hkms=float(np.interp(z,m["bg"].z,m["bg"]["H [1/Mpc]"]))*C
            DA=float(np.interp(z,m["bg"].z,m["bg"]["comov. dist."]))/(1+z)
            qpar=float(fid.efunc(z)/(Hkms/(100*m["h"])))
            qper=float(DA*m["h"]/fid.angular_diameter_distance(z))
            Dz=float(np.sqrt(pcb[-1,iz]/pcb[-1,0]))
            fk=np.sqrt(np.maximum(ptt[:-1,iz]/pcb[:-1,iz],0))
            pks=pt.compute_redshift_space_power_multipoles_tables(
                fk,apar=qpar,aperp=qper,ngauss=4,
                pcb=pcb[:-1,iz],pcb_nw=pnw[:-1,iz],Dz=Dz)[1:]
            basis=np.stack([interp1d(pt.kv,pks[j],kind="cubic",fill_value="extrapolate",
                                    axis=0,assume_sorted=True)(kobs) for j in range(3)],axis=0)
            ans[b["namespace"]]=dict(basis=basis,sigma8=float(sig8[iz]),fsigma8=float(fs8[iz]),
                                     qpar=qpar,qper=qper,Dz=Dz)
        return ans

    def poles(info,x,settings,sn):
        pkt=info["basis"]; sigma8=info["sigma8"]; f=info["fsigma8"]/sigma8
        b1p,b2p,bsp=x
        b1L=b1p/sigma8-1.; b2L=b2p/sigma8**2; bsL=bsp/sigma8**2
        b1=1+b1L; b2=8/21*b1L+b2L; bs=bsL-2/7*b1L; b3=b1L
        G=np.zeros((7,7))
        G[0,0]=b1*b1; G[1,0]=G[1,1]=f*b1
        G[2,1]=f*f; G[2,2]=f*b1; G[3,2]=f*f
        for jj,ell in enumerate([0,2,4]):
            G[jj-3,jj-3]=sn*(settings["fsat"] if ell>0 else 1.)*settings["sigv"]**ell
        base=[1,b1,b1*b1,b2,b1*b2,b2*b2,bs,b1*bs,b2*bs,bs*bs,b3,b1*b3]
        pp=np.sum(pkt*np.concatenate([np.asarray(base),G.dot(np.zeros(7))]),axis=-1)
        return pp,pkt[...,-7:].dot(G)

    def profile(label,m):
        tables=prepare(m); total=0.; rows=[]
        for b in blocks:
            settings=get_physical_stochastic_settings(tracer=b["tracer"].upper()[:3])
            def logp(x):
                pp,g=poles(tables[b["namespace"]],x,settings,b["shotnoise"])
                theory=b["window"].dot(pp.ravel())
                g=b["window"].dot(g[...,gi].reshape(-1,len(gi)))
                diff=theory-b["data"]; pg=b["precision"].dot(g)
                postgrad=-pg.T.dot(diff); lh=-g.T.dot(pg); ph=prior_hess+lh
                dx=-np.linalg.solve(ph,postgrad)
                val=-.5*diff.T.dot(b["precision"]).dot(diff)
                val+=.5*dx.dot(lh).dot(dx)+postgrad.dot(dx)+.5*dx.dot(prior_hess).dot(dx)
                val+=-.5*np.linalg.slogdet(-ph)[1]-.5*(x[1]/5.)**2-.5*(x[2]/5.)**2
                return float(val)
            seeds=[np.array(starts[b["namespace"]],float)]+[x.copy() for x in seeds_extra]
            trials=[]
            for iseed,seed in enumerate(seeds):
                opt=minimize(lambda x:-2*logp(x),seed,method="L-BFGS-B",
                             bounds=[(0,3),(-20,20),(-20,20)],
                             options={"maxiter":800,"ftol":1e-11,"gtol":3e-7,"maxls":50})
                trials.append(opt)
                print("STAGE2_TRIAL",label,b["namespace"],iseed,float(opt.fun),
                      bool(opt.success),[float(v) for v in opt.x],flush=True)
            finite=[o for o in trials if np.isfinite(o.fun)]
            opt=min(finite,key=lambda o:o.fun)
            total+=float(opt.fun)
            info=tables[b["namespace"]]
            rows.append(dict(model=label,namespace=b["namespace"],chi2=float(opt.fun),
                             b1p=float(opt.x[0]),b2p=float(opt.x[1]),bsp=float(opt.x[2]),
                             qpar=info["qpar"],qper=info["qper"],sigma8=info["sigma8"],
                             fsigma8=info["fsigma8"]))
        print("STAGE2_DESI_TOTAL",label,total,flush=True)
        return rows,total

    lcdm=model("lcdm",67.36,0.02237,0.1200,0.9649)
    lcdm_rows,lcdm_chi=profile("LCDM",lcdm)

    allrows=list(lcdm_rows); summary=[]
    for sr in screen.itertuples(index=False):
        label=str(sr.id)
        m=model(label,70.8514,0.02239952,0.12444227328918850,0.964)
        rows,chi=profile(label,m); allrows.extend(rows)
        dd=chi-lcdm_chi
        rec=dict(id=label,AF=float(sr.AF),zc=float(sr.zc),width=float(sr.width),
                 desi_chi2=chi,lcdm_desi_chi2=lcdm_chi,delta_desi=dd,
                 delta_planck=float(sr.delta_planck),
                 delta_pantheonplus=float(sr.delta_pantheonplus),
                 delta_union3=float(sr.delta_union3),
                 delta_desy5=float(sr.delta_desy5),
                 joint_planck_desi=float(sr.delta_planck)+dd,
                 joint_pp=float(sr.delta_planck)+dd+float(sr.delta_pantheonplus),
                 joint_union3=float(sr.delta_planck)+dd+float(sr.delta_union3),
                 joint_desy5=float(sr.delta_planck)+dd+float(sr.delta_desy5))
        summary.append(rec)
        print("STAGE2_JOINT_RESULT",rec,flush=True)

    pd.DataFrame(allrows).to_csv(out/"stage2_desi_blocks.csv",index=False)
    df=pd.DataFrame(summary).sort_values("joint_pp")
    df.to_csv(out/"stage2_joint_summary.csv",index=False)
    print("STAGE2_BEST_PP",df.nsmallest(4,"joint_pp").to_dict("records"),flush=True)
    print("STAGE2_BEST_UNION3",df.nsmallest(4,"joint_union3").to_dict("records"),flush=True)
    print("STAGE2_BEST_DESY5",df.nsmallest(4,"joint_desy5").to_dict("records"),flush=True)

if __name__=="__main__":
    main()
