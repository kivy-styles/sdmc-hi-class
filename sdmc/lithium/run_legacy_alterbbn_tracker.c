/* SDMC investigation L: independent historical AlterBBN nuclear-network check.
 *
 * Clone external AlterBBN source at pinned upstream SHA in CI; do not vendor it.
 * This uses a *radiationlike extra-energy proxy* for the SDMC BBN tracker
 * (ndd=4, rho_dark/rho_standard at T_BBN = fr/(1-fr)), not exact H(T)
 * from the full canonical+Horndeski microscopic action.
 *
 * The primordial *total* eventual lithium is direct Li7 + electron-capture
 * descendant Be7: Y_Li_total = Y(Li7) + Y(Be7).
 *
 * Upstream: Espen S. Jenssen's 2016 AlterBBN_master (2016 frozen legacy
 * network), https://github.com/espensem/AlterBBN_master
 * Upstream pinned SHA 8219dd4b578244e36d56ebda31541705ae382093
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "src/include.h"

int main(int argc, char **argv) {
    if (argc != 4) {
      fprintf(stderr,"USAGE: %s <label> <omega_b_h2> <radiation_tracker_fraction>\n",argv[0]);
      return 2;
    }
    const char *label=argv[1];
    const double omega_b=atof(argv[2]);
    const double f_r=atof(argv[3]);
    if (!(omega_b>0.015 && omega_b<0.03 && f_r>=0. && f_r<0.2)) {
       fprintf(stderr,"BAD_PHYSICAL_INPUT\n"); return 2;
    }
    /* Eta conversion taken from PDG-2025 combined reference omega_b h²=0.02233,
     * eta_10=6.115.  The scalar perturbation has NO independent baryon production.
     * Keep Neff and neutron lifetime equal in paired source trials.
     */
    double eta=6.115e-10*(omega_b/0.02233);
    double kappa=f_r/(1-f_r);
    struct relicparam p;
    Init_cosmomodel(&p);
    Init_cosmomodel_param(27.0,eta,3.046,0.0,879.4,0.,0.,0.,&p);
    if (f_r>0.) Init_dark_density(kappa,4.0,0.0,&p);
    double y[NNUC+1]={0};
    int status=nucl(0,p,y);
    if (status!=0) {
        fprintf(stderr,"BBN_SOLVER_FAILED %s %d\n",label,status);return 1;
    }
    double lithium_total=y[8]+y[9];
    if (!(isfinite(lithium_total) && lithium_total>1e-11 &&
          isfinite(y[6]) && y[6]>0.1 && y[6]<0.5)) {
        fprintf(stderr,"BAD_BBN_RESULT %s Yp=%g Li7=%g Be7=%g\n",
                label,y[6],y[8],y[9]);
        return 1;
    }
    printf("BBN_RESULT label=%s omega_b=%.12g eta10=%.12g fr=%.12g "
           "kappa=%.12g Yp=%.12e D_H=%.12e He3_H=%.12e "
           "Li7_direct=%.12e Be7=%.12e Li7_total=%.12e Hspeed_ideal=%.12e\n",
           label,omega_b,eta*1e10,f_r,kappa,y[6],y[3],y[5],
           y[8],y[9],lithium_total,1/sqrt(1-f_r));
    return 0;
}
