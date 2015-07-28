import sys
sys.path.append(".")
sys.path.append("../")
import libmpdata
import numpy as np
import libcloudphxx as libcl
#import analytic_blk_1m_pytest as an
#from constants_pytest import Rd, Rv, cp, p0
import pdb

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import pdb

def test():
    a = np.zeros((10,))
    a[1]=1
    print a
#.5 - l. courant, 3-liczba krokow                                                              
    libmpdata.mpdata(a, 1., 3);
#assert a[1] < 1 and a[2] > 0 #SA - po co ten assert?                                          
    print a


#test()

def plotting(dct, time = None, figname="plot_test.pdf", ylim_dic = {}):
    nrow = (len(dct)+1)/2
    fig, tpl = plt.subplots(nrows=nrow, ncols=2, figsize=(10,8.5))
    i=0
    for k,v in dct.iteritems():
      tpl[i%nrow,i/nrow].set_title(k+", " + time)
      if k in ylim_dic.keys():
          tpl[i%nrow,i/nrow].set_ylim(ylim_dic[k])
      tpl[i%nrow,i/nrow].plot(v)
      i+=1
    plt.savefig(figname)
    plt.show()


def libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt, aerosol):
    opts = libcl.blk_2m.opts_t()
    opts.acti = True
    opts.cond = True
    opts.acnv = False
    opts.accr = False
    opts.sedi = False
    distr = [{
      "mean_rd" : aerosol["meanr"], 
      "sdev_rd" : aerosol["gstdv"], 
      "N_stp"   : aerosol["n_tot"], 
      "chem_b"  : aerosol["chem_b"]
    }]
    opts.dry_distros = distr

    shp = th_d.shape
    dot_th = np.zeros(shp)
    dot_rv = np.zeros(shp)
    dot_rc = np.zeros(shp)
    dot_nc = np.zeros(shp)
    dot_rr = np.zeros(shp)
    dot_nr = np.zeros(shp)


    print "qc min, max przed mikro", rc.min(), rc.max()
    print "nc min, max przed mikro", nc.min(), nc.max()
    
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, th_d, rv, rc, nc, rr, nr, dt)
    
    print "rc min po mikro" ,(rc + dot_rc * dt).min(), (rc + dot_rc * dt).max()
    print "nc min po mikro" ,(nc + dot_nc * dt).min(), (nc + dot_nc * dt).max()

    th_d += dot_th * dt
    rv   += dot_rv * dt
    rc   += dot_rc * dt
    nc   += dot_nc * dt
    rr   += dot_rr * dt
    nr   += dot_nr * dt
    np.place(rc, rc<0, 0)
    np.place(nc, nc<0, 0)
    print "rc min, max po mikro i place", rc.min(), rc.max()
    print "nc min, max po mikro i place", nc.min(), nc.max()

def libcl_1mom(rho_d, th_d, rv, rc, rr, dt):
    opts = libcl.blk_1m.opts_t()
    opts.cond = True
    opts.cevp = True
    opts.revp = False
    opts.conv = False
    opts.accr = False
    opts.sedi = False

    print "1m rc max, min przed mikro", rc.max(), rc.min()
    libcl.blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    print "1m rc max, min po mikro", rc.max(), rc.min()

def libcl_spdr(rho_d, th_d, rv, dt, aerosol):
    opts_init = libcl.lgrngn.opts_init_t()
#    opts_init.dry_distros = {kappa : distros.lognormal(n_tot, meanr, gstdv)}
#    opts_init.dt = dt

#    backend = libcl.lgrngn.backend_t.serial
#    prtcls = libcl.lgrngn.factory(backend, opts_init)

#    opts = libcl.lgrngn.opts_t()
#    opts.adve = False
#    opts.sedi = False
#    opts.cond = True
#    opts.coal = False
#    opts.sstp_cond = sstp
    pass


def calc_RH(RH, Temp, rho_d, th_d, rv):
    for i in range(len(RH)):
	Temp[i] = libcl.common.T(th_d[i], rho_d[i])
	p = libcl.common.p(rho_d[i], rv[i], Temp[i])
        #pdb.set_trace()
	p_v = rho_d[i] * rv[i] * libcl.common.R_v * Temp[i]
	RH[i] = p_v / libcl.common.p_vs(Temp[i])

def main(scheme, 
  nx=300, sl_sg = slice(50,100), crnt=0.1, dt=0.2, nt=1501, outfreq=1500,
  aerosol={
    "meanr":.02e-6, "gstdv":1.4, "n_tot":60e6, 
    # ammonium sulphate aerosol parameters:
    "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
    "kappa":.61    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
  }
):
    th_d = np.ones((nx,))* 303.
    #th_d[sl_sg] += 1
    rv = np.ones((nx,))* 5.e-3
    rc = np.zeros((nx,))
    rv[sl_sg] += 8.e-3
    #rc[sl_sg] += 1.e-3 
    rr = np.zeros((nx,))
    rho_d = np.ones((nx,))*.97
    testowa = np.zeros((nx,))
    testowa[sl_sg]= 1.e4

    RH1 = np.empty((nx,))
    RH2 = np.empty((nx,))
    Temp = np.empty((nx,))
    
    var_adv = [th_d, rv, testowa]
 
    if scheme in ["1m", "2m"]:
           rc = np.zeros((nx,))
           rr = np.zeros((nx,))
           var_adv = var_adv + [rc, rr]

    if scheme == "2m":
           nc = np.zeros((nx,))
           #nc[sl_sg] = 3.e7
           nr = np.zeros((nx,))
           var_adv  = var_adv + [nc, nr]

    calc_RH(RH2, Temp, rho_d, th_d, rv)
    dic_var = {"rc":rc, "rv":rv, "th":th_d, "Temp":Temp, "S":RH2-1}
    if   scheme == "2m":
        dic_var["nc"] =nc

    plotting(dic_var, figname="plot_init.pdf", time="init") 
    for it in range(nt):
        print "it", it

        print "testowa min, max przed adv", testowa.min(), testowa.max()
        if scheme == "2m": print "qc min, max przed adv", rc.min(), rc.max()        
        for var in var_adv:
            libmpdata.mpdata(var, crnt, 1);
        if scheme == "2m": print "qc min, max po adv", rc.min(), rc.max()
        print "testowa min, max po adv", testowa.min(), testowa.max()

        #calc_RH(RH1, Temp, rho_d, th_d, rv)

        if   scheme == "1m":
            libcl_1mom(rho_d, th_d, rv, rc, rr,         dt, aerosol)
        elif scheme == "2m":
            libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt, aerosol)
        elif scheme == "sd":
            libcl_spdr(rho_d, th_d, rv,                 dt, aerosol)
        else: 
            assert(False)

        calc_RH(RH2, Temp, rho_d, th_d, rv) 
                
        print "testowa po it = ", it
        if it % outfreq == 0 or it in [100]:
            dic_var = {"rc":rc, "rv":rv, "th":th_d, "Temp":Temp, "S":RH2-1}
            if   scheme == "2m":
                dic_var["nc"] = nc
            plotting(dic_var, figname="plot_"+str(int(it*dt))+"s.pdf", 
            time=str(int(it*dt))+"s" )
            plotting(dic_var, figname="plot_"+str(int(it*dt))+"s_ylim.pdf",
                     time=str(int(it*dt))+"s", ylim_dic={"S":[-0.005, 0.015], "nc":[4.86e7, 4.92e7], "rv":[0.0119,0.0121], "rc":[0.00098, 0.00104]} )



#main("2m") 
#main("1m") 
main("sd")
