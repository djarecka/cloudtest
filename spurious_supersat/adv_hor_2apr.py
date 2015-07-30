import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math
#import analytic_blk_1m_pytest as an
#from constants_pytest import Rd, Rv, cp, p0

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

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


def libcl_2mom(rho_d, thd, rv, rc, rr, nc, nr, dt, aerosol):
    opts = libcl.blk_2m.opts_t()
    opts.acti = opts.cond = True
    opts.acnv = opts.accr = opts.sedi = False
    distr = [{
      "mean_rd" : aerosol["meanr"], 
      "sdev_rd" : aerosol["gstdv"], 
      "N_stp"   : aerosol["n_tot"], 
      "chem_b"  : aerosol["chem_b"]
    }]
    opts.dry_distros = distr

    shp = thd.shape
    dot_th, dot_rv = np.zeros(shp), np.zeros(shp)
    dot_rc, dot_nc = np.zeros(shp), np.zeros(shp)
    dot_rr, dot_nr = np.zeros(shp), np.zeros(shp)

    print "qv min, max przed mikro", rv.min(), rv.max()
    print "qc min, max przed mikro", rc.min(), rc.max()
    print "nc min, max przed mikro", nc.min(), nc.max()
    
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, thd, rv, rc, nc, rr, nr, dt)
    
    print "qv min, max po mikro", rv.min(), rv.max()
    print "rc min po mikro" ,(rc + dot_rc * dt).min(), (rc + dot_rc * dt).max()
    print "nc min po mikro" ,(nc + dot_nc * dt).min(), (nc + dot_nc * dt).max()


    thd  += dot_th * dt
    rv   += dot_rv * dt
    rc   += dot_rc * dt
    nc   += dot_nc * dt
    rr   += dot_rr * dt
    nr   += dot_nr * dt
    #np.place(rv, rv<0, 0)
    np.place(rc, rc<0, 0)
    np.place(nc, nc<0, 0)
    print "rc min, max po mikro i place", rc.min(), rc.max()
    print "nc min, max po mikro i place", nc.min(), nc.max()

def libcl_spdr_init(rho_d, th_d, rv, C, dt, aerosol, dx=2): #TODO dx
    opts_init = libcl.lgrngn.opts_init_t()
    opts_init.dt = dt
    opts_init.nx = rv.shape[0]
    opts_init.dx = dx
    opts_init.x1 = opts_init.nx * opts_init.dx

    def lognormal(lnr):
        from math import exp, log, sqrt, pi
        return aerosol["n_tot"] * exp(
      -(lnr - log(aerosol["meanr"]))**2 / 2 / log(aerosol["gstdv"])**2
    ) / log(aerosol["gstdv"]) / sqrt(2*pi);

    opts_init.sd_conc_mean = aerosol["sd_conc"]
    opts_init.dry_distros = {aerosol["kappa"]:lognormal}

    opts_init.coal_switch = opts_init.sedi_switch = False

    micro = libcl.lgrngn.factory(libcl.lgrngn.backend_t.serial, opts_init)

    C_arr = np.ones(rv.shape[0]+1) * C
    micro.init(th_d, rv, rho_d, C_arr)
    return micro

def libcl_spdr(rhod, thd, rv, rc, nc, na, sd, dt, aerosol, micro):
    libopts = libcl.lgrngn.opts_t()
    libopts.cond = True
    libopts.adve = True
    libopts.coal = libopts.sedi = False

    micro.step_sync(libopts, thd, rv, rhod)
    micro.step_async(libopts)

    # absolute number of super-droplets per grid cell
    micro.diag_sd_conc()
    sd[:] = np.frombuffer(micro.outbuf())
    
    # number of particles (per kg of dry air) with r_w < .5 um
    micro.diag_wet_rng(0, .5e-6)
    micro.diag_wet_mom(0)
    na[:] = np.frombuffer(micro.outbuf())

    # number of particles (per kg of dry air) with r_w > .5 um 
    micro.diag_wet_rng(.5e-6, 1)
    micro.diag_wet_mom(0)
    nc[:] = np.frombuffer(micro.outbuf())

    # cloud water mixing ratio [kg/kg] (same size threshold as above)
    micro.diag_wet_mom(3)
    rho_H2O = 1e3
    rc[:] = 4./3 * math.pi * rho_H2O * np.frombuffer(micro.outbuf())

def calc_S(S, Temp, rho_d, th_d, rv):
    for i in range(len(S)):
	Temp[i] = libcl.common.T(th_d[i], rho_d[i])
	p = libcl.common.p(rho_d[i], rv[i], Temp[i])
	p_v = rho_d[i] * rv[i] * libcl.common.R_v * Temp[i]
	S[i] = p_v / libcl.common.p_vs(Temp[i]) - 1

def rv2absS(del_S, rho_d, th_d, rv):
    for i in range(len(rho_d)):
        Temp = libcl.common.T(th_d[i], rho_d[i])
        pvs =  libcl.common.p_vs(Temp)
        rvs = pvs / (rho_d[i] * libcl.common.R_v * Temp)
        del_S[i] = rv[i] - rvs

def absS2rv(del_S, rho_d, th_d, rv):
    for i in range(len(rho_d)):
        Temp = libcl.common.T(th_d[i], rho_d[i])
        pvs =  libcl.common.p_vs(Temp)
        rvs = pvs / (rho_d[i] * libcl.common.R_v * Temp)
        rv[i] = del_S[i] + rvs


def main(scheme, 
  nx=300, sl_sg = slice(50,100), crnt=0.1, dt=0.2, nt=251, outfreq=50,
  aerosol={
    "meanr":.02e-6, "gstdv":1.4, "n_tot":550e6, 
    # ammonium sulphate aerosol parameters:
    "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
    "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
    "sd_conc":512. #TODO trzeba tu?
  }
):
    th_d = np.ones((nx,))* 303.
    #th_d[sl_sg] += 1
    rv = np.ones((nx,))* 5.e-3
    rc = np.zeros((nx,))
    rv[sl_sg] += 8.e-3
    #rc[sl_sg] += 1.e-3 
    rho_d = np.ones((nx,))*.97
    testowa = np.zeros((nx,))
    testowa[sl_sg]= 1.e4

    S     = np.empty((nx,))
    del_S = np.empty((nx,))
    Temp  = np.empty((nx,))
    
    var_adv = [th_d, del_S, testowa]
 
    if scheme == "2m":
           rr = np.zeros((nx,))
           nc = np.zeros((nx,))
           #nc[sl_sg] = 3.e7
           nr = np.zeros((nx,))
           var_adv = var_adv + [rc, rr, nc, nr]

    if scheme == "sd":
           na = np.zeros((nx,))
           nc = np.zeros((nx,))
           rc = np.zeros((nx,))
           sd = np.zeros((nx,))
           micro = libcl_spdr_init(rho_d, th_d, rv, crnt, dt, aerosol)
        

    calc_S(S, Temp, rho_d, th_d, rv)
    dic_var = {"nc":nc, "rc":rc, "rv":rv, "th":th_d, "Temp":Temp, "S":S}

    if scheme == "sd":
           dic_var["na"] = na
           dic_var["sd"] = sd

    plotting(dic_var, figname=scheme+"plot_init.pdf", time="init") 
    for it in range(nt):
        print "it", it
        rv2absS(del_S, rho_d, th_d, rv)
        if it==246: pdb.set_trace()
        print "testowa min, max przed adv", testowa.min(), testowa.max()
        if scheme == "2m": print "qc min, max przed adv", rc.min(), rc.max()        
        for var in var_adv:
            libmpdata.mpdata(var, crnt, 1);
        if scheme == "2m": print "qc min, max po adv", rc.min(), rc.max()
        print "testowa min, max po adv", testowa.min(), testowa.max()
        if it==246: pdb.set_trace()
        absS2rv(del_S, rho_d, th_d, rv)

        #TODO przerobic S na th_D
 
        if scheme == "2m":
            libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt, aerosol)
        elif scheme == "sd":
            libcl_spdr(rho_d, th_d, rv, rc, nc, na, sd, dt, aerosol, micro)
        else: 
            assert(False)

        calc_S(S, Temp, rho_d, th_d, rv) # tu juz sie przerabia na S 
                
        print "testowa po it = ", it
        if it % outfreq == 0 or it in [100]:
            plotting(dic_var, figname=scheme+"plot_"+str(int(it*dt))+"s.pdf", 
              time=str(int(it*dt))+"s" 
            )
            #plotting(dic_var, figname=scheme+"plot_"+str(int(it*dt))+"s_ylim.pdf",
            #         time=str(int(it*dt))+"s", ylim_dic={"S":[-0.005, 0.015], "nc":[4.86e7, 4.92e7], "rv":[0.0119,0.0121], "rc":[0.00098, 0.00104]} )



main("2m")  
#main("sd")
