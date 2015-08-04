import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math
import json 

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import init_WH as wh

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

def saving_state(dic, filename):
    dic_list = {}
    for key, value in dic.iteritems():
        dic_list[key] = value.tolist()
    f_w = open(filename, 'w')
    json.dump(dic_list, f_w)


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

    print "qc min, max przed mikro", rc.min(), rc.max()
    print "nc min, max przed mikro", nc.min(), nc.max()
    
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, thd, rv, rc, nc, rr, nr, dt)
    
    print "rc min po mikro" ,(rc + dot_rc * dt).min(), (rc + dot_rc * dt).max()
    print "nc min po mikro" ,(nc + dot_nc * dt).min(), (nc + dot_nc * dt).max()


    thd  += dot_th * dt
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
    opts.cond = opts.cevp = True
    opts.revp = opts.conv = opts.accr = opts.sedi = False
    print "1m rc max, min przed mikro", rc.max(), rc.min()
    libcl.blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    print "1m rc max, min po mikro", rc.max(), rc.min()


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

    opts_init.sd_conc = aerosol["sd_conc"]
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
	p = libcl.common.p(rho_d[i], rv[i], Temp[i]) #TODO needed?
	p_v = rho_d[i] * rv[i] * libcl.common.R_v * Temp[i] #TODO pomyslec o rho_D
	S[i] = p_v / libcl.common.p_vs(Temp[i]) - 1
        #pdb.set_trace()

def rv2absS(del_S, rho_d, th_d, rv): #czy trzeba updatowac roho_d
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

def micro_adj(rv, rc, th_d, rho_d, del_S): #TODO musze zrobic lepiej
    L = 2.5e6 #TODO
    #pdb.set_trace()
    for i in range(len(rho_d)):
        Temp = libcl.common.T(th_d[i], rho_d[i])
        p = libcl.common.p(rho_d[i], rv[i], Temp)
        pvs =  libcl.common.p_vs(Temp)
        rvs = pvs / (rho_d[i] * libcl.common.R_v * Temp)
        drs_dT = L * rvs / (rv[i] * Temp**2)  
        Gamma = 1  + drs_dT * L / libcl.common.c_pd
        epsilon = (rv[i] - rvs - del_S[i]) / Gamma
        #if i == 60:
        #    pdb.set_trace()
        rv[i] -= epsilon
        rc[i] += epsilon
        th = libcl.common.th_dry2std(th_d[i], rv[i])
        th += L/libcl.common.c_pd * (libcl.common.p_1000/p)**(libcl.common.R_d/libcl.common.c_pd) * epsilon
        th_d[i] = libcl.common.th_std2dry(th, rv[i])
    #pdb.set_trace()

def thermo_init(nx, sl_sg, scheme, apr):
    state = {}
    state["th_d"] = np.ones((nx,))* 303.
    #state["th_d"][sl_sg] += 1                                                              
    state["rv"] = np.ones((nx,))* 5.e-3
    state["rc"] = np.zeros((nx,))
    state["rv"][sl_sg] += 8.e-3
    #state["rc"][sl_sg] += 1.e-3                                           
    state["rho_d"] = np.ones((nx,))*.97
    state["testowa"] = np.zeros((nx,))
    state["testowa"][sl_sg]= 1.e4

    state["S"]    = np.empty((nx,))
    state["del_S"] = np.empty((nx,))
    state["Temp"] = np.empty((nx,))

    if apr == "trad":
        var_adv = ["th_d", "rv", "testowa"]
    elif apr in ["S_adv", "S_adv_adj"]:
        var_adv = ["th_d", "del_S", "testowa"]
    else:
        assert(False)

    if scheme in ["1m", "2m"]:
           state["rr"] = np.zeros((nx,))
           var_adv = var_adv + ["rc", "rr"]

    if scheme == "2m":
           state["nc"] = np.zeros((nx,))
           #state["nc"][sl_sg] = 3.e7                                            
           state["nr"] = np.zeros((nx,))
           var_adv  = var_adv + ["nc", "nr"]

    if scheme == "sd":
           state["na"] = np.zeros((nx,))
           state["nc"] = np.zeros((nx,))
           state["rc"] = np.zeros((nx,))
           state["sd"] = np.zeros((nx,))
 
    return state, var_adv



def main(scheme, apr="trad", setup="rhoconst", pl_flag = False, 
  nx=300, sl_sg = slice(50,100), crnt=0.1, dt=0.2, nt=501, outfreq=500,
  aerosol={
    "meanr":.02e-6, "gstdv":1.4, "n_tot":1000e6, 
    # ammonium sulphate aerosol parameters:
    "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
    "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
    "sd_conc":512 #TODO trzeba tu?
  }
):

    if setup == "rhoconst":
        state, var_adv = thermo_init(nx, sl_sg, scheme, apr)
    elif setup =="wh":
        state, var_adv = wh.thermo_init(nx, sl_sg, scheme, apr)
    else:
        assert(False)

    if scheme == "sd":
        micro = libcl_spdr_init(state["rho_d"], state["th_d"], state["rv"], crnt, dt, aerosol)
        
    #pdb.set_trace()
    calc_S(state["S"], state["Temp"], state["rho_d"], state["th_d"], state["rv"])

    if scheme == "1m":
        dic_var = dict((k, state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S"))
    elif scheme == "2m":
        dic_var = dict((k, state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc"))
    elif scheme == "sd":
        dic_var = dict((k, state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc", "na", "sd"))
    else:
        assert(False)
    
    if pl_flag: plotting(dic_var, figname=scheme+"_"+apr+"_"+setup+"_"+"plot_init.pdf", time="init") 
    saving_state(dic_var, filename=scheme+"_"+apr+"_"+setup+"_"+"data_init.txt")
    for it in range(nt):
        print "it", it
        if apr in ["S_adv", "S_adv_adj"]: rv2absS(state["del_S"], state["rho_d"], state["th_d"], state["rv"])
        print "testowa min, max przed adv", state["testowa"].min(), state["testowa"].max()
        if scheme == "2m": print "qc min, max przed adv", state["rc"].min(), state["rc"].max()        
        for var in var_adv:
            libmpdata.mpdata(state[var], crnt, 1);
        if scheme == "2m": print "qc min, max po adv", state["rc"].min(), state["rc"].max()
        print "testowa min, max po adv", state["testowa"].min(), state["testowa"].max()
        if setup == "wh": wh.rho_adjust(state, nx)


        if apr in ["S_adv", "S_adv_adj"]: absS2rv(state["del_S"], state["rho_d"], state["th_d"], state["rv"])
 

        if   scheme == "1m":
            libcl_1mom(state["rho_d"], state["th_d"], state["rv"], state["rc"], state["rr"], 
                       dt)
        elif scheme == "2m":
            libcl_2mom(state["rho_d"], state["th_d"], state["rv"], state["rc"], state["rr"], 
                       state["nc"], state["nr"], dt, aerosol)
            #pdb.set_trace()
        elif scheme == "sd":
            libcl_spdr(state["rho_d"], state["th_d"], state["rv"], state["rc"], 
                       state["nc"], state["na"], state["sd"], dt, aerosol, micro)
        else: 
            assert(False)


        if apr == "S_adv_adj": micro_adj(state["rv"], state["rc"], state["th_d"], state["rho_d"], state["del_S"])
        #pdb.set_trace()
        if setup == "wh": wh.rho_adjust(state, nx)
        calc_S(state["S"], state["Temp"], state["rho_d"], state["th_d"], state["rv"]) 
                
        print "testowa po it = ", it
        if it % outfreq == 0 or it in [100]:
            if pl_flag: plotting(dic_var, figname=scheme+"_"+apr+"_"+setup+"_"+"plot_"+str(int(it*dt))+"s.pdf", 
              time=str(int(it*dt))+"s" 
            )
            if pl_flag: plotting(dic_var, figname=scheme+"_"+apr+"_"+setup+"_"+"plot_"+str(int(it*dt))+"s_ylim.pdf",
                     time=str(int(it*dt))+"s", ylim_dic={"S":[-0.005, 0.015], "nc":[5.e8, 6.e8], "rv":[0.0108,0.0112], "rc":[0.00095, 0.0011]} )
            if it == nt-1:
                saving_state(dic_var, filename=scheme+"_"+apr+"_"+setup+"_"+"data_"+str(int(it*dt))+"s.txt")


if __name__ == '__main__':
    #main("2m", pl_flag=True) 
    main("2m", pl_flag=True, setup="wh")
    #main("1m")
    #main("sd",pl_flag=True)
    #main("2m", apr="S_adv", setup="wh", pl_flag=True)
    #main("2m", apr="S_adv", pl_flag=True)
    #main("sd", apr="S_adv", pl_flag=True)
    #main("2m", apr="S_adv_adj", pl_flag=True)
    #main("sd", apr="S_adv_adj", pl_flag=True)
