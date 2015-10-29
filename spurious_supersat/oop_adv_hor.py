import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math
import json

from adv_hor import plotting, saving_state

import init_WH as wh
import init_WH_rhoconst as wh_rho
import init_slow_act as sl_act
import init_default as dft

class Micro:
    def __init__(self, nx, sl_sg, apr, setup, scheme):
        self.nx = nx
        self.sl_sg = sl_sg
        self.apr = apr
        self.scheme = scheme
        self.setup = setup
        if self.setup=="rhoconst":
            self.state, self.var_adv = dft.thermo_init(nx=self.nx, sl_sg=self.sl_sg,
                                                          scheme=self.scheme, apr=self.apr)
        elif self.setup=="slow_act":            
            self.state, self.var_adv = sl_act.thermo_init(nx=self.nx, sl_sg=self.sl_sg,
                                                          scheme=self.scheme, apr=self.apr)
        
    def calc_S(self):
        for i in range(len(self.state["S"])):
            self.state["Temp"][i] = libcl.common.T(self.state["th_d"][i], self.state["rho_d"][i])
            p = libcl.common.p(self.state["rho_d"][i], self.state["rv"][i], self.state["Temp"][i]) #TODO needed?
            p_v = self.state["rho_d"][i] * self.state["rv"][i] * libcl.common.R_v * self.state["Temp"][i] #TODO pomyslec o rho_D
            self.state["S"][i] = p_v / libcl.common.p_vs(self.state["Temp"][i]) - 1
                                            

    def rv2absS(self): #czy trzeba updatowac roho_d
        for i in range(len(self.state["rho_d"])):
            Temp = libcl.common.T(self.state["th_d"][i], self.state["rho_d"][i])
            pvs =  libcl.common.p_vs(Temp)
            rvs = pvs / (self.state["rho_d"][i] * libcl.common.R_v * Temp)
            self.state["del_S"][i] = self.state["rv"][i] - rvs
            
    def absS2rv(self):
        for i in range(len(self.state["rho_d"])):
            Temp = libcl.common.T(self.state["th_d"][i], self.state["rho_d"][i])
            pvs =  libcl.common.p_vs(Temp)
            rvs = pvs / (self.state["rho_d"][i] * libcl.common.R_v * Temp)
            self.state["rv"][i] = self.state["del_S"][i] + rvs
                                                                        
    def advection(self):
        for var in self.var_adv:
            libmpdata.mpdata(self.state[var], self.C, 1);
                                    

class Superdroplet(Micro):
    def __init__(self, nx, sl_sg, apr, C, dt, nt, sl_act_it, aerosol, setup="rhoconst", dx=2):

        Micro.__init__(self, nx, sl_sg, apr, setup, scheme="sd")
         
        #TODO czy to sie updatuje??/
        self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc"))
                
        self.aerosol = aerosol
        self.dt = dt
        self.nt = nt
        self.sl_act_it = sl_act_it
        self.C = C
        self.sl_sg = sl_sg
        
        self.opts_init = libcl.lgrngn.opts_init_t()
        self.opts_init.dt = dt
        self.opts_init.nx = self.state["rv"].shape[0] #TODO to jest chyba nx??
        self.opts_init.dx = dx
        self.opts_init.x1 = self.opts_init.nx * self.opts_init.dx

        def lognormal(lnr):
            from math import exp, log, sqrt, pi
            return self.aerosol["n_tot"] * exp(
                -(lnr - log(self.aerosol["meanr"]))**2 / 2 / log(self.aerosol["gstdv"])**2
            ) / log(self.aerosol["gstdv"]) / sqrt(2*pi)
                
        self.opts_init.sd_conc_mean = self.aerosol["sd_conc"]
        self.opts_init.dry_distros = {self.aerosol["kappa"]:lognormal}
        
        self.opts_init.coal_switch = self.opts_init.sedi_switch = False
        
        self.opts_init.sstp_cond = 5
        self.micro = libcl.lgrngn.factory(libcl.lgrngn.backend_t.serial, self.opts_init)
        
        self.C_arr = np.ones(self.state["rv"].shape[0]+1) * C
        self.micro.init(self.state["th_d"], self.state["rv"], self.state["rho_d"], self.C_arr)


    def micro_step(self, adve=True):
        libopts = libcl.lgrngn.opts_t()
        libopts.cond = True
        libopts.adve = adve
        libopts.coal = libopts.sedi = False
        
        self.micro.step_sync(libopts, self.state["th_d"], self.state["rv"], self.state["rho_d"])
        self.micro.step_async(libopts)
        
        # absolute number of super-droplets per grid cell
        self.micro.diag_sd_conc()
        self.state["sd"][:] = np.frombuffer(self.micro.outbuf())
        
        # number of particles (per kg of dry air) with r_w < .5 um
        self.micro.diag_wet_rng(0, .5e-6)
        self.micro.diag_wet_mom(0)
        self.state["na"][:] = np.frombuffer(self.micro.outbuf())
        
        # number of particles (per kg of dry air) with r_w > .5 um
        self.micro.diag_wet_rng(.5e-6, 1)
        self.micro.diag_wet_mom(0)
        self.state["nc"][:] = np.frombuffer(self.micro.outbuf())
        
        # cloud water mixing ratio [kg/kg] (same size threshold as above)
        self.micro.diag_wet_mom(3)
        rho_H2O = 1e3
        self.state["rc"][:] = 4./3 * math.pi * rho_H2O * np.frombuffer(self.micro.outbuf())

    def slow_act(self):
        self.state["rv"][self.sl_sg] += self.state["rv_sl_act"] / self.sl_act_it


    def all_sym(self):
        self.calc_S()
        print "po S", self.state["S"].max()
        plotting(self.dic_var, figname="Testplot_init.pdf", time="init")

        for it in range(self.nt+self.sl_act_it):
            print "it", it
            if it < self.sl_act_it:
                self.slow_act()
                self.micro_step(adve=False)
            else:
                if self.apr in ["S_adv", "S_adv_adj"]: self.rv2absS()
                ss.advection()
                if self.apr in ["S_adv", "S_adv_adj"]: self.absS2rv()
                print "po adv rv", self.state["rv"].max()
                self.micro_step()
                #if apr == "S_adv_adj": self.micro_adj() #TODO dolaczyc metode micro_adjust
            self.calc_S()
            
                
        plotting(self.dic_var, figname="Testplot_init.pdf", time=str(it), )#ylim_dic={"S":[-0.005, 0.015]})
        saving_state(self.dic_var, filename="test_dane.txt")#scheme+"_"+apr+"_"+setup+"_"+"data_"+str(int(it*dt))+"s.txt")
        
        
        
ss  = Superdroplet(nx=300, sl_sg=slice(50,100), apr="trad", C=0.1, dt=0.4, nt=200,
                    aerosol={
                        "meanr":.02e-6, "gstdv":1.4, "n_tot":1000e6,
                        "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
                        "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
                        "sd_conc":512 #TODO trzeba tu?
                        }, sl_act_it=200, setup="slow_act")

ss.all_sym()
