import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math
import json
import os
from subprocess import call

from adv_hor import plotting, saving_state

import init_WH as wh
import init_WH_rhoconst as wh_rho
import init_slow_act as sl_act
import init_default as dft

class Micro:
    def __init__(self, nx, dx, time_adv_tot, dt, sl_sg, apr, setup, scheme, n_intrp, sl_act_time, dir_name, test):
        self.nx = nx
        self.dx = dx
        self.nt = int(time_adv_tot/dt)
        self.dt = dt
                                
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
        self.test = test    
        if not self.test:
            if setup=="slow_act": dir_name += "_sl_act_time=" + str(sl_act_time) + "s"
            if n_intrp>1: dir_name += "_nintrp=" + str(n_intrp)
            self.outputdir = os.path.join("output", dir_name)
            if os.path.exists(self.outputdir):
                call(["rm", "-r", self.outputdir])
            call(["mkdir", self.outputdir])

        # I only want to remove the dir if self.test=False
        self.plotdir = os.path.join("plots", dir_name)
        if os.path.exists(self.plotdir) and not self.test:
            call(["rm", "-r", self.plotdir])
            call(["mkdir", self.plotdir])
        elif not os.path.exists(self.plotdir):
            call(["mkdir", self.plotdir])
                                
        
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
    def __init__(self, nx, dx, sl_sg, apr, C, dt, time_adv_tot, sl_act_time, aerosol, scheme="sd", setup="rhoconst", n_intrp=1, test=True):

        dir_name = "Test_scheme=sd_conc="+str(int(aerosol["sd_conc"]))+"_setup="+setup+"_apr="+apr+"_dt="+str(dt)+"_dx="+str(dx)+"_C="+str(C)
        Micro.__init__(self, nx, dx, time_adv_tot, dt, sl_sg, apr, setup, scheme, n_intrp, sl_act_time, dir_name, test)
         
        #TODO czy to sie updatuje??/ nie jak jest n_intrp>1 - pomyslec!
        self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc"))
                
        self.aerosol = aerosol
        self.n_intrp = n_intrp
        if setup=="rhoconst": self.sl_act_it = 0
        if setup=="slow_act": self.sl_act_it = int(sl_act_time/self.dt)
        self.C = C
        
        self.opts_init = libcl.lgrngn.opts_init_t()
        self.opts_init.dt = dt
        self.opts_init.nx = self.nx * self.n_intrp 
        self.opts_init.dx = dx / float(self.n_intrp)
        self.opts_init.x1 = self.opts_init.nx * self.opts_init.dx

        def lognormal(lnr):
            from math import exp, log, sqrt, pi
            return self.aerosol["n_tot"] * exp(
                -(lnr - log(self.aerosol["meanr"]))**2 / 2 / log(self.aerosol["gstdv"])**2
            ) / log(self.aerosol["gstdv"]) / sqrt(2*pi)
                
        self.opts_init.sd_conc = self.aerosol["sd_conc"] / self.n_intrp
        self.opts_init.dry_distros = {self.aerosol["kappa"]:lognormal}
        
        self.opts_init.coal_switch = self.opts_init.sedi_switch = False
        
        self.opts_init.sstp_cond = 1
        self.micro = libcl.lgrngn.factory(libcl.lgrngn.backend_t.serial, self.opts_init)

        self.state_micro = {}
        if self.n_intrp == 1:
            for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
                self.state_micro[var] = self.state[var]
        elif self.n_intrp > 1:
            self.interp_adv2micro()
            
        self.C_arr = np.ones(self.state_micro["rv"].shape[0]+1) * self.C * float(self.n_intrp)
        self.micro.init(self.state_micro["th_d"], self.state_micro["rv"], self.state_micro["rho_d"], self.C_arr)
        
    def interp_adv2micro(self):
        x_adv = np.arange(self.nx) + 1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1 / 2.) / float(self.n_intrp)
        #pdb.set_trace()
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state_micro[var] = np.interp(x_micro, x_adv, self.state[var])
            
    def interp_micro2adv(self):
        x_adv = np.arange(self.nx) +1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1. /2) / float(self.n_intrp)
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state[var] = np.interp(x_adv, x_micro, self.state_micro[var])
                                                    

    def fake_interp_adv2micro(self):
        x_adv = np.arange(self.nx) + 1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1 / 2.) / float(self.n_intrp)
        #pdb.set_trace()
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state_micro[var] = np.repeat(self.state[var], self.n_intrp)
            
    def fake_interp_micro2adv(self):
        x_adv = np.arange(self.nx) +1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1. /2) / float(self.n_intrp)
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            for i in range(self.state["rv"].shape[0]):
                self.state[var][i] = self.state_micro[var][i*self.n_intrp:(i+1)*self.n_intrp].mean()
        

            
    def micro_step(self, adve=True):
        libopts = libcl.lgrngn.opts_t()
        libopts.cond = True
        libopts.adve = adve
        libopts.coal = libopts.sedi = False

        self.micro.step_sync(libopts, self.state_micro["th_d"], self.state_micro["rv"], self.state_micro["rho_d"])
        self.micro.step_async(libopts)
        
        # absolute number of super-droplets per grid cell
        self.micro.diag_sd_conc()
        self.state_micro["sd"][:] = np.frombuffer(self.micro.outbuf())
        
        # number of particles (per kg of dry air) with r_w < .5 um
        self.micro.diag_wet_rng(0, .5e-6)
        self.micro.diag_wet_mom(0)
        self.state_micro["na"][:] = np.frombuffer(self.micro.outbuf())
        
        # number of particles (per kg of dry air) with r_w > .5 um
        self.micro.diag_wet_rng(.5e-6, 1)
        self.micro.diag_wet_mom(0)
        self.state_micro["nc"][:] = np.frombuffer(self.micro.outbuf())
        
        # cloud water mixing ratio [kg/kg] (same size threshold as above)
        self.micro.diag_wet_mom(3)
        rho_H2O = 1e3
        self.state_micro["rc"][:] = 4./3 * math.pi * rho_H2O * np.frombuffer(self.micro.outbuf())

    def slow_act(self):
        self.state["rv"][self.sl_sg] += self.state["rv_sl_act"] / self.sl_act_it


    def all_sym(self):
        self.calc_S()
        print "po S", self.state["S"].max()
        plotting(self.dic_var, figname=os.path.join(self.plotdir, "TestplotRH95_init.pdf"), time="init")
        if not self.test: saving_state(self.dic_var, filename=os.path.join(self.outputdir, "init.txt"))

        it_output = [50, self.sl_act_it/4., self.sl_act_it/2]
        for Dx in [0, 5, 10, 20, 50, 100, 150, 200]:
            it_out = Dx * 1./ self.C
            for i in range(1,int(1/self.C)+2):
                it_output.append(self.sl_act_it + it_out + i)
                                                                                    
        
        for it in range(self.nt+self.sl_act_it):
            print "it", it, it <= self.sl_act_it
            if it <= self.sl_act_it:
                
                self.slow_act()
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step(adve=False)
                if self.n_intrp > 1: self.interp_micro2adv()
            else:
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step()
                if self.n_intrp > 1: self.interp_micro2adv()
                if self.apr in ["S_adv", "S_adv_adj"]: self.rv2absS()
                ss.advection()
                if self.apr in ["S_adv", "S_adv_adj"]: self.absS2rv()
                            
                #pdb.set_trace()
                #if apr == "S_adv_adj": self.micro_adj() #TODO dolaczyc metode micro_adjust
            self.calc_S()
                               
                
            if it in it_output:
                #pdb.set_trace()
                self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'th_d',"na", "nc", "S"))
                plotting(self.dic_var, figname=os.path.join(self.plotdir, "newplot_slowit"+str(self.sl_act_it)+"_Crr"+str(self.C)+"_nintrp"+str(self.n_intrp)+"_it="+str(int(it*self.dt))+"s.pdf"), time=str(int(self.dt*(it-self.sl_act_it))), ylim_dic={"S":[-0.005, 0.015]})
                if not self.test: saving_state(self.dic_var, filename=os.path.join(self.outputdir, "it="+str(int(self.dt*(it-self.sl_act_it)))+"s.txt"))#scheme+"_"+apr+"_"+setup+"_"+"data_"+str(int(it*dt))+"s.txt")
                #pdb.set_trace()
        
        
ss  = Superdroplet(nx=320, dx=2, sl_sg=slice(50,100), apr="trad", C=.1, dt=.1, time_adv_tot=1001,
                    aerosol={
                        "meanr":.02e-6, "gstdv":1.4, "n_tot":1e9,
                        "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
                        "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
                        "sd_conc":256 #TODO trzeba tu?
                        }, sl_act_time=60, n_intrp=1, setup="slow_act", scheme="sd",
                   test=False)

ss.all_sym()
