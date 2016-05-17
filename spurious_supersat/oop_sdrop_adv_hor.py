import sys
sys.path.append(".")
sys.path.append("../")

import libcloudphxx as libcl
import numpy as np
import math
import json
import os
import pdb

from oop_adv_hor import Micro, plotting_timeevol, plotting, saving_state 

class Superdroplet(Micro):
    """ Class for Lagrangian microphysics; inherits from Micro. """
    def __init__(self, nx, dx, sl_sg, apr, C, dt, time_adv_tot, sl_act_time, aerosol, RHenv, setup="rhoconst", n_intrp=1, sstp_cond=1, test=True, dirname_pre = "test", dirname=None, it_output_l = []):
        # creating a name for output dir 
        if dirname:
            dir_name = dirname
        else:
            dir_name = dirname_pre+"_scheme=sd_conc="+str(int(aerosol["sd_conc"]))+"_setup="+setup+"_apr="+apr+"_dt="+str(dt)+"_dx="+str(dx)+"_C="+str(C) + "_ntot="+str(int(aerosol["n_tot"]/1.e6))
        # calling Micro parent class
        Micro.__init__(self, nx, dx, time_adv_tot, dt, C, aerosol, RHenv, sl_sg, apr, setup, n_intrp, sl_act_time, dir_name, test, it_output_l, scheme="sd")
         
        #TODO czy to sie updatuje??/ nie jak jest n_intrp>1 - pomyslec!
        # creating a dictionary for model variables 
        self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc"))
        # calling a microphysics scheme from libcloudphxx library; defining opts  
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
        self.opts_init.n_sd_max = self.opts_init.sd_conc * self.opts_init.nx
        self.opts_init.dry_distros = {self.aerosol["kappa"]:lognormal}
        
        self.opts_init.coal_switch = self.opts_init.sedi_switch = False
        
        self.opts_init.sstp_cond = sstp_cond
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
        """ function for fields interpolation; smaller dx for microphysics """
        x_adv = np.arange(self.nx) + 1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1 / 2.) / float(self.n_intrp)
        #pdb.set_trace()
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state_micro[var] = np.interp(x_micro, x_adv, self.state[var])
            
    def interp_micro2adv(self):
        """ function for fields interpolation; biger dx for advection """
        x_adv = np.arange(self.nx) +1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1. /2) / float(self.n_intrp)
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state[var] = np.interp(x_adv, x_micro, self.state_micro[var])
                                                    

    def fake_interp_adv2micro(self):
        """ fake interpolation, only for testing """
        x_adv = np.arange(self.nx) + 1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1 / 2.) / float(self.n_intrp)
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            self.state_micro[var] = np.repeat(self.state[var], self.n_intrp)
            
    def fake_interp_micro2adv(self):
        """ fake interpolation, only for testing """
        x_adv = np.arange(self.nx) +1 / 2.
        x_micro = (np.arange(self.opts_init.nx) + 1. /2) / float(self.n_intrp)
        for var in ["rv", "th_d", "rho_d", "sd", "na", "nc", "rc"]:
            for i in range(self.state["rv"].shape[0]):
                self.state[var][i] = self.state_micro[var][i*self.n_intrp:(i+1)*self.n_intrp].mean()
            
    def micro_step(self, adve=True, cond=True):
        """ Defining microphysical step """
        libopts = libcl.lgrngn.opts_t()
        libopts.cond = cond
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

    def rc_adjust(self):
        self.micro.step_rc_adjust(self.state["eps"]) #TODO 
        # cloud water mixing ratio [kg/kg] (same size threshold as above)
        self.micro.diag_wet_mom(3)
        rho_H2O = 1e3
        self.state_micro["rc"][:] = 4./3 * math.pi * rho_H2O * np.frombuffer(self.micro.outbuf())
        
    def slow_act(self):
        if self.setup=="slow_act":
            self.state["rv"][self.sl_sg] += self.state["rv_sl_act"] / self.sl_act_it
        else:
            pass

    def all_sym(self):
        """ Full symulation """
        self.calc_S()
        print "po S", self.state["S"].max()
        plotting(self.dic_var, figname=os.path.join(self.plotdir, "TestplotRH5_init.pdf"), time="init")
        if not self.test: saving_state(self.dic_var, filename=os.path.join(self.outputdir, "init.txt"))

        it_output = self.it_output

        # calculating some statistics
        max_state = {}
        meancl_state = {}
        for vv in ["S", "nc", "rc"]:
            max_state[vv] = []
            meancl_state[vv] = []

        # time loop
        for it in range(self.nt+self.sl_act_it):
            print "it", it, it < self.sl_act_it
            # if slow_act setup, first part is without advection
            if it < self.sl_act_it:
                self.slow_act()
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step(adve=False)
                if self.n_intrp > 1: self.interp_micro2adv()
                if it==self.sl_act_it-1: all_water_i = self.state["rv"].sum() + self.state["rc"].sum()
            else:
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step()
                if self.n_intrp > 1: self.interp_micro2adv()
                if it==self.sl_act_it+self.nt-1: all_water_f = self.state["rv"].sum() + self.state["rc"].sum()
                print "przed adjust", it, "%.25f" % self.state["rc"].sum(), "%.25f" % self.state["rv"].sum()
                if self.apr in ["S_adv_adj"]:
                    self.epsilon_adj()

                print "po adjust", it, "%.25f" % self.state["rc"].sum(), "%.25f" % self.state["rv"].sum()
            #if apr == "S_adv_adj": self.micro_adj() #TODO dolaczyc metode micro_adjust

                if self.apr in ["S_adv", "S_adv_adj"]: self.rv2absS()
                self.advection()
                if self.apr in ["S_adv", "S_adv_adj"]: self.absS2rv()
                                                            

            self.calc_S()
            
            for vv in ["S", "nc", "rc"]:
                max_state[vv].append(self.state[vv].max())
                meancl_state[vv].append(self.state[vv][np.where(self.state["rc"]>0)].mean())
            # plotting and saving variables
            if it in it_output:
                #pdb.set_trace()
                self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'sd', "na", "nc", "th_d", "S"))
                plotting(self.dic_var, figname=os.path.join(self.plotdir, "newplot_slowit"+str(self.sl_act_it)+"_Crr"+str(self.C)+"_nintrp"+str(self.n_intrp)+"_it="+str((it+1)*self.dt)+"s.pdf"), time=str(self.dt*(it-self.sl_act_it)), ylim_dic={"S":[-0.005, 0.015]})
                if not self.test:
                    saving_state(self.dic_var, filename=os.path.join(self.outputdir, "it="+str(int(self.dt*(it+1-self.sl_act_it)))+"s.txt"))
    
        print "Nc_max", max_state["nc"]

        #print "all_water, init, final, rel_diff", all_water_i, all_water_f, (all_water_i-all_water_f)/all_water_i
        plotting_timeevol(max_state, meancl_state, figname=os.path.join(self.plotdir, "ewolucja_max.pdf"))
        plotting_timeevol(max_state, meancl_state, figname=os.path.join(self.plotdir, "ewolucja_max_fullstep.pdf"), it0=1, it_step=int(1/self.C))
        plotting_timeevol(max_state, meancl_state, figname=os.path.join(self.plotdir, "ewolucja_max_halfstep.pdf"), it0=1+int(0.5/self.C), it_step=int(1/self.C), show=False)
        
    

if __name__ == '__main__':
    ss  = Superdroplet(nx=300, dx=2, sl_sg=slice(50,100), apr="trad", C=.2, dt=.1, time_adv_tot=16, #76,
                       aerosol={
                           "meanr":.02e-6, "gstdv":1.4, "n_tot":1e9,
                           "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
                           "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
                           "sd_conc":256 #TODO trzeba tu?
                       }, RHenv=.95, sl_act_time=60, n_intrp=1, setup="slow_act",
                       test=False, it_output_l=[100, 649, 699, 849])
                           #600, 850, 851, 1100,1101, 1350, 1351,1700,1701])

    ss.all_sym()
