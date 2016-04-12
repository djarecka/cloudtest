import libmpdata
import libcloudphxx as libcl
import numpy as np
import math
import json
import os
from subprocess import call


import init_WH_rhoconst as wh_rho
import init_slow_act as sl_act
import init_default as dft

from oop_adv_hor import Micro, plotting, saving_state

import pdb

class Eul_2mom(Micro):
    def __init__(self, nx, dx, sl_sg, apr, C, dt, time_adv_tot, sl_act_time, aerosol, RHenv, setup="rhoconst", n_intrp=1, test=True, dirname_pre = "test", dirname=None, it_output_l=[]):
        if dirname:
            dir_name = dirname
        else:
            dir_name = dirname_pre+"_scheme=2mom_setup="+setup+"_apr="+apr+"_dt="+str(dt)+"_dx="+str(dx)+"_C="+str(C) + "_ntot="+str(int(aerosol["n_tot"]/1.e6))
        #pdb.set_trace()
        Micro.__init__(self, nx, dx, time_adv_tot, dt, C, aerosol, RHenv, sl_sg, apr, setup, n_intrp, sl_act_time, dir_name, test, it_output_l, scheme="2m")

        self.opts = libcl.blk_2m.opts_t()
        self.opts.acti = self.opts.cond = True
        self.opts.acnv = self.opts.accr = self.opts.sedi = False
        self.opts.dry_distros = [{
            "mean_rd" : aerosol["meanr"],
            "sdev_rd" : aerosol["gstdv"],
            "N_stp"   : aerosol["n_tot"],
            "chem_b"  : aerosol["chem_b"]
        }]
        
        self.state_dot = dict((k, np.zeros(self.nx)) for k in ('dot_rc', 'dot_rv', 'dot_th_d', "dot_nc", "dot_rr", "dot_nr"))
        self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', 'th_d', "Temp", "S", "nc"))
                
        
    def micro_step(self):
        print "qc min, max przed mikro", self.state["rc"].min(), self.state["rc"].max()
        print "nc min, max przed mikro", self.state["nc"].min(), self.state["nc"].max()
        for k in ("dot_th_d", "dot_rv", "dot_rc", "dot_nc", "dot_rr", "dot_nr"):
            self.state_dot[k] *= 0.
        #pdb.set_trace()
        libcl.blk_2m.rhs_cellwise(self.opts,
                                  self.state_dot["dot_th_d"], self.state_dot["dot_rv"],
                                  self.state_dot["dot_rc"], self.state_dot["dot_nc"],
                                  self.state_dot["dot_rr"], self.state_dot["dot_nr"],
                                  self.state["rho_d"], self.state["th_d"],
                                  self.state["rv"], self.state["rc"], self.state["nc"],
                                  self.state["rr"], self.state["nr"], self.dt)
        print "qc min, max po mikro", self.state["rc"].min(), self.state["rc"].max()
        print "nc min, max po mikro", self.state["nc"].min(), self.state["nc"].max()
                        

        for k in ("th_d", "rv", "rc", "nc", "rr", "nr"):
            self.state[k] += self.state_dot["dot_"+k] * self.dt

        #pdb.set_trace()
        np.place(self.state["rc"], self.state["rc"]<0, 0)
        np.place(self.state["nc"], self.state["nc"]<0, 0)
        print "qc min, max po place", self.state["rc"].min(), self.state["rc"].max()
        print "nc min, max po place", self.state["nc"].min(), self.state["nc"].max()
                        
    def rc_adjust(self):
        self.state["rc"] += self.state["eps"]
        
    def all_sym(self):
        self.calc_S()
        print "po S", self.state["S"].max()
        plotting(self.dic_var, figname=os.path.join(self.plotdir, "TestplotRH5_init.pdf"), time="init")
        if not self.test: saving_state(self.dic_var, filename=os.path.join(self.outputdir, "init.txt"))
        it_output = self.it_output

        max_state = {}
        meancl_state = {}
        for vv in ["S", "nc", "rc"]:
            max_state[vv] = []
            meancl_state[vv] = []
            
        for it in range(self.nt+self.sl_act_it):
            print "it", it, it < self.sl_act_it
                
            if it < self.sl_act_it:
                self.slow_act()
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step(adve=False)
                if self.n_intrp > 1: self.interp_micro2adv()
                if it==self.sl_act_it: all_water_i = self.state["rv"].sum() + self.state["rc"].sum()
            else:
                if self.apr in ["S_adv", "S_adv_adj"]: self.rv2absS()
                self.advection()
                if self.apr in ["S_adv", "S_adv_adj"]: self.absS2rv()
                if self.n_intrp > 1: self.interp_adv2micro()
                self.micro_step()
                if self.n_intrp > 1: self.interp_micro2adv()

                print "przed adjust", it, "%.25f" % self.state["rc"].sum(), "%.25f" % self.state["rv"].sum()
                if self.apr in ["S_adv_adj"]:
                    self.epsilon_adj()
                print "po adjust", it, "%.25f" % self.state["rc"].sum(), "%.25f" % self.state["rv"].sum()
                                
                if it==self.sl_act_it+self.nt-1: all_water_f = self.state["rv"].sum() + self.state["rc"].sum()
                    
            self.calc_S()

            for vv in ["S", "nc", "rc"]:
                max_state[vv].append(self.state[vv].max())
                meancl_state[vv].append(self.state[vv][np.where(self.state["rc"]>0)].mean())

            if it in it_output:
                #pdb.set_trace()
                self.dic_var = dict((k, self.state[k]) for k in ('rc', 'rv', "nc", "th_d", "S"))
                plotting(self.dic_var, figname=os.path.join(self.plotdir, "newplot_slowit"+str(self.sl_act_it)+"_Crr"+str(self.C)+"_nintrp"+str(self.n_intrp)+"_it="+str((it+1)*self.dt)+"s.pdf"), time=str(self.dt*(it-self.sl_act_it)), ylim_dic={"S":[-0.005, 0.015]})
                if not self.test: saving_state(self.dic_var, filename=os.path.join(self.outputdir, "it="+str(int(self.dt*(it+1-self.sl_act_it)))+"s.txt"))



if __name__ == '__main__':
    micro_2mom = Eul_2mom(nx=300, dx=2, sl_sg=slice(50,100), apr="S_adv_adj",
                          C=.2, dt=.1, time_adv_tot=21,
                          aerosol={
                              "meanr":.02e-6, "gstdv":1.4, "n_tot":1e9,
                              "chem_b":.505, 
                              "kappa":.61,    
                          },
                          RHenv=.95, sl_act_time=60, n_intrp=1, setup="rhoconst",
                          test=False, it_output_l=[50,100,200])
    
    micro_2mom.all_sym()
                                                                                     
