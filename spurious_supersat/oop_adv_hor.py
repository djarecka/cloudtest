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

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


import init_WH as wh
import init_WH_rhoconst as wh_rho
import init_slow_act as sl_act
import init_default as dft

#TODO przeniesc?
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

def plotting_timeevol(dct_max, dct_mean, figname="evol_test.pdf", it0=0, it_step=1, show=False):
    nrow = (len(dct_max))
    fig, tpl = plt.subplots(nrows=nrow, ncols=1, figsize=(10,8.5))
    i=0
    for k,v in dct_max.iteritems():
        tpl_el = tpl[i%nrow]
        tpl_el.set_title(k)
        tpl_el.plot(v[it0::it_step], "b")
        if k in dct_mean.keys():
            tpl_el.plot(dct_mean[k][it0::it_step], "r")
        i+=1
    plt.savefig(figname)
    if show: plt.show()
    

class Micro:
    def __init__(self, nx, dx, time_adv_tot, dt, C, aerosol, RHenv, sl_sg, apr, setup, n_intrp, sl_act_time, dir_name, test, it_output_l, scheme):
        self.nx = nx
        self.dx = dx
        self.nt = int(time_adv_tot/dt)
        self.dt = dt
        self.C = C
        self.aerosol = aerosol
        self.RHenv = RHenv
        
        self.sl_sg = sl_sg
        self.apr = apr
        self.scheme = scheme
        self.setup = setup
        
        if self.setup=="rhoconst":
            self.sl_act_it = 0
            self.state, self.var_adv = dft.thermo_init(nx=self.nx, sl_sg=self.sl_sg,
                                                          scheme=self.scheme, apr=self.apr)
        elif self.setup=="slow_act":
            self.sl_act_it = int(sl_act_time/self.dt)
            self.state, self.var_adv = sl_act.thermo_init(RHenv=self.RHenv, nx=self.nx,
                                                          sl_sg=self.sl_sg,
                                                          scheme=self.scheme, apr=self.apr)
        self.n_intrp = n_intrp
        self.it_output = it_output_l
        self.test = test    
        if not self.test:
            if setup=="slow_act": dir_name += "_sl_act_time=" + str(sl_act_time) + "s"
            if self.n_intrp>1: dir_name += "_nintrp=" + str(self.n_intrp)
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
                                    
    def epsilon_adj(self): #TODO musze zrobic lepiej
        L = 2.5e6 #TODO
        for i in range(len(self.state["rho_d"])):
            Temp = libcl.common.T(self.state["th_d"][i], self.state["rho_d"][i])
            p = libcl.common.p(self.state["rho_d"][i], self.state["rv"][i], Temp)
            pvs =  libcl.common.p_vs(Temp)
            rvs = pvs / (self.state["rho_d"][i] * libcl.common.R_v * Temp)
            drs_dT = L * rvs / (libcl.common.R_v * Temp**2)
            Gamma = 1  + drs_dT * L / libcl.common.c_pd
            print "UWAGA: EPS ZMNIEJSZONY"
            eps = 1.e0 * 1.*(self.state["rv"][i] - rvs - self.state["del_S"][i]) / Gamma
            #eps = max(eps, -self.state["rc"][i])
            #if self.state["rv"][i] / rvs < 1:
            #    eps = min(0., eps)
            self.state["eps"][i] = eps
            if self.state["rc"][i]>0: pdb.set_trace()
            #TODO czy to tu??
            self.state["rv"][i] -= self.state["eps"][i]
            ##rc[i] += epsilon
            th = libcl.common.th_dry2std(self.state["th_d"][i], self.state["rv"][i])
            th += L/libcl.common.c_pd * (libcl.common.p_1000/p)**(libcl.common.R_d/libcl.common.c_pd) * self.state["eps"][i]
            self.state["th_d"][i] = libcl.common.th_std2dry(th, self.state["rv"][i])
        self.rc_adjust()
