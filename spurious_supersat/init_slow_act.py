import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math

def thermo_init(nx, sl_sg, scheme, apr):
    state = {}
    state["th_d"] = np.empty((nx,))
    state["rv"] = np.empty((nx,))
    state["rc"] = np.empty((nx,))

    state["testowa"] = np.zeros((nx,))
    state["testowa"][sl_sg]= 1.e4

    state["S"]    = np.empty((nx,))
    state["del_S"] = np.empty((nx,))
    state["Temp"]  = np.empty((nx,))
    state["nc"] = np.zeros((nx,))

    #daje podobne ustawienia jak WH, ale zakladam stale rho_d
    # przyjmuje, ze poczatkowo podane sa th jest bliska thd
    # korzystam z cisnienia i wartosci temperatury w "suchej czesci"
    th_0 =  302.8
    press_0 = 8.e4
    Temp_0 = th_0 / (libcl.common.p_1000/press_0)**(libcl.common.R_d/libcl.common.c_pd)
    rho_d_0 = press_0 / (Temp_0 * libcl.common.R_d)
    state["rho_d"] = rho_d_0 * np.ones((nx,))


    for ii in range(nx):
        if ii in range(sl_sg.start, sl_sg.stop):
            #th_0 = 303.8
            RH_0    = 1.
            rc_0 = 0
            nc_0 = 0
        else:
            th_0 =  302.8
            RH_0    = 0.5
            rc_0 = 0.
            nc_0 = 0

#        if th_0==303.8:
#            pdb.set_trace()
        state["Temp"][ii] = libcl.common.T(th_0, state["rho_d"][ii])
        #pdb.set_trace()
        pvs = libcl.common.p_vs(state["Temp"][ii])
        rvs = pvs / (state["rho_d"][ii] * libcl.common.R_v * state["Temp"][ii])
        state["rv"][ii] = RH_0 * rvs #+ 2.e-3
        state["th_d"][ii] = th_0
        #if th_0==303.8:                                                               
        #    pdb.set_trace()                                                           


        state["rc"][ii] = rc_0
        state["nc"][ii] = nc_0

    state["rv_sl_act"] = .5e-3

    #pdb.set_trace()
    if apr == "trad":
        var_adv = ["th_d", "rv", "testowa"]
    elif apr in ["S_adv", "S_adv_adj"]:
        var_adv = ["th_d", "del_S", "testowa"]
    else:
        assert(False)

    if scheme in ["1m", "2m", "sd"]:
           state["rr"] = np.zeros((nx,))
           var_adv = var_adv + ["rc", "rr"]

    if scheme == "2m":
        state["nr"] = np.zeros((nx,))
        var_adv  = var_adv + ["nc", "nr"]

    if scheme == "sd":
           state["na"] = np.zeros((nx,))
           state["nc"] = np.zeros((nx,))
           state["rc"] = np.zeros((nx,))
           state["sd"] = np.zeros((nx,))
 
    return state, var_adv




