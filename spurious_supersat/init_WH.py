import sys
sys.path.append(".")
sys.path.append("../")

import pdb

import libmpdata
import libcloudphxx as libcl
import numpy as np
import math

def thermo_init(nx, sl_sg, scheme, apr, press=8.e4):
    state = {}
    state["th_d"] = np.empty((nx,))
    state["rho_d"] = np.empty((nx,))
    state["rv"] = np.empty((nx,))
    state["rc"] = np.empty((nx,))

    state["testowa"] = np.zeros((nx,))
    state["testowa"][sl_sg]= 1.e4

    state["S"]    = np.empty((nx,))
    state["del_S"] = np.empty((nx,))
    state["Temp"]  = np.empty((nx,))
    state["nc"] = np.zeros((nx,))

    # to wyliczenia rho zakladam, ze rv jest zero, nie wiem chwilowo jak inaczej
    for ii in range(nx):
        if ii in range(sl_sg.start, sl_sg.stop):
            th_std = 303.8
            RHr    = 1.0015
            rc = 1.e-3
            nc = 550.e6
        else:
            th_std =  302.8
            RHr    = 0.5
            rc = 0.
            nc = 0

        state["Temp"][ii] = th_std / (libcl.common.p_1000/press)**(libcl.common.R_d/libcl.common.c_pd)
        #pdb.set_trace()
        p_vs = libcl.common.p_vs(state["Temp"][ii])
        rvs = libcl.common.eps * p_vs / (press - p_vs)
        state["rv"][ii] = RHr * rvs
        state["th_d"][ii] = libcl.common.th_std2dry(th_std, state["rv"][ii])
        state["rho_d"][ii] = libcl.common.rhod(press,
           libcl.common.th_dry2std(state["th_d"][ii], state["rv"][ii]),
                                               state["rv"][ii]) #TODO zwykle th?

        state["rc"][ii] = rc
        state["nc"][ii] = nc



    #pdb.set_trace()
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
        state["nr"] = np.zeros((nx,))
        var_adv  = var_adv + ["nc", "nr"]

    if scheme == "sd":
           state["na"] = np.zeros((nx,))
           state["nc"] = np.zeros((nx,))
           state["rc"] = np.zeros((nx,))
           state["sd"] = np.zeros((nx,))
 
    return state, var_adv


def rho_adjust(state, nx, press=8.e4):
    for ii in range(nx):
        state["rho_d"][ii] = libcl.common.rhod(8.e4,
                       libcl.common.th_dry2std(state["th_d"][ii], state["rv"][ii]),
                                               state["rv"][ii])



