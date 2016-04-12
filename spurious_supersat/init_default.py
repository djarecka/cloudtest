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
    state["th_d"] = np.ones((nx,))* 303.
    #state["th_d"][sl_sg] += 1                                                              
    state["rv"] = np.ones((nx,))* 5.e-3
    state["rc"] = np.zeros((nx,))
    state["rv"][sl_sg] += 7.e-3
    #state["rc"][sl_sg] += 1.e-3                                           
    state["rho_d"] = np.ones((nx,))*.97
    state["testowa"] = np.zeros((nx,))
    state["testowa"][sl_sg]= 1.e4

    state["S"]    = np.empty((nx,))
    state["del_S"] = np.empty((nx,))
    state["Temp"] = np.empty((nx,))

    if apr in ["S_adv_adj"]: state["eps"] = np.zeros((nx,))
    
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



