import sys
import oop_adv_hor as ah
import numpy as np
import pdb

rho_d = np.ones(5)
th_d = 295 * np.ones(5) 
th_d[1:3] += 5
rv =  np.arange(7, 12, 0.5) * 1.e-3
del_S = np.empty(5)

def test_supersat():
    rho_d0 = rho_d.copy()
    th_d0  = th_d.copy()
    rv0    = rv.copy()
    ah.rv2absS(del_S, rho_d, th_d, rv)
    ah.absS2rv(del_S, rho_d, th_d, rv)
#    pdb.set_trace()
    assert (rv0 == rv).all()
