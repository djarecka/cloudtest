import sys
sys.path.append(".")
sys.path.append("wrf_microphys/kessler/")
sys.path.append("libcloudphxx/")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import inspect
import pdb
import numpy as np
from numpy import array as arr_t
import analytic_blk_1m_pytest as an_blk
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# typical values as an example
press_0 = arr_t([900.e2  ])
th_0   = arr_t([291.8])
T_0    = arr_t([283.15])
rv_0   = arr_t([8.e-3])
rc_0   = arr_t([5.e-4])
rr_0   = arr_t([0.  ])
dt_0   = 1

def condensation(lib, press = None, T = None,
                 rv = None, rc = None, rr = None, dt=dt_0):

    import importlib
    lib_adj = importlib.import_module(lib)
    
    print "\n na poczatku srawdzam, czy ma qv", rv
    #pdb.set_trace()
    T = T if T!=None else T_0.copy()
    rv = rv if rv!=None else rv_0.copy()
    rc = rc if rc!=None else rc_0.copy()
    rr = rr if rr!=None else rr_0.copy()
    press = press if press!=None else press_0.copy()
    print "\n In condensation. Who is calling..?", inspect.stack()[1][3]
    rv, rc = lib_adj.adj_cellwise(press, T, rv, rc, rr, dt)
    return rv, rc




#TODO tylko dla przypadku gdy rc=0...
def cond_all(libname, libname_wrf, rv_0_min, rv_0_max,
                                  rc_0=0.):
    print "Arange", np.arange(rv_0_min, rv_0_max, 1.e-4)
    rc_an_l = []
    rc_xx_l = []
    rc_wrf_l = []
    r_rsat_an_l = []
    rsat_an = an_blk.mixrat_sat(T_0[0], press_0[0]) 
    for rv_0 in np.arange(rsat_an, rv_0_max, 2.e-4):
        print "\n w test_expected value przed",rv_0, rc_0
        rv_xx, rc_xx = condensation(lib=libname, rv=arr_t(rv_0), rc=arr_t(rc_0))
        rv_wrf, rc_wrf = condensation(lib=libname_wrf, rv=arr_t(rv_0), rc=arr_t(rc_0))
        delta_an = an_blk.delta_r(rv_0, T_0[0], press_0[0]) 
        rv_an, rc_an = rv_0 -  delta_an,  delta_an
        r_rsat_an = rv_0 - rsat_an
        #TODO bardzo dziwne sa te tablice...0-d array nie mozna indeksowac
        #import pdb
        #pdb.set_trace()
        rc_an_l.append(rc_an)
        rc_xx_l.append(rc_xx)
        rc_wrf_l.append(rc_wrf)
        r_rsat_an_l.append(r_rsat_an)
    print "rc po", rc_xx_l, rc_wrf_l, rc_an_l
    return r_rsat_an_l, rc_xx_l, rc_wrf_l, rc_an_l 
#    for key, value in expected.items():
#        print "\n key, valuu, eval(key)", key, value, eval(key)
#        assert abs(eval(key) - value) <= epsilon * abs(value)


def figure_plot():
    Libname = "libcloudphxx_blk_1m_pytest"
    Libname_wrf = "wrfkessler_blk_1m_pytest" #TODO: nie widzi libkessler.so (musialam przekopiowac)
    r_rsat, c_xx, c_wrf, c_an = cond_all(Libname, Libname_wrf, rv_0_min = 6.e-3, rv_0_max = 12.e-3)

    fig = plt.figure(1, figsize = (6.5,5))
    plt.plot(r_rsat, c_an, "b", r_rsat, c_xx, "r", r_rsat, c_wrf, "g" )
    plt.legend(["analytic","libcloudphxx", "wrf_kessler"],
               loc=2, prop = FontProperties(size=10))
    plt.xlabel(r'$rv_0 - rv_{sat0}$', fontsize=14)
    plt.ylabel(r'$\Delta rc$', fontsize=14)
    plt.savefig("cond_comp.pdf")
    plt.show()

figure_plot()
