import sys
sys.path.append(".")
sys.path.append("wrf_microphys/kessler/")
sys.path.append("libcloudphxx/")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import inspect
import pdb
import numpy as np
import analytic_blk_1m_pytest as an
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([5.e-4])
rr_0   = np.array([0.  ])
nc_0   = np.array([1.e5]) #TODO
nr_0   = np.array([1.e5])#TODO
dt_0   = 1

def condensation(lib, lib_type, press = None, T = None,
                 rv = None, rc = None, rr = None,
                 nc = None, nr = None, dt=dt_0):

    import importlib
    lib_adj = importlib.import_module(lib)
    
    print "\n na poczatku srawdzam, czy ma qv", rv
    #pdb.set_trace()
    T = T if T!=None else T_0.copy()
    rv = rv if rv!=None else rv_0.copy()
    rc = rc if rc!=None else rc_0.copy()
    rr = rr if rr!=None else rr_0.copy()
    nc = nc if nc!=None else nc_0.copy()
    nr = nr if nr!=None else nr_0.copy()
    press = press if press!=None else press_0.copy()
    print "\n In condensation. Who is calling..?", inspect.stack()[1][3]
    if lib_type == 1:
        rv, rc, rr = lib_adj.adj_cellwise(press, T, rv, rc, rr, dt)
    elif lib_type == 2:
        rv, rc, nc, rr, nr = lib_adj.adj_cellwise(press, T, rv, rc, nc, rr, nr, dt)
    return rv, rc


#TODO tylko dla przypadku gdy rc=0...
def cond_all(libname, libname_2m, libname_wrf, rc_0, sup_lim, sup_step, 
             temp_0, press_0):
    rc_an_l = []
    rc_xx_l = []
    rc_2xx_l = []
    rc_wrf_l = []
    r_rsat_an_l = []
    rsat_an = an.mixrat_sat(temp_0, press_0) 
    for rv_0 in np.arange(rsat_an - sup_lim, rsat_an + sup_lim + sup_step, sup_step):
        print "\n w test_expected value przed",rv_0, rc_0
        rv_xx, rc_xx = condensation(lib=libname, lib_type=1, T=np.array(temp_0), 
                                    press = np.array(press_0),
                                    rv=np.array(rv_0), rc=np.array(rc_0))
        rv_2xx, rc_2xx = condensation(lib=libname_2m, lib_type=2, T=np.array(temp_0),
                                    press = np.array(press_0),
                                    rv=np.array(rv_0), rc=np.array(rc_0))
        rv_wrf, rc_wrf = condensation(lib=libname_wrf, lib_type=1, T=np.array(temp_0),
                                      press = np.array(press_0),
                                      rv=np.array(rv_0), rc=np.array(rc_0))
        #pdb.set_trace()
        delta_an = an.delta_r(rv_0, temp_0, press_0) 
        rv_an, rc_an = rv_0 - delta_an, rc_0 + delta_an
        r_rsat_an = rv_0 - rsat_an
        #pdb.set_trace()
        #TODO bardzo dziwne sa te tablice...0-d array nie mozna indeksowac
        #import pdb
        #pdb.set_trace()
        rc_an_l.append(rc_an - rc_0)
        rc_xx_l.append(rc_xx - rc_0)
        rc_2xx_l.append(rc_2xx - rc_0)
        rc_wrf_l.append(rc_wrf - rc_0)
        r_rsat_an_l.append(r_rsat_an)
    print "rc_xx po", rc_xx_l,"\n", "rc_erf", rc_wrf_l, "\n", "rc_an", rc_an_l
    return r_rsat_an_l, rc_xx_l, rc_2xx_l, rc_wrf_l, rc_an_l 
#    for key, value in expected.items():
#        print "\n key, valuu, eval(key)", key, value, eval(key)
#        assert abs(eval(key) - value) <= epsilon * abs(value)


def figure_plot(temp_list, libname_1m, libname_2m, libname_wrf, 
                rc_0, sup_lim, sup_step, plotname):    
    fig = plt.figure(1, figsize = (14.5,6))
    for nr, temp in enumerate(temp_list):
        ax = plt.subplot(1, len(temp_list), nr+1)
        r_rsat, c_xx, c_2xx, c_wrf, c_an = cond_all(libname_1m, libname_2m, 
                                                    libname_wrf, rc_0,
                                                    sup_lim, sup_step,
                                                    temp_0 = temp, press_0 = 900.e2)

        fig = plt.figure(1, figsize = (6.5,5.5))
        plt.plot(r_rsat, c_an, "b", r_rsat, c_xx, "r", r_rsat, c_xx, "r .", 
                 r_rsat, c_wrf, "g",
                 r_rsat, len(r_rsat)*[0], 'k--', len(c_an)*[0], c_an, 'k--' )
        plt.legend(["analytic","libcloudphxx_1m", "libcloudphxx_2m", "wrf_kessler"],
                   loc=2, prop = FontProperties(size=10))
        plt.title(str(temp),  fontsize=14)
        if nr == 0:
            plt.ylabel(r'$\Delta rc$', fontsize=14)
        plt.xlabel(r'$rv_0 - rv_{sat0}$', fontsize=14)
        ax.set_xlim(-1.*sup_lim,  sup_lim)
        ax.set_ylim(min(c_an), max(c_an))

        for item in plt.xticks()[1] + plt.yticks()[1]:
            item.set_fontsize(10)

    plt.savefig(plotname + ".pdf")
    plt.show()

def figure_plot_all(temp_list):
    Libname_1m = "libcloudphxx_blk_1m_pytest"
    Libname_2m = "libcloudphxx_blk_2m_pytest"
    Libname_wrf = "wrfkessler_blk_1m_pytest" #TODO: nie widzi libkessler.so (musialam przekopiowac)

    figure_plot(temp_list, 
                libname_1m=Libname_1m, libname_2m=Libname_2m, libname_wrf=Libname_wrf,
                rc_0=5.e-4, sup_lim=.6e-3, sup_step=0.1e-3, plotname="cond_evap_comp")


figure_plot_all([278.15, 283.15, 288.15])
