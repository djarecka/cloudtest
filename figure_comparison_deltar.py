import sys
sys.path.append(".")
sys.path.append("wrf_microphys/kessler/")
sys.path.append("wrf_microphys/morrison_2mom/")
sys.path.append("libcloudphxx/")
#sys.path.append("/Users/dorota/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages")
#sys.path.append("/Users/dorota/libcloudphxx/build/tests/python")
import inspect
import pdb
import numpy as np
import analytic_blk_1m_pytest as an
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import importlib

# library used in comparison
Libxx_1m = "libcloudphxx_blk_1m_pytest"
Libxx_2m = "libcloudphxx_blk_2m_pytest"
Libwrf_kes = "wrfkessler_blk_1m_pytest" #TODO: nie widzi libkessler.so (musialam przekopiowac)               
Lib_2mor = "wrfmorrison_blk_2m_pytest"


# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([5.e-4])
rr_0   = np.array([0.  ])
nc_0   = np.array([1.e8]) #TODO
nr_0   = np.array([1.e7])#TODO
dt_0   = 1.

def condensation(lib, lib_type, press = None, T = None,
                 rv = None, rc = None, rr = None,
                 nc = None, nr = None, dt=dt_0):

    lib_adj = importlib.import_module(lib)

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

#condensation/evaporation if no rain
def cond_all(libname, libname_2m, libname_wrf, libname_2mor, rc_0, sup_lim, sup_step, 
             temp_0, press_0):
    rc_an_l = []
    rc_xx_l = []
    rc_2xx_l = []
    rc_wrf_l = []
    rc_2mor_l = []
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
        rv_2mor, rc_2mor = condensation(lib=libname_2mor, lib_type=2, T=np.array(temp_0),
                                      press = np.array(press_0),
                                      rv=np.array(rv_0), rc=np.array(rc_0))

        delta_an = an.delta_r(rv_0, temp_0, press_0) 
        rv_an, rc_an = rv_0 - delta_an, rc_0 + delta_an
        r_rsat_an = rv_0 - rsat_an

        rc_an_l.append(rc_an - rc_0)
        rc_xx_l.append(rc_xx - rc_0)
        rc_2xx_l.append(rc_2xx - rc_0)
        rc_wrf_l.append(rc_wrf - rc_0)
        rc_2mor_l.append(rc_2mor - rc_0)
        r_rsat_an_l.append(r_rsat_an)

    print "rc_xx po", rc_xx_l,"\n", "rc_2xx po", rc_2xx_l,"\n","rc_wrf", rc_wrf_l, "\n", "rc_2mor", rc_2mor_l, "\n", "rc_an", rc_an_l
    return r_rsat_an_l, rc_xx_l, rc_2xx_l, rc_wrf_l, rc_2mor_l, rc_an_l 


def figure_plot(ax, temp, libxx_1m, libxx_2m, libwrf_kes, lib_2mor, 
                rc_0, sup_lim, sup_step, label_y):    
    r_rsat, c_xx, c_2xx, c_wrf, c_2mor, c_an = cond_all(libxx_1m, libxx_2m, 
             libwrf_kes, lib_2mor, rc_0, sup_lim, sup_step,
             temp_0 = temp, press_0 = 900.e2)

    plt.plot(r_rsat, c_an, "b", r_rsat, c_xx, "r", r_rsat, c_2xx, "r .", 
             r_rsat, c_wrf, "g", r_rsat, c_2mor, "g.",
             r_rsat, len(r_rsat)*[0], 'k--', len(c_an)*[0], c_an, 'k--' )
    plt.legend(["analytic","libcloudphxx_1m", "libcloudphxx_2m", 
                "wrf_kessler", "morrison_2m",],
               loc=2, prop = FontProperties(size=10))
    plt.title(str(temp),  fontsize=14)
    if label_y:
        plt.ylabel(r'$\Delta rc$', fontsize=14)
    plt.xlabel(r'$rv_0 - rv_{sat0}$', fontsize=14)
    ax.set_xlim(-1.*sup_lim,  sup_lim)
    #ax.set_ylim(min(c_an), max(c_an))

    for item in plt.xticks()[1] + plt.yticks()[1]:
        item.set_fontsize(10)



def figure_plot_all(temp_list):

    fig = plt.figure(1, figsize = (14.5,6))
    label_y = True
    for nr, temp in enumerate(temp_list):
        ax = plt.subplot(1, len(temp_list), nr+1)
        if nr > 0:
            label_y=False
        figure_plot(ax, temp, libxx_1m=Libxx_1m, libxx_2m=Libxx_2m, 
                    libwrf_kes=Libwrf_kes, lib_2mor=Lib_2mor,
                    rc_0=5.e-4, sup_lim=.6e-3, sup_step=0.1e-3, 
                    label_y=label_y)

    plt.savefig("cond_evap_comp_nc1e8.pdf")
    plt.show()


figure_plot_all([278.15, 283.15, 288.15])
