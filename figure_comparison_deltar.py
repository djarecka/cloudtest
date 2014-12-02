import sys
sys.path.append(".")
sys.path.append("wrf_microphys/kessler/")
sys.path.append("wrf_microphys/morrison_2mom/")
sys.path.append("libcloudphxx/")
import inspect
import pdb
import numpy as np
import analytic_blk_1m_pytest as an
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import importlib

Libs_comp = {
             "wrf_kessler":["wrfkessler_blk_1m_pytest", "1mom", 'g'], 
             "wrf_morrison":["wrfmorrison_blk_2m_pytest", "2mom", 'g.'],
             "libcloudphxx_1m":["libcloudphxx_blk_1m_pytest", "1mom", 'r'], 
             "libcloudphxx_2m":["libcloudphxx_blk_2m_pytest", "2mom", 'r.'],
             "libcloudph_lgr" : ["libcloudphxx_lgr_pytest", "lgr", 'b*']
             }

# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([5.e-4])
rr_0   = np.array([0.  ])
nc_0   = np.array([1.e8]) #TODO
nr_0   = np.array([1.e7])#TODO
dt_0   = 16.
sstp = 8


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
    if lib_type == "1mom":
        rv, rc, rr = lib_adj.adj_cellwise(press, T, rv, rc, rr, dt)
    elif lib_type == "2mom":
        rv, rc, nc, rr, nr = lib_adj.adj_cellwise(press, T, rv, rc, nc, rr, nr, dt)
    elif lib_type == "lgr":
        rv0 = rv.copy()
        print "rv0 before", rv0
        rv = lib_adj.adj_cellwise(press, T, rv, dt, sstp)
        print "rv0 after", rv0
        rc = rv0 - rv + rc
    return rv, rc

def figure_plot(ax, temp, r_rsat, c_lib, color, xy_lim=(None, None)):
    #pdb.set_trace()
    plt.plot(r_rsat, c_lib, color, markersize=10)
    plt.title("T="+str(temp)+",   dt="+str(dt_0)+",   substeps="+str(sstp), fontsize=14)
    plt.ylabel(r'$\Delta rc$', fontsize=14)
    plt.xlabel(r'$rv_0 - rv_{sat0}$', fontsize=14)
    ax.set_xlim(0, xy_lim[0])
    ax.set_ylim(0, xy_lim[1])

    for item in plt.xticks()[1] + plt.yticks()[1]:
        item.set_fontsize(10)


def main(temp_0, press_0, libs_comp_d, rc_0=.1e-4, sup_lim=1.e-3, sup_step=0.1e-3):
    fig = plt.figure(1, figsize = (9.5,6))
    ax = plt.subplot(1, 1, 1)
    legend_l = []
    
    # plotting the analytical solution
    r_rsat_an_l = []
    rc_an_l = []
    rsat_an = an.mixrat_sat(temp_0, press_0)
    for rv_0 in np.arange(rsat_an, rsat_an + sup_lim + sup_step, sup_step):
        r_rsat_an = rv_0 - rsat_an
        r_rsat_an_l.append(r_rsat_an)
        delta_an = an.delta_r(rv_0, temp_0, press_0)
        rv_an, rc_an = rv_0 - delta_an, rc_0 + delta_an
        rc_an_l.append(rc_an - rc_0)
    
    figure_plot(ax, temp_0, r_rsat_an_l, rc_an_l, color='k--', 
                xy_lim=(sup_lim, max(rc_an_l)))
    legend_l.append("analytic")

    # plotting the various model results
    for name, lib_info in libs_comp_d.iteritems():
        rc_lib_l = []
        for rv_0 in np.arange(rsat_an, rsat_an + sup_lim + sup_step, sup_step):
            rv_lib, rc_lib = condensation(lib=lib_info[0], lib_type=lib_info[1], 
                                          T=np.array(temp_0), press = np.array(press_0),
                                          rv=np.array(rv_0), rc=np.array(rc_0))

            rc_lib_l.append(rc_lib - rc_0)
        figure_plot(ax, temp_0, r_rsat_an_l, rc_lib_l, color=lib_info[2])
        legend_l.append(name)
    plt.legend(legend_l, loc=2, prop = FontProperties(size=10)) 
    plt.savefig("porownanie_time="+str(dt_0)+"_sstp="+str(sstp)+".pdf")
    plt.show()



                     
main(temp_0=273.15, press_0= 900.e2, libs_comp_d = Libs_comp)
