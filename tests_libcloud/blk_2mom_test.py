import sys
sys.path.append(".")
sys.path.append("../../")
sys.path.append("wrf_microphys/morrison_2momNc/")
import pytest
import inspect
import pdb
import numpy as np
import analytic_blk_1m_pytest as an_blk

# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([5.e-4])
nc_0   = np.array([0.  ])
rr_0   = np.array([0.  ])
nr_0   = np.array([0.  ])
dt_0   = 2

def condensation(lib, press = None, T = None, rv = None,
                 rc = None, nc = None, rr = None, nr = None, dt=dt_0):

    import importlib
    lib_adj = importlib.import_module(lib)

    print "\n initial value of qv", rv
    #TODO: in some loop?
    T = T if T!=None else T_0.copy()
    rv = rv if rv!=None else rv_0.copy()
    rc = rc if rc!=None else rc_0.copy()
    nc = nc if nc!=None else nc_0.copy()
    rr = rr if rr!=None else rr_0.copy()
    nr = nr if nr!=None else nr_0.copy()
    press = press if press!=None else press_0.copy()
    print "\n In condensation. Who is calling..?", inspect.stack()[1][3]
    rv, rc, nc, rr, nr = lib_adj.adj_cellwise(press, T, rv, rc, nc, rr, nr, dt)
    return rv, rc

def analytic_condensation(press = None, T = None, rv = None, rc = None):
     T = T if T!=None else T_0.copy()
     rv = rv if rv!=None else rv_0.copy()
     rc = rc if rc!=None else rc_0.copy()
     press = press if press!=None else press_0.copy()
     r_an = {}
     #calling analytical function calculating del_r for supersaturation adjustment
     del_r = an_blk.delta_r(rv, rc, T, press)
     r_an["rv"] = rv - del_r
     r_an["rc"] = rc + del_r
     return r_an

#@pytest.mark.skipif
@pytest.mark.parametrize("arg", [
    {'T':np.array([255.])},    pytest.mark.xfail({'T':np.array([500. ])}),
    {'rv':np.array([-1.e-5])}, pytest.mark.xfail({'rv':np.array([0.1 ])}),
    {'rc':np.array([-1.e-5])}, pytest.mark.xfail({'rc':np.array([0.01])}),
    {'rr':np.array([-1.e-5])}, pytest.mark.xfail({'rr':np.array([0.01])})
    ])
def test_exeptions_wrongvalue(libname, arg):
    print "\n in test_exeption", arg
    with pytest.raises(Exception) as excinfo:
        condensation(lib=libname, **arg)
    #the line below can give you information about the exception 
    #print "exception info:", excinfo.value
    
#@pytest.mark.skipif
@pytest.mark.parametrize("arg", [
        ({"rv" :  np.array([0.]),      "rc" :  np.array([0.])}), # no water
        ({"rv" :  np.array([7.e-3]),   "rc" :  np.array([0.])}), # no cl water and subsat.
        ({"rv" :  np.array([11.5e-3]), "rc" :  np.array([0.])}), # no cl water and supersat.
        ({"rv" :  np.array([5.e-3]), "rc" :  np.array([1.e-3]), "nc" :  np.array([1.e8])}), # subsat. leads to coplete evap.
        ({"rv" :  np.array([8.e-3]), "rc" :  np.array([1.e-3]), "nc" :  np.array([1.e8])}), # subsat. leads to some evap.
        ({"rv" :  np.array([9.e-3]), "rc" :  np.array([1.e-3]), "nc" :  np.array([1.e8])}), # supersat. leads to cond.
        ({"rv" :  np.array([0.005077]), "rc" :  np.array([1.1958382416435485e-08]), "nc" :  np.array([241337.]), "T" : np.array([284.1]), "dt" : 0.1 }), #assert error 
])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(libname, arg, epsilon = 0.15):
    print "\n in test_expected value before cond.", arg
    r_an = analytic_condensation(rv=arg["rv"], rc=arg["rc"])
    rv, rc = condensation(lib=libname, **arg)
    #print "rv, rc po", rv, rc
    for key in ["rv", "rc"]:
        print "\n key, valuu, eval(key)", key, arg[key], r_an[key]
        assert abs(arg[key] - r_an[key]) <= epsilon * abs(r_an[key])
            
