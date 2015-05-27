import sys
sys.path.append(".")
sys.path.append("wrf_microphys/morrison_2momNc/")
import pytest
import inspect
import pdb
import numpy as np

# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([5.e-4])
nc_0   = np.array([0.  ])
rr_0   = np.array([0.  ])
nr_0   = np.array([0.  ])
dt_0   = 1

def condensation(lib, press = None, T = None, rv = None,
                 rc = None, nc = None, rr = None, nr = None, dt=dt_0):

    import importlib
    lib_adj = importlib.import_module(lib)

    print "\n initial value of qv", rv
    #pdb.set_trace()
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


@pytest.mark.skipif
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
    

#TODO: wypisac znane outputy, moze tez theta?
#TODO: polaczyc z plikiem analytic_blk_1m_pytest??
#@pytest.mark.skipif
@pytest.mark.parametrize("arg, expected", [
        ({"rv" :  np.array([0.]),      "rc" :  np.array([0.])},
         {"rv" :  np.array([0.]),      "rc" :  np.array([0.])}), # no water
        ({"rv" :  np.array([7.e-3]),   "rc" :  np.array([0.])},
         {"rv" :  np.array([7.e-3]),   "rc" :  np.array([0.])}), # no cl water and subsat.
        ({"rv" :  np.array([10.e-3]),  "rc" :  np.array([0.])},
         {"rv" :  np.array([9.44e-3]), "rc" :  np.array([.56e-3])}), # no cl water and supersat.
        ({"rv" :  np.array([5.e-3]),   "rc" :  np.array([1.e-3])},
         {"rv" :  np.array([6.e-3]),   "rc" :  np.array([0.])}), # subsat. leads to coplete evap.
        ({"rv" :  np.array([8.e-3]),   "rc" :  np.array([1.e-3])},
         {"rv" :  np.array([8.26e-3]), "rc" :  np.array([0.74e-3])}), # subsat. leads to some evap.
        ({"rv" :  np.array([9.e-3]),   "rc" :  np.array([1.e-3])},
         {"rv" :  np.array([8.85e-3]), "rc" :  np.array([1.15e-3])}), # supersat. leads to cond.
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(libname, arg, expected, epsilon = 0.1):
    print "\n in test_expected value before cond.", arg
    rv, rc = condensation(lib=libname, **arg)
    #print "rv, rc po", rv, rc
    for key, value in expected.items():
        print "\n key, valuu, eval(key)", key, value, eval(key)
        assert abs(eval(key) - value) <= epsilon * abs(value)

