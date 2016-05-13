import sys
sys.path.append(".")
sys.path.append("wrf_microphys/kessler/")
import pytest
import inspect
import pdb
import numpy as np

# typical values as an example
press_0 = np.array([900.e2  ])
th_0   = np.array([291.8])
T_0    = np.array([283.15])
rv_0   = np.array([8.e-3])
rc_0   = np.array([0. ])
rr_0   = np.array([0.  ])
dt_0   = 20
sstp_0 = 10

#works for various library, during calling has to be added --libname=....
def condensation(lib, press = None, T = None,
                 rv = None, rc = None, rr = None, dt=dt_0, sstp=sstp_0):

    import importlib
    lib_adj = importlib.import_module(lib)
    
    print "\n initial value of qv", rv
    #pdb.set_trace()
    T = T if T!=None else T_0.copy()
    rv = rv if rv!=None else rv_0.copy()
    rc = rc if rc!=None else rc_0.copy()
    rr = rr if rr!=None else rr_0.copy()
    press = press if press!=None else press_0.copy()
    print "\n In condensation. Who is calling..?", inspect.stack()[1][3]
    rv0 = rv.copy() 
    rv = lib_adj.adj_cellwise(press, T, rv, dt, sstp)
    rc = rv0 - rv + rc
    return rv, rc

#@pytest.mark.skipif
@pytest.mark.parametrize("arg", [
    {'T':np.array([255.])},    pytest.mark.xfail({'T':np.array([500. ])}),
    {'rv':np.array([-1.e-5])}, pytest.mark.xfail({'rv':np.array([0.1 ])}),
    ])
def test_exeptions_wrongvalue(libname, arg):
    print "\n in test_exeption", arg
    with pytest.raises(Exception) as excinfo:
        condensation(lib=libname, **arg)
    #the line below can give you information about the exception 
    #print "exception info:", excinfo.value
    

#@pytest.mark.skipif
@pytest.mark.parametrize("arg, expected", [
        ({"rv" : np.array([0.]),     "rc" :  np.array([0.])},
         {"rv" : np.array([0.]),     "rc" :  np.array([0.])}), # no water
        ({"rv" :  np.array([7.e-3]),  "rc" :  np.array([0.])},
         {"rv" :  np.array([7.e-3]),  "rc" :  np.array([0.])}), # no cl water and subsat.
        ({"rv" :  np.array([10.e-3]), "rc" :  np.array([0.])},
         {"rv" :  np.array([9.44e-3]), "rc" : np.array([.56e-3])}), # no cl water and supersat.
    ])
#TODO zastanowic sie nad epsilonem
def test_expected_output_evapcond(libname, arg, expected, epsilon = 0.2):
    print "\n in test_expected value before cond.", arg
    rv, rc = condensation(lib=libname, **arg)
    #print "rv, rc po", rv, rc
    for key, value in expected.items():
        print "\n key, valuu, eval(key)", key, value, eval(key)
        assert abs(eval(key) - value) <= epsilon * abs(value)

