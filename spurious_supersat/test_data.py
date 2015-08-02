import sys, glob, os
import subprocess
sys.path.insert(0, "./")
import numpy as np
import pytest
import pdb
import json

import adv_hor_wh as ah
#import adv_hor as ah

@pytest.mark.parametrize("arg", [
  {"scheme":"2m", "apr":"trad",  "setup":"rhoconst", "nt":501, "outfreq":500, "dt":.2},
  {"scheme":"2m", "apr":"S_adv", "setup":"rhoconst", "nt":501, "outfreq":500, "dt":.2},
  {"scheme":"sd", "apr":"trad",  "setup":"rhoconst", "nt":501, "outfreq":500, "dt":.2},
  {"scheme":"sd", "apr":"S_adv", "setup":"rhoconst", "nt":501, "outfreq":500, "dt":.2}
  ])
def test_data_compare(arg, eps = 0.1):
#    filename = arg["scheme"]+"_"+arg["apr"]+"_"+arg["setup"]+"_"+"data_init.txt" 
    filename = arg["scheme"]+"_"+arg["apr"]+"_"+arg["setup"]+"_"+"data_"+str(int((arg["nt"]-1)*arg["dt"]))+"s.txt"

    ah.main(**arg)

    f_test = open(filename)
    data_test = json.load(f_test)

    f_ref = open("ref_data/"+filename)
    data_ref = json.load(f_ref)
   

    for key in data_ref:
        assert np.isclose(np.array(data_test[key]), np.array(data_ref[key]), atol=0, rtol=eps).all()
