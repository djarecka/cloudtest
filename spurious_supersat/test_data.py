import sys, glob, os
import subprocess
sys.path.insert(0, "./")
import numpy as np
import pytest
import pdb
import json

import adv_hor as ah
import oop_adv_hor as oah

Aerosol={"meanr":.02e-6, "gstdv":1.4, "n_tot":1000e6,
         # ammonium sulphate aerosol parameter:
         "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
         "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
         "sd_conc":128 
}
Dt, Nt, It_output_l, Nx, Dx, Sl_sg = 0.2, 51, [50], 300, 2, slice(50,100)


@pytest.mark.parametrize("arg", [
   {"scheme":"sd", "apr":"trad", "setup":"rhoconst", "pl_flag":False, "nx":Nx,
    "sl_sg":Sl_sg, "crnt":1., "dt":Dt, "nt":Nt, "it_output_l":It_output_l, "aerosol":Aerosol},
   {"scheme":"sd", "apr":"S_adv", "setup":"rhoconst", "pl_flag":False, "nx":Nx,
    "sl_sg":Sl_sg, "crnt":1., "dt":Dt, "nt":Nt, "it_output_l":It_output_l, "aerosol":Aerosol},
   {"scheme":"sd", "apr":"trad", "setup":"rhoconst", "pl_flag":False, "nx":Nx,
    "sl_sg":Sl_sg, "crnt":.2, "dt":Dt, "nt":Nt, "it_output_l":It_output_l, "aerosol":Aerosol},
   {"scheme":"sd", "apr":"S_adv", "setup":"rhoconst", "pl_flag":False, "nx":Nx,
    "sl_sg":Sl_sg, "crnt":.2, "dt":Dt, "nt":Nt, "it_output_l":It_output_l, "aerosol":Aerosol}
    
  ])
def test_data_compare(arg, eps = 0.05):
    filename = arg["scheme"]+"_"+arg["apr"]+"_"+arg["setup"]+"_C"+str(arg["crnt"])+"_data_"+"10s.txt"

    ah.main(**arg)
    
    f_test = open(filename)
    data_test = json.load(f_test)

    f_ref = open("ref_data/"+filename)
    data_ref = json.load(f_ref)

    for key in data_ref:
        assert np.isclose(np.array(data_test[key]), np.array(data_ref[key]), atol=0, rtol=eps).all()



@pytest.mark.parametrize("arg", [
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"trad", "C":1., "dt":Dt, "time_adv_tot":Nt*Dt,
      "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "dirname_pre":"pytest"},
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"S_adv", "C":1., "dt":Dt, "time_adv_tot":Nt*Dt,
     "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "dirname_pre":"pytest"},
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"trad", "C":.2, "dt":Dt, "time_adv_tot":Nt*Dt,
     "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "dirname_pre":"pytest"},
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"S_adv", "C":.2, "dt":Dt, "time_adv_tot":Nt*Dt,
     "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "dirname_pre":"pytest"}
])
def test_data_compare_oop(arg, eps = 0.05):
    filename_ref = arg["scheme"]+"_"+arg["apr"]+"_"+arg["setup"]+"_C"+str(arg["C"])+"_data_"+str(int((Nt-1)*Dt))+"s.txt"
        
    dir_name = "output/pytest_scheme=sd_conc="+str(int(Aerosol["sd_conc"]))+"_setup="+arg["setup"]+"_apr="+arg["apr"]+"_dt="+str(Dt)+"_dx="+str(Dx)+"_C="+str(arg["C"]) + "_ntot="+str(int(Aerosol["n_tot"]/1.e6))
    filename=os.path.join(dir_name, "it="+str(int((Nt-1)*Dt))+"s.txt")
    ss = oah.Superdroplet(**arg)
    ss.all_sym()
                    

    f_test = open(filename)
    data_test = json.load(f_test)
    
    f_ref = open("ref_data/"+filename_ref)
    data_ref = json.load(f_ref)


    for key in data_ref:
        if key in ('rv', 'sd', "th_d", "na", "nc"):
            #pdb.set_trace()
            assert np.isclose(np.array(data_test[key]), np.array(data_ref[key]), atol=0, rtol=eps).all()
                                    
