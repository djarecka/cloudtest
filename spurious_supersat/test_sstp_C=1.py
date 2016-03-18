import sys, glob, os
import subprocess
sys.path.insert(0, "./")
import numpy as np
import pytest
import pdb
import json

import oop_adv_hor as oah

Aerosol={"meanr":.02e-6, "gstdv":1.4, "n_tot":1000e6,
         # ammonium sulphate aerosol parameter:
         "chem_b":.505, # blk_2m only (sect. 2 in Khvorosyanov & Curry 1999, JGR 104)
         "kappa":.61,    # lgrngn only (CCN-derived value from Table 1 in Petters and Kreidenweis 2007)
         "sd_conc":128 
}
Dt, Nt, It_output_l, Nx, Dx, Sl_sg = 0.2, 51, [20, 50], 300, 2, slice(50,100)


@pytest.mark.parametrize("arg", [
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"trad", "C":1., "dt":Dt, "time_adv_tot":Nt*Dt,
      "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "sstp_cond":5},
    {"nx":Nx, "dx":Dx, "sl_sg":Sl_sg, "apr":"S_adv", "C":1., "dt":Dt, "time_adv_tot":Nt*Dt,
     "aerosol":Aerosol, "sl_act_time":0, "n_intrp":1, "setup":"rhoconst", "scheme":"sd",
     "test":False, "it_output_l":It_output_l, "sstp_cond":5},
])
# checking if there are no oscilations in cloud water for courant=1
def test_sstp_courant1(tmpdir, arg, eps = 0.01):
    ss = oah.Superdroplet(dirname=str(tmpdir), **arg)
    ss.all_sym()

    for filename in os.listdir(str(tmpdir)):
        if filename.startswith("it"):
            f_test = open(os.path.join(str(tmpdir), filename))
            data_test = json.load(f_test)
            rc_array = np.array(data_test["rc"])
            rc_array_nozero = rc_array[np.where(rc_array>0)]
            assert rc_array_nozero.max() / rc_array_nozero.mean() <= 1. + eps
                                    
