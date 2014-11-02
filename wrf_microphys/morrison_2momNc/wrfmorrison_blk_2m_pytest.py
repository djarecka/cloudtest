import sys
sys.path.append(".")

import numpy as np
from cffi_morrison_2momNc import morrison_2momNc_simplewarm
import pdb

nx = 1
ny = 1
nz = 1
dz = np.ones((nx,nz,ny)) * 20.

def adj_cellwise(press_in, T_in, qv_in, qc_in, nc_in, qr_in, nr_in, dt):
    print "wrfmorrison qv before simplewarm", qv_in, qc_in, nc_in, press_in, T_in
    morrison_2momNc_simplewarm(nx, ny, nz, qc_in, qr_in, nc_in, nr_in,
                               T_in, qv_in, press_in, dz, dt)
    print "wrfmorrison qv after simplewarm", qv_in, qc_in, qr_in, T_in
    return qv_in, qc_in, nc_in, qr_in, nr_in
