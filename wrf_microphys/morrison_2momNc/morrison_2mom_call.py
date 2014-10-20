import numpy as np
from cffi_morrison_2mom import  morrison_2mom_simplewarm

nx = 1
ny = 1
nz = 1
dt_in = 1.

temp_np = np.ones((nx,nz,ny)) * 283.15
qv_np = np.ones((nx,nz,ny)) * 10.e-3
qc_np = np.ones((nx,nz,ny)) * 0.e-3
nr_np = np.ones((nx,nz,ny)) * 0.e-3
qr_np = np.ones((nx,nz,ny)) * 0.e-3
press_np = np.ones((nx,nz,ny)) * 900.e2
dz_np = np.ones((nx,nz,ny)) * 20.

print "w call przed", qv_np, qc_np, temp_np, qr_np
morrison_2mom_simplewarm(nx, ny, nz, qc_np, qr_np, nr_np,
                         temp_np, qv_np, press_np, dz_np, dt_in)
print "w call po", qv_np, qc_np, temp_np, qr_np
