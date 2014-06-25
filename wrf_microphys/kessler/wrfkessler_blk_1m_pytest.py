"""prepares kessler-python function for tests """

import numpy as np
from cffi_kessler import kessler
from stale import Rd, cp #TODO ujednolicic stale
import pdb

#TODO: general constants
Rv = 461.4
p0 = 1000.e2

nx = 1
ny = 1
nz = 1

dz8w = np.ones((nx,nz,ny)) * 20.
z = np.ones((nx,nz,ny)) * 700.
rainnc = np.zeros((nx,ny))
rainncv = np.zeros((nx,ny))

def exner(press):
    return np.array((press/1.e5)**(Rd/cp))

def density(rv, press, T):
    #TODO check, if it's full rho
    R_tot = (Rd + rv*Rv) / (1. + rv)
    rho = press / R_tot / T
    return np.array(rho)

def pottemp(press, T):
     theta = T * (p0 / press)**(Rd/cp) 
     return np.array(theta)

def adj_cellwise(press_in, T_in, qv_in, qc_in, qr_in, dt):
    pii = exner(press_in)
    rho = density(qv_in, press_in, T_in)
    th = pottemp(press_in, T_in)
    kessler(nx, ny, nz, dt,
            th, qv_in, qc_in, qr_in, rho, pii, dz8w, z,
            rainnc, rainncv)
    return qv_in, qc_in

