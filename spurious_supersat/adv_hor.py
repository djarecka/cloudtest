import sys
sys.path.append(".")
sys.path.append("../")
import libmpdata
import numpy as np
import libcloudphxx as libcl
#import analytic_blk_1m_pytest as an
#from constants_pytest import Rd, Rv, cp, p0
import pdb

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


def test():
    a = np.zeros((10,))
    a[1]=1
    print a
#.5 - l. courant, 3-liczba krokow                                                              
    libmpdata.mpdata(a, 1., 3);
#assert a[1] < 1 and a[2] > 0 #SA - po co ten assert?                                          
    print a


#test()

def plotting(prof):

    plt.figure(1)
    plt.plot(prof)
    plt.show()


def libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt, nx):
    opts = libcl.blk_2m.opts_t()
    opts.acti = False
    opts.cond = True
    opts.acnv = False
    opts.accr = False
    opts.sedi = False

    dot_th = np.zeros((nx,))
    dot_rv = np.zeros((nx,))
    dot_rc = np.zeros((nx,))
    dot_nc = np.zeros((nx,))
    dot_rr = np.zeros((nx,))
    dot_nr = np.zeros((nx,))

    
    print "nc min, max przed mikro", nc.min(), nc.max()
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, th_d, rv, rc, nc, rr, nr, dt)
        
    th_d += dot_th * dt
    rv   += dot_rv * dt
    rc   += dot_rc * dt
    nc   += dot_nc * dt
    rr   += dot_rr * dt
    nr   += dot_nr * dt
    print "nc min, max po mikro", nc.min(), nc.max()



def libcl_1mom(rho_d, th_d, rv, rc, rr, dt):
    opts = libcl.blk_1m.opts_t()
    opts.cond = True
    opts.cevp = True
    opts.revp = False
    opts.conv = False
    opts.accr = False
    opts.sedi = False

    print "1m rc max, min przed mikro", rc.max(), rc.min()
    libcl.blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    print "1m rc max, min po mikro", rc.max(), rc.min()



def main(scheme, nx=300, sl_sg = slice(50,100), crnt=0.1, dt=0.2):
    th_d = np.ones((nx,))* 287.
    rv = np.ones((nx,))* 2.e-3
    rc = np.zeros((nx,))
    rc[sl_sg] = 1.e-3
    rr = np.zeros((nx,))
    rho_d = np.ones((nx,))
    testowa = np.zeros((nx,))
    testowa[sl_sg]= 1.e4
    
    var_adv = [th_d, rv, rc, rr, testowa]

    if scheme == "2m":
           nc = np.zeros((nx,))
           nc[sl_sg] = 1.e8
           nr = np.zeros((nx,))
           var_adv  = var_adv + [nc, nr]

    for it in range(500):
        print "it", it

        print "testowa min, max przed adv", testowa.min(), testowa.max()
        if scheme == "2m": print "nc min, max przed adv", nc.min(), nc.max()        
        for var in var_adv:
            libmpdata.mpdata(var, crnt, 1);
        if scheme == "2m": print "nc min, max po adv", nc.min(), nc.max()
        print "testowa min, max po adv", testowa.min(), testowa.max()

        if scheme == "1m":
            libcl_1mom(rho_d, th_d, rv, rc, rr, dt)

        if scheme == "2m":
            libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt, nx)
        

        print "testowa po it = ", it
        if it/100 * 100 == it:
            plotting(nc)



main("2m") 
#main("1m") 
