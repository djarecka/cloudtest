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


def libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt):
    opts = libcl.blk_2m.opts_t()

    dot_th = np.zeros((100,))
    dot_rv = np.zeros((100,))
    dot_rc = np.zeros((100,))
    dot_nc = np.zeros((100,))
    dot_rr = np.zeros((100,))
    dot_nr = np.zeros((100,))

    
    print "th przed mikro", th_d
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, th_d, rv, rc, nc, rr, nr, dt)
    print "th po mikro", th_d + dot_th*dt # daje to samo - cos jest zle!
     
    th_d = th_d + dot_th * dt
    rv   = rv   + dot_rv * dt
    rc   = rc   + dot_rc * dt
    nc   = nc   + dot_nc * dt
    rr   = rr   + dot_rr * dt
    nr   = nr   + dot_nr * dt
    print "max qc po mikro", rc.max()
    return th_d, rv, rc, nc, rr, nr

#libcl_2mom()

def libcl_1mom(rho_d, th_d, rv, rc, rr, dt):
    opts = libcl.blk_1m.opts_t()

    libcl.blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    

#libcl_1mom()

def main(scheme, nx=100):
    th_d = np.ones((nx,))* 287.
    rv = np.ones((nx,))* 2.e-3
    rc = np.zeros((nx,))
    rc[10:20] = 1.e-3
    rr = np.zeros((nx,))
    rho_d = np.ones((nx,))
    testowa = np.zeros((nx,))
    testowa[10:20]= 1.
    
    var_adv = [th_d, rv, rc, rr, testowa]

    if scheme == "2m":
           nc = np.ones((nx,))*1.e8
           nr = np.zeros((nx,))
           var_adv  = var_adv + [nc, nr]

    dt = 1.

    plotting(rc)
    for it in range(100):
        print "it", it
        #print "qc przed adv", rc
        for var in var_adv:
            libmpdata.mpdata(var, .5, 1);
        #print "qc po adv", rc

        if scheme == "1m":
            libcl_1mom(rho_d, th_d, rv, rc, rr, dt)

            #TODO spr dlaczego musze przekazywac!!!
        if scheme == "2m":
            th_d, rv, rc, nc, rr, nr = libcl_2mom(rho_d, th_d, rv, rc, rr, nc, nr, dt)


        #print " qc po mikro", rv

        print "testowa po it = ", it
        if it/10 * 10 == it:
            plotting(rc)



#main("2m") #dziala, ale nie przesuwa
#main("2m", nx=50) # nie dziala! mowi, ze th<0!
main("1m", nx=50) # dziala (tylko, ze periodycznosci nie ma, ale to mniejsza z tym)
#main("1m") # dziala
