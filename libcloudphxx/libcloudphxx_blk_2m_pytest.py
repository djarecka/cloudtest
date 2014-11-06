import sys
sys.path.append(".")
sys.path.append("../")
import numpy as np
from constants_pytest import Rd, Rv, cp, p0 #TODO
import libcloudphxx as libcl
import analytic_blk_1m_pytest as an

#TODO: do I need it??
#def opts_cr(cond = True, cevp = True, revp = True, conv = True,
#            accr = True, sedi = False):
#    opts = blk_1m.opts_t()
#    opts.cond = cond
#    opts.cevp = cevp
#    opts.revp = revp
#    opts.conv = conv
#    opts.accr = accr
#    opts.sedi = sedi
#    return opts

#TODO: not needed for cond/evap!!!
#distr = [{"mean_rd":.04e-6 / 2, "sdev_rd":1.4, "N_stp":60e6, "chem_b":.55},
#         {"mean_rd":.15e-6 / 2, "sdev_rd":1.6, "N_stp":40e6, "chem_b":.55}]


def adj_cellwise(press, T, rv, rc, nc, rr, nr, dt):
    opts = libcl.blk_2m.opts_t()
#    opts.dry_distros = distr

    #TODO: should be always zero?
    dot_th = np.array([0.])
    dot_rv = np.array([0.])
    dot_rc = np.array([0.])
    dot_nc = np.array([0.])
    dot_rr = np.array([0.])
    dot_nr = np.array([0.])

    rho_d = np.array(an.density_dry(rv, press, T))
    th_d = np.array(an.pottemp_dry(rv, press, T))
    print "libcloudphxx_blk_2m, before calling", rv, rc, dot_nc, dt
    libcl.blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                              rho_d, th_d, rv, rc, nc, rr, nr, dt)
    print "libcloudphxx_blk_2m, after calling", rv, rc, dot_rc, dot_nc, dot_th
    return rv + dot_rv*dt, rc + dot_rc*dt, nc + dot_nc*dt, rr + dot_rr*dt, nr + dot_nr*dt

