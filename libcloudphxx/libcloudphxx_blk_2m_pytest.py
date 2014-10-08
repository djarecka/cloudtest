import sys
sys.path.append(".")

from numpy import array as arr_t
from constants_pytest import Rd, Rv, cp, p0
from libcloudphxx import blk_2m
import analytic_blk_1m_pytest as anpy

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

#TODO: should be an argument?
distr = [{"mean_rd":.04e-6 / 2, "sdev_rd":1.4, "N_stp":60e6, "chem_b":.55},
         {"mean_rd":.15e-6 / 2, "sdev_rd":1.6, "N_stp":40e6, "chem_b":.55}]


#TODO: ta f-cja jest tylko po to, aby byly keword arg., konieczna?
def adj_cellwise(press, T, rv, rc, nc, rr, nr, dt):
    opts = blk_2m.opts_t()
    opts.dry_distros = distr

    #TODO: should be always zero?
    dot_th = arr_t([0.])
    dot_rv = arr_t([0.])
    dot_rc = arr_t([0.])
    dot_nc = arr_t([0.])
    dot_rr = arr_t([0.])
    dot_nr = arr_t([0.])

    rho_d = arr_t(anpy.density_dry(rv, press, T))
    th_d = arr_t(anpy.pottemp_dry(rv, press, T))
    print "libcloudphxx_blk_2m, before calling", rv, rc, dot_nc
    blk_2m.rhs_cellwise(opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
                        rho_d, th_d, rv, rc, nc, rr, nr, dt)
    print "libcloudphxx_blk_2m, after calling", rv, rc, dot_rc, dot_nc, dot_th
    return rv + dot_rv, rc + dot_rc, nc + dot_nc, rr + dot_rr, nr + dot_nr

