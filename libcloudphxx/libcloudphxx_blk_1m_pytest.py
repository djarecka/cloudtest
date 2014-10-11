import sys
sys.path.append(".")
sys.path.append("../")
import numpy as np
from constants_pytest import Rd, Rv, cp, p0 #TODO
import libcloudphxx as libcl
import analytic_blk_1m_pytest as an

def opts_cr(cond = True, cevp = True, revp = True, conv = True,
            accr = True, sedi = False):
    opts = libcl.blk_1m.opts_t()
    opts.cond = cond
    opts.cevp = cevp
    opts.revp = revp
    opts.conv = conv
    opts.accr = accr
    opts.sedi = sedi
    return opts

#TODO: think about opts_cr, now it can't be changed calling adj_cellwise
def adj_cellwise(press, T, rv, rc, rr, dt):
    opts = opts_cr()
    rho_d = np.array(an.density_dry(rv, press, T))
    th_d = np.array(an.pottemp_dry(rv, press, T))
    libcl.blk_1m.adj_cellwise(opts, rho_d, th_d, rv, rc, rr, dt)
    return rv, rc

