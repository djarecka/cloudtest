import sys
sys.path.append(".")
sys.path.append("../")
import numpy as np
from constants_pytest import Rd, Rv, cp, p0 #TODO                                        
import libcloudphxx as libcl
import analytic_blk_1m_pytest as an
from drops_py import distros

#TODO: should be an argument?
kappa = .61
n_tot = [60e6, 40e6]
meanr = [.04e-6 / 2, .15e-6 / 2]
gstdv = [1.4, 1.6]


def adj_cellwise(press, T, rv, dt, sstp=1):
    opts_init = libcl.lgrngn.opts_init_t()
    opts_init.dry_distros =  {kappa : distros.lognormal(n_tot, meanr, gstdv)}
    opts_init.dt = dt

    backend = libcl.lgrngn.backend_t.serial
    prtcls = libcl.lgrngn.factory(backend, opts_init)

    opts = libcl.lgrngn.opts_t()
    opts.adve = False
    opts.sedi = False
    opts.cond = True
    opts.coal = False
    opts.sstp_cond = sstp

    rho_d = np.array(an.density_dry(rv, press, T))
    th_d = np.array(an.pottemp_dry(rv, press, T))

    print "before init", rv, th_d
    prtcls.init(th_d, rv, rho_d)
    prtcls.step_sync(opts, th_d, rv, rho_d)
    print "after sync", rv, th_d
    return rv

