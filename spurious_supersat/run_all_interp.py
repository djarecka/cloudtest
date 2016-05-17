import os
from subprocess import call
from PyPDF2 import PdfFileReader, PdfFileMerger
import numpy as np
from oop_sdrop_adv_hor import Superdroplet
import json
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import pdb

# test for interpolation

Arg = {"nx":300, "dx":2, "sl_sg":slice(50,100), "dt":.1, "time_adv_tot":26,
       "aerosol":{"meanr":.02e-6, "gstdv":1.4, "n_tot":1e9, "chem_b":.505, "kappa":.61, "sd_conc":256},
              "sl_act_time":60, "setup":"rhoconst", "apr":"trad", "RHenv":0.95,
       "test":False, "it_output_l":[249],#, 499],
       "dirname_pre":"Interp"}


def simulat(numb_intrp, arg, Crnt):
    dirname_full = "plots/"+arg["dirname_pre"]+"_scheme=sd_conc="+str(int(arg["aerosol"]["sd_conc"]))+"_setup="+arg["setup"]+"_apr="+arg["apr"]+"_dt="+str(arg["dt"])+"_dx="+str(arg["dx"])+"_C="+str(Crnt) + "_ntot="+str(int(arg["aerosol"]["n_tot"]/1.e6))
    if arg["setup"]=="slow_act":
        dirname_full += "_sl_act_time=" + str(arg["sl_act_time"]) + "s"
    if numb_intrp > 1:
        dirname_full += "_nintrp=" + str(numb_intrp)
                                    
    ss  = Superdroplet(C=Crnt, n_intrp=numb_intrp, **arg)
    ss.all_sym()
    return dirname_full


def main(Crnt, arg=Arg, num_intrp_l = [1, 2, 4, 8]):
    dirname_l = []
    for num in num_intrp_l:
        dirname = simulat(num, arg, Crnt)
        dirname_l.append(dirname)

    dirname_allplots = os.path.join("plots", "testinterpolation_C="+str(Crnt)+"_num_interp="+str(num_intrp_l))
    if os.path.exists(dirname_allplots):
        call(["rm", "-r", dirname_allplots])
    call(["mkdir", dirname_allplots])
                                            
    # joining files from all directories for every timestep
    for it in arg["it_output_l"]:
        time = str((it+1)*arg["dt"])
        filename_l = []
        for dir in dirname_l:
            filename_l.append(os.path.join(dir, "profiles_time="+time+"s.pdf"))
            
        merger = PdfFileMerger()

        for filename in filename_l:
            merger.append(PdfFileReader(filename, "rb"))

        merger.write(os.path.join(dirname_allplots, "all_plots_time="+time+".pdf"))

    plotting_comparison(arg, dirname_l, Crnt, num_intrp_l)

def plotting_comparison(arg, dirname_l, Crnt, num_intrp_l):
    # comparison on one plot
    dirname_allplots = os.path.join("plots", "testinterpolation_C="+str(Crnt)+"_num_interp="+str(num_intrp_l))
    it_last = arg["it_output_l"][-1]
    nrow, ncol = 2, 2
    fig, tpl = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12,6))
    Colors = ["k", "r", "g", "b"]
    legend_l = []
    for nn, dirname in enumerate(dirname_l):
        #pdb.set_trace()
        legend_l.append("interp_num="+ str(num_intrp_l[nn]))
        dirname_o = dirname.replace("plots","output")
        f_r = open(os.path.join(dirname_o, "it="+str(int(arg["dt"]*(it_last+1)))+"s.txt"))
        dict_r = json.load(f_r)
        i=0
        for k, var in enumerate(["rv", "rc", "th_d", "S"]):
            tpl[i%nrow,i/nrow].set_title(var+", ")
            #if k in ylim_dic.keys():
            #    tpl[i%nrow,i/nrow].set_ylim(ylim_dic[k])
            tpl[i%nrow,i/nrow].plot(dict_r[var], Colors[nn])
            i+=1
    plt.legend(legend_l, prop = FontProperties(size=10)) 
    plt.savefig(os.path.join(dirname_allplots, "compar_plots_time="+str((it_last+1)*arg["dt"])+".pdf"))
    plt.show()
            

        
main(Crnt=0.)
main(Crnt=0.2)

