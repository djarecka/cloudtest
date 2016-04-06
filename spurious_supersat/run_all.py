import os
from PyPDF2 import PdfFileReader, PdfFileMerger
import numpy as np
import oop_adv_hor as oah


Arg = {"nx":300, "dx":2, "sl_sg":slice(50,100), "C":.2, "dt":.1, "time_adv_tot":26,
       "aerosol":{"meanr":.02e-6, "gstdv":1.4, "n_tot":1e9, "chem_b":.505, "kappa":.61, "sd_conc":256},
       "sl_act_time":60, "n_intrp":1, "setup":"slow_act", "RHenv":0.95,
       "test":False, "it_output_l":[600, 850, 851, 1100,1101, 1350, 1351,1700,1701],
       "dirname_pre":"Sym"}

def main(aproach):
    dirname_full = "plots/"+Arg["dirname_pre"]+"_scheme=sd_conc="+str(int(Arg["aerosol"]["sd_conc"]))+"_setup="+Arg["setup"]+"_apr="+aproach+"_dt="+str(Arg["dt"])+"_dx="+str(Arg["dx"])+"_C="+str(Arg["C"]) + "_ntot="+str(int(Arg["aerosol"]["n_tot"]/1.e6))
    if Arg["setup"]=="slow_act":
        dirname_full += "_sl_act_time=" + str(Arg["sl_act_time"]) + "s"
    if Arg["n_intrp"]>1:
        dirname_full += "_nintrp=" + str(Arg["n_intrp"])
                                    
    ss  = oah.Superdroplet(apr=aproach, **Arg)
    ss.all_sym()

    files_dir = dirname_full
    pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
    merger = PdfFileMerger()

    for filename in pdf_files:
        merger.append(PdfFileReader(os.path.join(files_dir, filename), "rb"))

    merger.write(os.path.join(files_dir, "all_plots.pdf"))
        

main(aproach="trad")
main(aproach="S_adv")
