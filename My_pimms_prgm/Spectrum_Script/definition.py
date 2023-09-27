import os
import numpy as np
import matplotlib.pyplot as plt 
def get_spectrum(T, Nh, filterr, band, step,header): #header=[choix,T2,relative %,at which energy]
    lo = float(band.split("-")[0])
    hi = float(band.split("-")[1])
    if int(header[0])==1:
        modele=f"model plasma {T} logt 1 {Nh}"
    elif int(header[0])==2:
        modele=f"model plasma {T} logt 1 {Nh} plasma {header[1]} logt 1 {Nh} {header[2]} {header[3]}"
    else:
        print("ERROR, va va pas marcher car header ni 0 ni 1")
    # Gen d'input
    ####################################################################################################
    f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/Intermediate_file_spectrum/spectrum.txt", "w")
    f.write(f"from xmm pn {filterr} {band}\ninstrument flux ergs {band} unabsorbed")
    f.write(f"\n{modele}")
    f.write(f"\noutput spectrum.mdl {lo} {hi} {step}")
    f.close()
    ####################################################################################################

    # PIMMS Running
    ####################################################################################################
    os.system("/home/thomas/Documents/0_Stage/PIMMS/pimms < /home/thomas/Documents/0_Stage/My_pimms_prgm/Intermediate_file_spectrum/spectrum.txt > /home/thomas/Documents/0_Stage/My_pimms_prgm/Intermediate_file_spectrum/output_spectrum.txt")
    ####################################################################################################

    # Output reading
    ####################################################################################################
    f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/Spectrum_Script/spectrum.mdl", "r")
    liste = f.readlines()
    f.close()
    x = []
    y = []
    for ligne in liste:
        temp = ligne.split()
        x.append(float(temp[0]))
        y.append(float(temp[1]))
    return x,y
'''
T = 7
Nh = 1e19
filterr = "Medium"
band = "0.2-2"
step = 0.01
x,y=get_spectrum(T, Nh, filterr, band, step)
plt.plot(x,y)
plt.show()
'''
