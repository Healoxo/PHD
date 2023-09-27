import os
import numpy as np
import matplotlib.pyplot as plt 

# Params
####################################################################################################
Nh_list = [0, 1e18, 1e19, 1e20, 1e21, 1e22]
T = 7 #np.arange(6,8+0.5,0.5) #va de 6 a 8 log(T)
filterr = "Medium"
band = "0.2-2"
lo = float(band.split("-")[0])
hi = float(band.split("-")[1])
step = 0.01
####################################################################################################

# PIMMS Running and Output Reading for each Nh
####################################################################################################
plt.figure(figsize=(8, 6))
for Nh in Nh_list:
    # Gen d'input
    ####################################################################################################
    f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/Intermediate_file_spectrum/spectrum.txt", "w")
    f.write(f"from xmm pn {filterr} {band}\ninstrument flux ergs {band} unabsorbed")
    f.write(f"\nmodel plasma {T} logt 1 {Nh}")
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
    plt.plot(x, y, label=f"Nh={Nh:.0e}")
    ####################################################################################################

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Energy (keV)")
plt.ylabel("Flux")
plt.legend()
plt.show()
####################################################################################################

