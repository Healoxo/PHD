import os
import sys
import numpy as np
from definition import *
import matplotlib.pyplot as plt 

try :
    name=str(sys.argv[1])
except :
    name=input("Forgot to put argument Soft | Hard | all\n> ")

try :
    comp=int(sys.argv[2])
except :
    comp=int(input("Forgot to put specify nbr of composant: 1 or 2> "))

if comp==1:
    header=[comp]
elif comp==2:
    try :
        T2=float(sys.argv[3])
    except :
        T2=float(input("Forgot to put specify 2nd temperature in log(T)> "))
    
    try :
        percent=float(sys.argv[4])
    except :
        percent=float(input("Forgot to put specify flux relative to comp. 1\n> "))
    
    try:
        energy=float(sys.argv[5])
    except:
        energy=float(input("Forgot to put specify at which energy comp 1 and 2 are compared (in keV)\n> "))
    
    header=[comp,T2,percent,energy]
else:
    print("Mauvais nombre de comp, fin du prgm")
#print("Le header ressemble a ",header)
if name=="Soft": band="0.2-2"
if name=="Hard": band="2-12"
if name=="all": band="0.2-12"
filterr = "Medium"
step=0.01
Nh_list = [1e0, 1e18, 1e19, 1e20, 1e21, 1e22] 
T_list = np.arange(6, 8.51, 0.5)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))

for i, Nh in enumerate(Nh_list):
    row = i // 3
    col = i % 3
    ax = axs[row, col]
    for T in T_list:
        x,y=get_spectrum(T, Nh, filterr, band, step,header)
        if T==T_list[0]:
            ax.plot(x, y, label=f"T={10**T:.0e}", zorder=0)
        elif T==T_list[-1]:
            ax.plot(x, y, label=f"T={10**T:.0e}", zorder=0)
        else:
            ax.plot(x, y, label=f"T={10**T:.0e}", zorder=1,linestyle='dashed',linewidth=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Total model flux")
    ax.set_title(f"log(Nh)={round(np.log10(Nh))} | step={step}")
    ax.legend()

plt.tight_layout()
if comp==1:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/{name}_Spectrum_Nh_fixed_1_comp.pdf",format='pdf')
if comp==2:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/{name}_Spectrum_Nh_fixed_1_comp_T2_{T2}.pdf",format='pdf')
#plt.show()
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
for i, T in enumerate(T_list):
    row = i // 3
    col = i % 3
    ax = axs[row, col] # get the current subplot
    for Nh in Nh_list:
        x, y = get_spectrum(T, Nh, filterr, band, step,header)
        if Nh == Nh_list[0]:
            ax.plot(x, y, label=f"Nh={Nh:.0e}", zorder=0)
        elif Nh == Nh_list[-1]:
            ax.plot(x, y, label=f"Nh={Nh:.0e}", zorder=0)
        else:
            ax.plot(x, y, label=f"Nh={Nh:.0e}", zorder=1, linestyle='dashed', linewidth=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Total model flux")
    ax.set_title(f"log(T)={T}. step={step}")
    ax.legend()

plt.tight_layout() # adjust spacing between subplots
if comp==1:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/{name}_Spectrum_T_fixed.pdf", format='pdf')
if comp==2:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/{name}_Spectrum_T_fixed_1_comp_T2_{T2}.pdf",format='pdf')
#plt.show()
