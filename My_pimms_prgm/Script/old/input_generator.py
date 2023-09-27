import re
import numpy as np
import matplotlib.pyplot as plt
import sys
Nh=float(sys.argv[1])
cts=float(sys.argv[2])
filterr=str(sys.argv[3])
band=str(sys.argv[4])


f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt","w")
f.write(f"from xmm pn {filterr}  {band}\ninstrument flux ergs {band} unabsorbed") #En tete XMM de la camera pn
#f.write(f"from rosat pspc open {band}\ninstrument flux ergs {band} unabsorbed") #En tete du check de ROSAT
Temp=np.arange(5.6,8.5,0.05) #Je couvre toute la grille existante
for T in Temp:
    f.write(f"\nmodel plasma {T} logt 1  {Nh}\ngo {cts}")
f.close()

print("Input created")

