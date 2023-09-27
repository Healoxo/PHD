import re
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def run_pimms(Nh, T, filterr, band, cts,flux_type,header):
    f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt","w")
    f.write(f"from xmm pn {filterr} {band}\ninstrument flux ergs {band} {flux_type}") #En tete XMM de la camera pn
    #f.write(f"from rosat pspc open {band}\ninstrument flux ergs {band} unabsorbed") #En tete du check de ROSAT
    Temp=np.array([round(x,2) for x in np.arange(5.6,8.5,0.05)]) #Je couvre toute la grille existante
    if header[0]==1:
        for T in Temp:
            f.write(f"\nmodel plasma {round(T,2)} logt 1 {Nh}\ngo {cts}")
        f.close()
    elif header[0]==2:
        for T in Temp:
            f.write(f"\nmodel plasma {round(T,2)} logt 1 {Nh} plasma {header[1]} logt 1 {Nh} {header[2]} {header[3]}\ngo {cts}")
        f.close()
    else:
        print("Mauvaise composante de comp")
    print("Input created")

    os.system("/home/thomas/Documents/0_Stage/PIMMS/pimms < /home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt > /home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt ")
    f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt", "r")
    liste=f.readlines()
    f.close()
    liste_flux=[]
    for i in range(len(liste)):
        if re.search(".*PIMMS predicts a.*",str(liste[i])) != None:
            liste_flux.append((liste[i].split()[-2]))

    print("Ouput readable")

    f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt", "a")
    f.write(' '.join(liste_flux) + '\n')
    f.close()

def ECF(Nh,T1:list,T2:list,filterr, band, cts,flux_type,percent,energy):
    f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt","w")
    f.write(f"from xmm pn {filterr} {band}\ninstrument flux ergs {band} {flux_type}") #En tete XMM de la camera pn
    #f.write(f"from rosat pspc open {band}\ninstrument flux ergs {band} unabsorbed") #En tete du check de ROSAT
    if len(T1)==1 and len(T2)==1:
        f.write(f"\nmodel plasma {round(T1[0],2)} logt 1 {Nh} plasma {round(T2[0],2)} logt 1 {Nh} {percent} {energy}\ngo {cts}")
        f.close()
    elif len(list(T1))!=1 and len (list(T2))==1:
        print("T1 est une liste")
        for T in T1:
            f.write(f"\nmodel plasma {round(T,2)} logt 1 {Nh} plasma {round(T2[0],2)} logt 1 {Nh} {percent} {energy}\ngo {cts}")
        f.close()
    elif len(list(T1))==1 and len (list(T2))!=1:
        print("T2 est une liste")
        for T in T2:
            f.write(f"\nmodel plasma {round(T1[0],2)} logt 1 {Nh} plasma {round(T,2)} logt 1 {Nh} {percent} {energy}\ngo {cts}")
        f.close()
    else:
        print("Mauvaise composante de comp")
    print("Input created")

    os.system("/home/thomas/Documents/0_Stage/PIMMS/pimms < /home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt > /home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt ")
    f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt", "r")
    liste=f.readlines()
    f.close()
    liste_flux=[]
    for i in range(len(liste)):
        if re.search(".*PIMMS predicts a.*",str(liste[i])) != None:
            liste_flux.append((liste[i].split()[-2]))

    print("Ouput readable")

    f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt", "a")
    f.write(' '.join(liste_flux) + '\n')
    f.close()
