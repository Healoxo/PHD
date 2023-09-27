import re
import numpy as np
import matplotlib.pyplot as plt
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
