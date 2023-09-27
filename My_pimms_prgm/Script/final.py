import numpy as np
import matplotlib.pyplot as plt
from definition import *
os.system("rm /home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt")
# Paramètres d'entrée
try :
    name=str(sys.argv[1])
except :
    name=str(input("Forgot to put argument: Soft | Hard | all\n"))
try :
    flux_type=str(sys.argv[2])
except :
    flux_type=input("Forgot to put argument: blank | unabsorbed\n")


try :
    comp=int(sys.argv[3])
except :
    comp=int(input("Forgot to specify nbr of composant: 1 or 2\n> "))

if comp==1:
    header=[comp]
elif comp==2:
    try :
        T2=float(sys.argv[4])
    except :
        T2=float(input("Forgot to specify 2nd temperature in log(T)\n> "))
    
    try :
        percent=float(sys.argv[5])
    except :
        percent=float(input("Forgot to specify flux relative to comp. 1\n> "))
    
    try:
        energy=float(sys.argv[6])
    except:
        energy=float(input("Forgot to specify at which energy comp 1 and 2 are compared (in keV)\n> "))
    
    header=[comp,T2,percent,energy]










Nh = Nh=[0,1e18,1e19,1e20,1e21,1e22]  # nombre courbes
T = list(np.arange(5.6, 8.51, 0.05))  # température. le 0.01 en plus est pour forcer np.arange() a fonctioner correctement
if name=="Soft": band="0.2-2"
if name=="Hard": band="2-12"
if name=="all": band="0.2-12"
system="XMM"
filterr = 'Medium'  # filtre
#band = '0.2-2'  # bande passante
cts = 1e-2  # nombre de photons par seconde

# Appel de la fonction run_pimms() pour générer les données
for i in range(len(Nh)):
    run_pimms(Nh[i], T, filterr, band, cts,flux_type,header)

# Lecture des données dans le fichier all.txt
f = open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt","r")
liste = f.readlines()
f.close()

# Calcul de l'ECF
count = cts
ecf = []
for i in range(len(liste)):
    temp = count / np.array([float(y) for y in liste[i].split()])
    ecf.append(temp / (10**11))
#Condition speciale pour la bande Hard
if len(T) != len(ecf[0]):
    T.pop(0)
# Tracé du graphique
for i in range(len(Nh)):
    if i==0:
        plt.plot(T, ecf[i], label=f"Nh={Nh[i]}")
    elif i==(len(Nh)-1):
        plt.plot(T, ecf[i], label=f"Nh={Nh[i]}")
    else:
        plt.plot(T, ecf[i], label=f"Nh={Nh[i]}",linestyle='dashed',linewidth=0.5)
    
plt.legend()
plt.title(f"ECF=cts/Flux | cts={count} | band:{band}")
#plt.ylim(0, 10)
plt.xlabel("Temperature (log(T))")
plt.ylabel(r"ECF $(10^{11} \; erg^{-1} cm^{2})$")
if comp==1:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/ECF_vs_T/ecf_vs_T_{system}_{name}.pdf",format='pdf')
if comp==2:
    plt.savefig(f"/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/ECF_vs_T/ecf_vs_T_{system}_{name}_T2_{T2}.pdf",format='pdf')
plt.close()





































