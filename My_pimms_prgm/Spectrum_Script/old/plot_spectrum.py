import sys
import numpy as np
import matplotlib.pyplot as plt
'''
band=str(sys.argv[1])
lo,hi=band.split('-')
lo=float(lo)
hi=float(hi)
step=float(sys.argv[2])
Nh= list(sys.argv[3])

f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all_spectrum.txt","r")
liste=f.readlines()
f.close()
x=np.arange(lo,hi+step,step)
y=[]
i=0
for ligne in liste:
    plt.plot(x,ligne.split(),label=f"Nh={Nh[i]}")
    i+=1
plt.yscale('log')
#plt.xscale('log')
plt.title("Echelle log pour y")
plt.xlabel("Energy (keV)")
plt.ylabel("Flux?")
plt.show()
'''

f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/Spectum_Script/spectrum.mdl","r")
liste=f.readlines()
f.close()
x=[]
y=[]
for ligne in liste:
    X=float(ligne.split()[0])
    Y=float(ligne.split()[1])
    x.append(X)
    y.append(Y)
print(x[4],y[4])
plt.plot(x,y)
plt.yscale('log')
plt.xscale('log')
plt.show()
