import numpy as np
import matplotlib.pyplot as plt
import sys
Nh=[0,18,19,20,21,22,30]
count=float(sys.argv[1])
f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt","r")
liste=f.readlines()
x=np.arange(5.6,8.5,0.05)
for i in range(len(liste)):
    temp=count/np.array([float(y) for y in liste[i].split()])
    plt.plot(np.linspace(5.6,8.5,len(temp)),temp/(10**11),label=f"Nh={Nh[i]}")
#plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0)
plt.legend()
plt.title(f"On a fait ECF=cts/Flux. Ici cts={count}")
#plt.yscale('log')
plt.ylim(0,10)
plt.xlabel("Temperature (log(T))")
plt.ylabel(r"ECF $(10^{11} \; erg^{-1} cm^{2})$")
#plt.savefig("/home/thomas/Documents/0_Stage/My_pimms_prgm/Results/Flux_vs_T_XMM_Soft.pdf",format='pdf', bbox_inches='tight')
plt.show()
