f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/Spectum_Script/spectrum.mdl","r")
liste=f.readlines()
f.close()

f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/all_spectrum.txt","a")
liste_flux=[]
for ligne in liste:
    liste_flux.append(ligne.split()[1])
f.write(' '.join(liste_flux) + '\n')
