import pyvo as vo
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import re
import astropy
from astropy import table
from astropy.io import fits,ascii
from astropy.table import QTable,Table,vstack
from astropy.coordinates import SkyCoord,ICRS, Galactic
import astropy.units as u
from astropy.units import Quantity
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch
from scipy.stats import norm


size = Quantity(19, unit="arcmin")
big_radius=Quantity(3, unit="deg")
nom=["Gaia","Xmm"]

url_gaia="https://gea.esac.esa.int/tap-server/conesearch?TABLE=gaiadr3.gaia_source&IDCOL=source_id&RACOL=ra&DECCOL=dec&"
url_gaia2="http://vizier.cds.unistra.fr/viz-bin/conesearch/I/355/gaiadr3?"
url_xmm="http://vizier.cds.unistra.fr/viz-bin/conesearch/IX/65/xmm4d11s?"
url_obs_xmm="http://vizier.cds.unistra.fr/viz-bin/conesearch/IX/65/summary?"

dic={"Gaia":[url_gaia,"designation","ra","dec"],"Xmm":[url_xmm,"Source","RA_ICRS","DE_ICRS"],"Gaia2":[url_gaia2,"DR3Name","RA_ICRS","DE_ICRS"]}

def Sturge(liste):
    return int(1+3.322*np.log(len(liste)))


#########################################################################################################################################################
def conesearch(url_service , position:list , radius , units:str): 
    #### Transformation des positions et rayons ####
    coords=SkyCoord(position[0],position[1],unit="deg",frame="icrs")
    R=Quantity(radius, unit=units)
    
    #### Call du services ####
    service = vo.dal.SCSService(url_service)
    
    #### Conesearch in the service ####
    results = service.search(coords,radius=R)
    table = results.to_table() #Tranforme l'output en table Astropy

    return table
    
def Subcatalog(name:str , position:list , radius , units1:str ,sub_radius , units2): #Prend en entree le nom (Gaia ou Xmm), les coordonee et le rayon du gros conesearch
    liste_obs=conesearch(url_obs_xmm,[position[0] , position[1]],radius,units1)
    for i in range(len(liste_obs)):
        if i==0:
            subcatalog = conesearch(dic[name][0],[liste_obs["RAJ2000"][i] , liste_obs["DEJ2000"][i]],sub_radius,units2)
            if name=="Xmm" or name=="Gaia2":
                subcatalog.remove_column('_r')
        else:
            temp = conesearch(dic[name][0],[liste_obs["RAJ2000"][i] , liste_obs["DEJ2000"][i]],sub_radius,units2)
            if name=="Xmm" or name=="Gaia2":
                temp.remove_column('_r')
            subcatalog = vstack([subcatalog,temp])
            subcatalog = table.unique(subcatalog,keys=dic[name][1]) 
            print(f'\r{round(i/(len(liste_obs)-1) * 100)} % pour {name}',end=' ')
    c=SkyCoord(ra=subcatalog[dic[name][2]],dec=subcatalog[dic[name][3]],frame='icrs')
    glon=c.galactic.l
    glat=c.galactic.b
    subcatalog.add_columns([glon,glat],names=["glon","glat"])
    return subcatalog
    
def Crossing(names:list , position:list , radius , units1:str , sub_radius , units2): #names:nom des 2 catalogues
    colors=['orange','blue']
    marks=['+',"x"]
    final=[]
    a=float(input("show results?\n0:No\n1:Yes\n"))
    for i in range(len(names)):
        final.append(Subcatalog(names[i],position,radius,units1, sub_radius , units2))
        if a==1:
            plt.scatter(final[i]["glon"],final[i]["glat"],color=colors[i],label=names[i],marker=marks[i],zorder=i)
    if a==1:
        plt.legend()
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.title("Orion region")
        plt.gca().invert_xaxis()
        plt.show()
    return final[0],final[1] #Renvois les 2 sous catalogues
#########################################################################################################################################################


#def PMmodulus(catalog):
    #return np.sqrt(catalog['pmra']**2+catalog['pmdec']**2)


def crossmatch_gaia(xmm_cat,selec,distance_max=5): #selec= all ou best
    return XMatch.query(cat1=xmm_cat,
                     cat2='vizier:I/355/gaiadr3',
                     max_distance=distance_max * u.arcsec, colRA1='RA_ICRS',
                     colDec1='DE_ICRS',selection=selec)

def object_to_str(table):
    for colname in table.colnames:
        if table[colname].dtype=='object':
            table[colname]=table[colname].astype('str')

def gaussian(x, mu, sig):
    return 1./(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

def Gaussian_fit(xmin,xmax,data):
    temp=data[(data>xmin) & (data<xmax)]
    mean,std=norm.fit(temp)
    x=np.linspace(xmin,xmax,1000)
    p=norm.pdf(x,mean,std)
    return x,p,mean,std

def Median_scatter(x,y,space,plot=0):
    mini=int(np.min(x))
    maxi=int(np.max(x))+1
    interval=int((maxi-mini))
    median=[]
    stdd=[]
    for i in range(int(interval*(1/space))):
        g=mini + i*space
        d=mini + (i+1)*space
        #plt.axvline(mini + i*step)
        mean,std=norm.fit(y[(x>g) & (x<d)])
        stdd.append(std)
        median.append(np.median(np.sort(y[(x>g) & (x<d)])))
    if plot==1:
        X=np.linspace(mini+space/2,maxi-space/2,len(median))
        plt.scatter(x,y,zorder=0)
        plt.plot(X,median,color='black',zorder=1)
        plt.scatter(X,median,color='black',zorder=1)
        plt.gca().invert_xaxis()
        return median
    else:
        return np.array(median),np.array(stdd),mini+space/2,maxi-space/2

#########################################################################################################################################################
#Pour lire les isochrone que me renvois le site 
def readfile(path):       #Renvois la table ASCII en table Astropy
    f=open(f"{path}")
    liste=f.readlines()
    f.close()
    a=0
    while a==0:
        ligne=liste[1]
        if ligne.split()[0]=="#":
            liste.pop(0)
        else:
            a=1
    cols=liste[0].split()
    cols.pop(0)
    hehe=table.Table.read(path,format='ascii',delimiter=' ',comment='#')
    hehe=Table(hehe,names=cols)
    return hehe


def XY(table,names:list): #La list est de la forme [Gmag,Bp,Rp]
    x=table[names[1]]-table[names[2]]
    y=table[names[0]]
    return x,y
#########################################################################################################################################################
def Comparaison(data1,data2,BIN): #Input 2 list, return the value of their histogram, their ration (2/1), the errorbar of the ratio
    r=[np.nanmin(list(data1)+list(data2)),np.nanmax(list(data1)+list(data2))]
    hist1=np.histogram(data1,bins=BIN,range=r)
    hist2=np.histogram(data2,bins=BIN,range=r)
    errors=hist2[0]/hist1[0]*np.sqrt(1/hist2[0] + 1/hist1[0])
    ratio=hist2[0]/hist1[0]
    poserrors=[]
    for i in range(len(hist1[1])-1):
        poserrors.append((hist1[1][i]+hist1[1][i+1])/2)
    return hist1,hist2,ratio,errors,poserrors



def Plot_Comparaison(data1:list,data2:list,BIN,labelX,save=0,name=""): #data contient la liste et le label [donnee,label]
    hist1,hist2,ratio,errors,poserrors=Comparaison(data1[0],data2[0],BIN)
    plt.figure(figsize=(10,5))
    ax1=plt.subplot(2,1,1)
    ax1.stairs(hist1[0],hist1[1],label=data1[1],color='orange')
    ax1.stairs(hist2[0],hist2[1],label=data2[1],color='blue')
    ax1.set_ylabel("Numb. of stars", fontsize=16)
    ax1.legend()

    ax2=plt.subplot(2,1,2,sharex=ax1)
    ax2.stairs(ratio,hist1[1],color='red')
    ax2.errorbar(poserrors,ratio,yerr=errors,linestyle="",c="black")
    ax2.set_ylim(0,1.1)
    ax2.set_xlabel(labelX, fontsize=16)
    ax2.set_ylabel("Ratio",fontsize=16)
    plt.setp(ax1.get_xticklabels(), visible=False)
    if save==1:
        while name=="":
            name=str(input("Select name of the file:\n"))
        plt.savefig(f"/home/thomas/Documents/0_Stage/Code/Images/{name}.pdf")
#########################################################################################################################################################
def MeanDist(curveXY,dataX,dataY,BIN):
    bins=np.linspace(np.min(curveXY[0]),np.max(curveXY[0]),BIN)
    midbins=[(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)]
    
    meandist=[]
    sigma=[]
    for i in range(len(bins)-1):
        mask_bins1=(curveXY[0]>bins[i]) & (curveXY[0]<bins[i+1])
        mask_bins2=(dataX>bins[i]) & (dataX<bins[i+1])
        
        meany=np.mean(curveXY[1][mask_bins1])
        meandist.append(np.mean(np.abs(dataY[mask_bins2]-meany)))
        sigma.append(np.sqrt(len(dataY[mask_bins2])))
    return meandist,midbins,sigma

parsec2_gaia_dr2=["G_fSBmag","G_BP_fSBmag","G_RP_fSBmag"]
parsec1_2s_gaia_dr2=["Gmag","G_BPmag","G_RPmag"]

def Iso(filename,age,metal): #Lis le fichier contenant tout les isochrones, et renvois celle avec l'age/metalicite demande
    TABLE=readfile(filename)
    mask_age=(TABLE["logAge"]==age)
    mask_metal=(TABLE["Zini"]==metal)
    mastermask=mask_age & mask_metal
    x,y=XY(TABLE,parsec1_2s_gaia_dr2) #Pour l'instant je travail qu'avec parsec 1.2, mais le nom des colones est a generaliser
    x=x[mastermask]
    y=y[mastermask]
    return x,y




























