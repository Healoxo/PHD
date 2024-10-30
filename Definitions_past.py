############################################################################################################################################################################

import numpy as np
import warnings
from astropy.table import Table, join, Column, vstack
from astropy.coordinates import SkyCoord,ICRS, Galactic
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.spatial import cKDTree
from astropy.time import Time
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter,convolve
from scipy.stats import linregress
import healpy as hp
from collections import Counter

# %matplotlib widget

warnings.filterwarnings('ignore')
############################################################################################################################################################################










































############################################################################################################################################################################
#1
def kick_n_print(condition , str_condition, col_name , table):
    table[f"{col_name}"] = False
    table[f"{col_name}"][condition] = True

    table[f"not_{col_name}"] = False
    table[f"not_{col_name}"][~condition] = True

    table["Keep"][condition] = False
    print(f"Proportion perdu si on retire {str_condition}:",f"{round(len(table[condition])/len(table) * 100)}%")

############################################################################################################################################################################


































############################################################################################################################################################################
#2
def lb_to_xyz(Glon, Glat, d):
    x_sun = 8330  # Distance from the galactic center to the Sun in kpc
    y_sun = 0     # y-coordinate of the Sun, assumed to be 0 kpc
    z_sun = 25 # z-coordinate of the Sun in kpc

    x = d * np.cos(np.radians(Glon)) * np.cos(np.radians(Glat))
    y = d * np.sin(np.radians(Glon)) * np.cos(np.radians(Glat))
    z = d * np.sin(np.radians(Glat))

    x -= x_sun
    y -= y_sun
    z -= z_sun
    return x, y, z

def add_or_replace_column(table, column_data, column_name):
    try:
        table.add_column(column_data, name=column_name)
    except:  # Column already exists
        table.replace_column(column_name, column_data)

def haversine_arcsec(ra1, dec1, ra2, dec2):
    # Convertir les angles en radians
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2]) #map applique fonction à éléments de la listeµ
    
    # Différences de coordonnées
    delta_ra = ra2 - ra1
    delta_dec = dec2 - dec1
    
    # Formule de Haversine
    a = np.sin(delta_dec / 2.0)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(delta_ra / 2.0)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    # Convertir en arcsecondes
    return np.degrees(c) * 3600
############################################################################################################################################################################











































############################################################################################################################################################################
#3


def pos_to_pixel(order, glon, glat):
    nside = 2**order
    pixels = hp.ang2pix(nside, glon, glat, lonlat=True)
    return pixels

def calculate_density_map(order, glon, glat):
    nside = 2**order
    """Calculates the density map and returns the pixel-density pairs for non-zero densities."""
    # Convert angular coordinates to pixel indices
    pixels = hp.ang2pix(nside, glon, glat, lonlat=True)
    
    density_map = np.bincount(pixels, minlength=hp.nside2npix(nside))
    
    non_zero_pixels = np.nonzero(density_map)[0]
    non_zero_densities = density_map[non_zero_pixels]
    
    # Return pairs of (pixel, density) where density > 0
    return list(zip(non_zero_pixels, non_zero_densities))

def line(x , a , b):
    return a * x + b

def Rayleigh(x, c , s):
    s=1
    return c * (x / s**2) * np.exp(-x**2 / (2 * s**2))

def pT(x):
    return x*np.exp(-x**2 / 2)

def pF(x , Density):
    return x**2 * np.pi * Density 

def LR(x , Density):
    return 1/(2*np.pi * Density) * np.exp(-x**2 / 2)
    # return pT(x) / pF(x , Density)

def Rayleigh_plus_line(x, a , b , c , s):
    return Rayleigh(x, c, s) + line(x, a, b)




def MSE(distribution1, distribution2):
    if len(distribution1) != len(distribution2):
        raise ValueError("Les distributions doivent avoir la même longueur")
    
    mse = np.mean((np.array(distribution1) - np.array(distribution2)) ** 2)
    return mse


# def compute_a_b(data, bins):
#     midbins = 0.5 *(bins[1:] + bins[:-1])
#     mask_bins =  (midbins < 3.3)
#     a,b = np.polyfit(np.insert(midbins[mask_bins],0,0) ,np.insert(np.histogram( data , bins = bins)[0][mask_bins],0,0) 
#                      , deg = 1 , w = np.insert(np.ones(len(midbins[mask_bins])),0,1e3))
#     return a,b

# def compute_sigma_and_norm(data, bins, a, b):
#     midbins = np.insert(0.5 * (bins[1:] + bins[:-1]), 0, 0)
#     histo = np.insert(np.histogram(data, bins=bins)[0], 0, 0)
#     fr = histo - line(midbins , a , b) #Here I take my data and substract the line of spurious associations
#     fr[fr < 0] = 0
#     params, _ = curve_fit(Rayleigh, midbins, fr, p0=[np.max(fr), np.std(data)])
#     c, s = params
#     r2 = r_squared(fr, Rayleigh(midbins, *params))
#     if (c>10000) & (s>10): #Constante trop grande ( to be replaced by np.max(histo)???)
#         # print("c and s too big:\nc=",c,"\ns=",s)
#         c = np.nan
#         s = np.nan
#         r2 = np.nan
#     return c , s , r2


def compute_a_b(data, bins):
    midbins = 0.5 *(bins[1:] + bins[:-1])
    mask_bins =  (midbins < 3.3)
    a,b = np.polyfit(np.insert(midbins[mask_bins],0,0) ,np.insert(np.histogram( data , bins = bins , density = True)[0][mask_bins],0,0) 
                     , deg = 1 , w = np.insert(np.ones(len(midbins[mask_bins])),0,1e3))
    return a,b

def compute_sigma_and_norm(data, bins, a, b):
    midbins = np.insert(0.5 * (bins[1:] + bins[:-1]), 0, 0)
    histo = np.insert(np.histogram(data, bins=bins,density = True)[0], 0, 0)
    histemp = np.insert(np.histogram(data, bins=bins,density = False)[0], 0, 0)
    fr = histo - line(midbins , a , b) #Here I take my data and substract the line of spurious associations
    fr[fr < 0] = 0
    erreurs = np.sqrt(histemp) * (np.max(histo)/np.max(histemp))
    params, _ = curve_fit(Rayleigh, midbins, fr, p0=[np.max(fr), np.std(data)], maxfev = 6000)# , sigma = erreurs)
    c, s = params
    mse = 1 - MSE(fr, Rayleigh(midbins, c, s))
    if (c>10000) & (s>10): #Constante trop grande ( to be replaced by np.max(histo)???)
        # print("c and s too big:\nc=",c,"\ns=",s)
        c = np.nan
        s = np.nan
        mse = np.nan
    return c , s , mse, fr

def CR(x , a , b , c , s): #Inégrale analytique
    return (a/2 * x**2 + b*x) / (a/2 * x**2 + b*x + c * (1 - np.exp(-(x**2/(2*s**2)))) )

def all_in_one(data, data_shifted, bins):
    a,b = compute_a_b(data_shifted, bins)
    c,s,mse, fr = compute_sigma_and_norm(data, bins,a,b)
    return a , b , c , s , 1 - CR(data,a,b,c,s) , mse

 


def data_to_hist2D(x,y,minmax_x,minmax_y,N):
    x_bins = np.linspace(minmax_x[0], minmax_x[1], N)
    y_bins = np.linspace(minmax_y[0], minmax_y[1], N)
    counts = np.histogram2d(x, y, bins=(x_bins, y_bins))[0].T
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    y_centers = 0.5 * (y_bins[:-1] + y_bins[1:])
    x_mesh, y_mesh = np.meshgrid(x_centers, y_centers)
    return x_mesh, y_mesh, counts



def find_nearest(grille, x, y,name:str):
    # Conversion des colonnes en numpy array
    grille_x = np.array(grille['x'])
    grille_y = np.array(grille['y'])
    
    # Création de l'arbre KD
    tree = cKDTree(np.vstack((grille_x, grille_y)).T)
    
    # Recherche des points les plus proches
    _ , indices = tree.query(np.c_[x, y])
    
    # Récupération des probabilités associées aux indices trouvés
    proba = np.array(grille[f"{name}"])[indices]
    
    return proba
############################################################################################################################################################################

#For Proba angDist

def initialize_table(table, order):
    try:
        table.add_columns([Column(np.zeros(len(table)), dtype='i4'), Column(np.zeros(len(table)), dtype='i4')], names=["order", "HealPIX"])
    except:
        pass
    pixels = pos_to_pixel(order, table["glon"], table["glat"])
    table["order"] = [order] * len(pixels)
    table["HealPIX"] = pixels
    return table

def create_healpix_table(order, table):
    pairs = calculate_density_map(order, table["glon"], table["glat"])
    pixels, No = zip(*pairs) ; pixels = list(pixels) 
    No = list(No)
    Nx = np.array([len(np.unique(table["IAUName"][table["HealPIX"] == pixel])) for pixel in pixels])
    Density = No / hp.nside2pixarea(order ** 2, degrees=True)
    healpix_table = Table(
        [[order] * len(pixels), pixels, np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), No, Nx , Density ], #, np.zeros(len(pixels),dtype = bool)],
        names=('order'        , 'HealPIX', 'a'               , 'b'                  , 'c'                  , 's'                  , 'r2'                 ,'No','Nx', 'Density'), #, "Outliers"),
        dtype=('i4'   , 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8') #, 'bool')
    )
    # healpix_table["Outliers"][healpix_table["Nx"] < 100] = True
    return healpix_table

    # return Table(
        # [[order] * len(pixels), pixels, np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)), No, Nx , Density],
        # names=('order'        , 'HealPIX', 'a'               , 'b'                  , 'c'                  , 's'                  , 'r2'                 ,'No','Nx', 'Density'),
        # dtype=('i4'   , 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')
    # )

def update_healpix_table(order, Healpix_table, table, table_shifted, threshold_x, threshold_ratio):
    while True:
        mask = (
            (Healpix_table["Nx"] > threshold_x) &
            (Healpix_table["No"] / Healpix_table["Nx"] > threshold_ratio) &
            (Healpix_table["order"] == order)
        )
        pix_to_change = np.array([pixel for pixel in Healpix_table["HealPIX"]])[mask]
        if len(pix_to_change) == 0:
            break
        else:
            pass

        mask_pix = (np.isin(table["HealPIX"], pix_to_change)) & (table["order"] == order)
        glon_to_change = table["glon"][mask_pix]
        glat_to_change = table["glat"][mask_pix]

        mask_shifted_pix = (np.isin(table_shifted["HealPIX"], pix_to_change)) & (table_shifted["order"] == order)
        glon_shifted_to_change = table_shifted["glon"][mask_shifted_pix]
        glat_shifted_to_change = table_shifted["glat"][mask_shifted_pix]

        Healpix_table.remove_rows([i for i, x in enumerate(mask) if x])

        order += 1
        table["order"][mask_pix] = order
        table_shifted["order"][mask_shifted_pix] = order

        pixels = pos_to_pixel(order, glon_to_change, glat_to_change)
        table["HealPIX"][mask_pix] = pixels
        del pixels

        pixels = pos_to_pixel(order, glon_shifted_to_change, glat_shifted_to_change)
        table_shifted["HealPIX"][mask_shifted_pix] = pixels
        del pixels

        pairs = calculate_density_map(order, glon_to_change, glat_to_change)
        pixels, No = zip(*pairs)
        pixels = list(pixels)
        No = list(No)
        temp = table[mask_pix]
        Nx = np.array([len(np.unique(temp["IAUName"][temp["HealPIX"] == pixel])) for pixel in pixels])
        Density = No / hp.nside2pixarea(order ** 2, degrees=True)
        temp = Table(
            [[order] * len(pixels), pixels, np.zeros(len(pixels)), np.zeros(len(pixels)), np.zeros(len(pixels)),
             np.zeros(len(pixels)), np.zeros(len(pixels)), No, Nx, Density],
            names=('order', 'HealPIX', 'a', 'b', 'c', 's', 'r2', 'No', 'Nx', 'Density'),
            dtype=('i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')
        )
        Healpix_table = vstack([Healpix_table, temp])
    return Healpix_table

def initialize_probability_columns(table):
    table["Probability_angDist"] = np.zeros(len(table), dtype=float)
    table["a"] = np.zeros(len(table), dtype=float)
    table["b"] = np.zeros(len(table), dtype=float)
    table["c"] = np.zeros(len(table), dtype=float)
    table["s"] = np.zeros(len(table), dtype=float)
    table["r2"] = np.zeros(len(table), dtype=float)
    table["No"] = np.zeros(len(table), dtype=float)
    table["Nx"] = np.zeros(len(table), dtype=float)
    table["Density"] = np.zeros(len(table), dtype=float)
    return table

def calculate_probabilities(Healpix_table, table, table_shifted, bins):
    # liste_order = list(np.unique(Healpix_table["order"]))
    liste_order = [int(x) for x in np.unique(Healpix_table["order"])]
    print("Liste des ordres:", liste_order)
    i = 0
    for order in liste_order:
        mask_table_order = (table["order"] == order)
        pixels = Healpix_table["HealPIX"][Healpix_table["order"] == order]
        for pixel in pixels:
            Nx = Healpix_table["Nx"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)][0]
            No = Healpix_table["No"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)][0]
            Density = Healpix_table["Density"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)][0]
            mask = (table["HealPIX"] == pixel) & (mask_table_order)
            data = table["angDist"][mask] / table["s_Pos"][mask]
            mask_shifted = (table_shifted["HealPIX"] == float(pixel)) & (table_shifted["order"] == order)
            data_shifted = table_shifted["angDist"][mask_shifted] / table_shifted["s_Pos"][mask_shifted]
            try:
                a, b, c, s, proba, r2 = all_in_one(data, data_shifted, bins)
            except:
                a,b,c,s,proba,r2 = np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
            table["Probability_angDist"][mask] = proba
            Healpix_table["a"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)] = float(a)
            table["a"][mask] = float(a)
            Healpix_table["b"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)] = b
            table["b"][mask] = b
            Healpix_table["c"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)] = c
            table["c"][mask] = c
            Healpix_table["s"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)] = s
            table["s"][mask] = s
            Healpix_table["r2"][(Healpix_table["HealPIX"] == pixel) & (Healpix_table["order"] == order)] = r2
            table["r2"][mask] = r2
            table["No"][mask] = No
            table["Nx"][mask] = Nx
            table["Density"][mask] = Density

            print(f"Complétion: {round(i / len(Healpix_table) * 100, 2)} %", end="\r")
            i += 1
    print("Complétion: 100 %")

def rank_probabilities(table,proba_str:str):
    table[f"Rank_{proba_str}"] = np.zeros(len(table), dtype=int)
    sorted_index = np.lexsort((-table[proba_str], table["IAUName"]))
    table_copy = Table(np.copy(table[sorted_index]))

    _, counts = np.unique(table_copy['IAUName'], return_counts=True)
    ranks = np.concatenate([np.arange(1, count + 1) for count in counts])
    table[f"Rank_{proba_str}"][sorted_index] = ranks
    table[f"Top_1_{proba_str}"] = False
    table[f"Top_1_{proba_str}"][table[f"Rank_{proba_str}"] == 1] = True
    del counts, ranks, sorted_index, table_copy
    return table


'''
table[f"Probability_total"] = 0.
table[f"Probability_total"] = table["Probability_angDist"] * table[f"Probability_color_vs_flux_ratio"]



# Tri en fonction de 'Probability_total' (ordre décroissant) et 'IAUName'
sorted_index = np.lexsort((-table[f"Probability_total"], table["IAUName"]))  # le "-" est là pour simuler un ordre décroissant
table_copy = Table(np.copy(table[sorted_index]))  # Créer une copie triée

# Initialiser la colonne des rangs dans 'table'
table['ProbabilityRank'] = 0

# Utilisation de numpy pour obtenir les tailles de chaque groupe d'IAUName
unique_iau, counts = np.unique(table_copy['IAUName'], return_counts=True)

# Générer les rangs par groupe d'IAUName
ranks = np.concatenate([np.arange(1, count + 1) for count in counts])

# Assigner ces rangs dans la table triée
# table_copy['ProbabilityRank'] = ranks #Inutile ???

# Remplacer les rangs dans la table originale en utilisant les index triés
table['ProbabilityRank'][sorted_index] = ranks

# Flag de la probabilité la plus élevée
table["HighestProbability"] = False
table["HighestProbability"][table['ProbabilityRank'] == 1] = True



del counts , ranks , sorted_index
'''

############################################################################################################################################################################
































############################################################################################################################################################################
#4

def CR_computation(table , proba_str:str , liste_threshold ): #Completion and reliability diagram
    total_sum = np.sum(table[proba_str]) * 1/0.9973 #Je calcul la somme totale, et divise par le completness factor
    sorted_probs = np.sort(table[proba_str])[::-1] # je trie du plus grand au plus petit


    total_sum_unique = np.sum(table[proba_str][(table[f"Top_1_{proba_str}"])]) * 1/0.9973 #Je calcul la somme totale, et divise par le completness factor
    sorted_probs_unique = np.sort(table[proba_str][(table[f"Top_1_{proba_str}"])])[::-1] # je trie du plus grand au plus petit

    completeness = [] ; completeness_unique = []
    reliability = [] ; reliability_unique = []
    neighbours = [] ; at_least_two_neighbours = []

    for threshold_proba in liste_threshold:
        above_threshold = sorted_probs[sorted_probs >= threshold_proba]
        above_threshold_unique = sorted_probs_unique[sorted_probs_unique >= threshold_proba]
        
        if len(above_threshold) > 0:
            completeness_value = np.sum(above_threshold) / total_sum
            reliability_value = np.sum(above_threshold) / len(above_threshold)
            completeness_value_unique = np.sum(above_threshold_unique) / total_sum_unique
            reliability_value_unique = np.sum(above_threshold_unique) / len(above_threshold_unique)
            Neighbours_value = len(above_threshold)/len(above_threshold_unique)
        else:
            completeness_value = np.nan
            reliability_value = np.nan
            completeness_value_unique = np.nan
            reliability_value_unique = np.nan
            Neighbours_value = np.nan

        completeness.append(completeness_value)
        reliability.append(reliability_value)
        completeness_unique.append(completeness_value_unique)
        reliability_unique.append(reliability_value_unique)
        neighbours.append(Neighbours_value)

    return completeness, reliability, completeness_unique, reliability_unique, neighbours

def CR_computation_2(table , proba_str:str , liste_threshold ): #Completion and reliability diagram Pa concluant celui là
    total_sum = np.sum(table[proba_str]) * 1/0.9973 #Je calcul la somme totale, et divise par le completness factor
    # sorted_probs = np.sort(table[proba_str])[::-1] # je trie du plus grand au plus petit
    idx = np.lexsort( (table['IAUName'] , -table[proba_str]))
    sorted_probs = table["IAUName" , proba_str][idx] # je trie du plus grand au plus petit


    total_sum_unique = np.sum(table[proba_str][(table[f"Top_1_{proba_str}"])]) * 1/0.9973 #Je calcul la somme totale, et divise par le completness factor
    sorted_probs_unique = np.sort(table[proba_str][(table[f"Top_1_{proba_str}"])])[::-1] # je trie du plus grand au plus petit

    completeness = [] ; completeness_unique = []
    reliability = [] ; reliability_unique = []
    neighbours = [] ; at_least_two_neighbours = []

    for threshold_proba in liste_threshold:
        above_threshold = sorted_probs[proba_str][sorted_probs[proba_str] >= threshold_proba]
        above_threshold_unique = sorted_probs_unique[sorted_probs_unique >= threshold_proba]
        
        if len(above_threshold) > 0:
            completeness_value = np.sum(above_threshold) / total_sum
            reliability_value = np.sum(above_threshold) / len(above_threshold)
            completeness_value_unique = np.sum(above_threshold_unique) / total_sum_unique
            reliability_value_unique = np.sum(above_threshold_unique) / len(above_threshold_unique)
            Neighbours_value = len(above_threshold)/len(above_threshold_unique)
            iau_counts = Counter(sorted_probs["IAUName"][sorted_probs[proba_str] >= threshold_proba])
            at_least_two_neighbours_value = len([iau for iau, count in iau_counts.items() if count >= 2]) #Compte le nombre d'occurence des nom, garde ceux qui sont supérieur à 2 dans une table temporaire, puis renvois la taille de la table temporaire
        else:
            completeness_value = np.nan
            reliability_value = np.nan
            completeness_value_unique = np.nan
            reliability_value_unique = np.nan
            Neighbours_value = np.nan
            at_least_two_neighbours_value = np.nan

        completeness.append(completeness_value)
        reliability.append(reliability_value)
        completeness_unique.append(completeness_value_unique)
        reliability_unique.append(reliability_value_unique)
        neighbours.append(Neighbours_value)
        at_least_two_neighbours.append(at_least_two_neighbours_value)

    return completeness, reliability, completeness_unique, reliability_unique, neighbours, at_least_two_neighbours








def plot_data(table, liste, liste_chr, choix , folder , folder_cross , Healpix_folder): #0 = elements, 1 = intervals
    P = ["Probability_angDist" , "Probability_color_vs_flux_ratio" , "Probability_total"] # just for copy paste, "P" not used in code
    colors = ["blue" , "orange" , "green"]
    liste_threshold = np.linspace(0,1,101)

    if choix == 0:
        for order in liste:
            mask = (table[liste_chr] == order)
            fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))  # A4 size in inches (landscape format), with 2 rows and 3 columns
            axes = axes.ravel()  # Flatten the array of axes for easier indexing
            fig.suptitle(f'{liste_chr} = {order} | N={len(table[mask])} | N_unique={len(np.unique(table[mask]["IAUName"]))}', fontsize=20)

            for idx, (Proba, color) in enumerate(zip(P, colors)):
                C, R, C2, R2, N = CR_computation(table[mask], Proba, liste_threshold)

                # Plotting C, C2, R, and R2 in the first row
                ax1 = axes[idx]
                ax1.plot(liste_threshold, C, label=Proba[12:], color=color)
                ax1.plot(liste_threshold, C2, label="Unique", color=color, linestyle="dashed")
                ax1.plot(liste_threshold, R, color=color)
                ax1.plot(liste_threshold, R2, color=color, linestyle="dashed")
                ax1.set_xlim(0, 1)
                ax1.set_ylim(0, 1)
                ax1.legend(fontsize=10)
                ax1.set_xlabel("Probability threshold", fontsize=12)
                ax1.set_ylabel("Percentage", fontsize=12)
                ax1.grid(True)

                # Plotting N in the second row
                ax2 = axes[idx + 3]
                ax2.plot(liste_threshold, N, color=color)
                ax2.set_xlim(0, 1)
                ax2.set_xlabel("Probability threshold", fontsize=12)
                ax2.set_ylabel("Number of neighbours", fontsize=12)
                ax2.grid(True)
            plt.tight_layout()
            plt.savefig(f"/home/t.oliveira/Documents/1_Thesis/1_Code/Table/General_method_{folder}/Cross_{folder_cross}/Healpix_{Healpix_folder}/{liste_chr}/C&R_{order}.pdf")
            plt.close()
    elif choix == 1:
        for i in range(len(liste) - 1):
            mask =   (table[liste_chr] > liste[i]) & (table[liste_chr] < liste[i+1]) 
            fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))  # A4 size in inches (landscape format), with 2 rows and 3 columns
            axes = axes.ravel()  # Flatten the array of axes for easier indexing
            fig.suptitle(f'{liste[i]} < {liste_chr} < {liste[i+1]} | N={len(table[mask])} | N_unique={len(np.unique(table[mask]["IAUName"]))}', fontsize=20)

            for idx, (Proba, color) in enumerate(zip(P, colors)):
                C, R, C2, R2, N = CR_computation(table[mask], Proba, liste_threshold)

                # Plotting C, C2, R, and R2 in the first row
                ax1 = axes[idx]
                ax1.plot(liste_threshold, C, label=Proba[12:], color=color)
                ax1.plot(liste_threshold, C2, label="Unique", color=color, linestyle="dashed")
                ax1.plot(liste_threshold, R, color=color)
                ax1.plot(liste_threshold, R2, color=color, linestyle="dashed")
                ax1.set_xlim(0, 1)
                ax1.set_ylim(0, 1)  
                ax1.legend(fontsize=10)
                ax1.set_xlabel("Probability threshold", fontsize=12)
                ax1.set_ylabel("Percentage", fontsize=12)
                ax1.grid(True)

                # Plotting N in the second row
                ax2 = axes[idx + 3]
                ax2.plot(liste_threshold, N, color=color)
                ax2.set_xlim(0, 1)
                ax2.set_xlabel("Probability threshold", fontsize=12)
                ax2.set_ylabel("Number of neighbours", fontsize=12)
                ax2.grid(True)
            plt.tight_layout()
            plt.savefig(f"/home/t.oliveira/Documents/1_Thesis/1_Code/Table/General_method_{folder}/Cross_{folder_cross}/Healpix_{Healpix_folder}/{liste_chr}/C&R_{liste[i]}-{liste[i+1]}.pdf")
            plt.close()
    elif (choix !=0) & (choix !=1):
        choix = input("0 or 1\n0: elements\n1: intervals")
        plot_data(table, liste, liste_chr, choix , folder , folder_cross , Healpix_folder)