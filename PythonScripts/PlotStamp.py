#!/bin/env python

import sys
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.constants import N_A, calorie
import seaborn as sns


file = "Stamp.dat"
ensamble = "NPT" # sys.argv[1]

# Figures styles
sns.set(rc={'figure.figsize':(10,8), 'axes.labelsize': 16})
sns.set_theme("talk")
sns.set_style("white")

print(file, ensamble)

# Names of columns:

# NVT, NVE
###########

#Colonne     1 : Iteration 
#Colonne     2 : Energie totale par particule (en J)
#Colonne     3 : Energie potentielle par particule (en J)
#Colonne     4 : Energie potentielle intramoleculaire par particule (en J)
#Colonne     5 : Energie potentielle intermoleculaire par particule (en J)
#Colonne     6 : Energie cinetique par particule (en J)
#Colonne     7 : Temperature totale (en K)
#Colonne     8 : Composante x de la temperature (en K)
#Colonne     9 : Composante y de la temperature (en K)
#Colonne    10 : Composante z de la temperature (en K)
#Colonne    11 : Pression totale (en Pa)
#Colonne    12 : Composante x de la pression (en Pa)
#Colonne    13 : Composante y de la pression (en Pa)
#Colonne    14 : Composante z de la pression (en Pa)
#Colonne    15 : Composante x de la vitesse du centre de gravite (en m/s)
#Colonne    16 : Composante y de la vitesse du centre de gravite (en m/s)
#Colonne    17 : Composante z de la vitesse du centre de gravite (en m/s)
#Colonne    18 : Derive (en %) 
#Colonne    19 : Temps CPU par iteration (en s)


# NPT
##########

#Colonne     1 : Iteration 
#Colonne     2 : Energie totale par particule (en J)
#Colonne     3 : Energie potentielle par particule (en J)
#Colonne     4 : Energie potentielle intramoleculaire par particule (en J)
#Colonne     5 : Energie potentielle intermoleculaire par particule (en J)
#Colonne     6 : Energie cinetique par particule (en J)
#Colonne     7 : Temperature totale (en K)
#Colonne     8 : Composante x de la temperature (en K)
#Colonne     9 : Composante y de la temperature (en K)
#Colonne    10 : Composante z de la temperature (en K)
#Colonne    11 : Pression totale (en Pa)
#Colonne    12 : Composante x de la pression (en Pa)
#Colonne    13 : Composante y de la pression (en Pa)
#Colonne    14 : Composante z de la pression (en Pa)
#Colonne    15 : Composante x de la vitesse du centre de gravite (en m/s)
#Colonne    16 : Composante y de la vitesse du centre de gravite (en m/s)
#Colonne    17 : Composante z de la vitesse du centre de gravite (en m/s)
#Colonne    18 : Longueur de maille A (en m)
#Colonne    19 : Longueur de maille B (en m)
#Colonne    20 : Longueur de maille C (en m)
#Colonne    21 : d(A)/dt[0] (en m/s)
#Colonne    22 : d(B)/dt[1] (en m/s)
#Colonne    23 : d(C)/dt[2] (en m/s)
#Colonne    24 : Densite (en kg/m3)
#Colonne    25 : Derive (en %) 
#Colonne    26 : Temps CPU par iteration (en s)


names = {"NVT" : [
    "I",
    "Etot", "Epot", "Epot_intra", "Epot_inter", "Ekin",
    "T", "Tx", "Ty", "Tz",
    "P", "Px", "Py", "Pz",
    "Vx", "Vy", "Vz",
    "D", "cpu"],
         "NPT" : [
    "I",
    "Etot", "Epot", "Epot_intra", "Epot_inter", "Ekin",
    "T", "Tx", "Ty", "Tz",
    "P", "Px", "Py", "Pz",
    "Cmvx", "Cmvy", "Cmvz",
    "Lx", "Ly", "Lz",
    "Vx", "Vy", "Vz",
    "Dens",
    "D", "cpu"
    
         ]
}

def load_data(file, t="NVT"):
    data = pd.read_csv(file,
                       sep="\s+",
                       header=None,
                       names=names[t],
                       comment="#"
                      )
    data["Etot"] = data["Etot"] * N_A / 1000  # to kJ/mol
    return data


# Load data
tab = load_data(file, t=ensamble)

print(tab)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

axs[0].plot(tab["I"] / 1e-12, tab["T"], color="red", alpha=0.6)
axs[0].set_ylabel("Temperature (K)")
axs[0].set_xlabel("time (ps)")

axs[1].plot(tab["I"] / 1e-12, tab["Dens"], color="olive", alpha=0.6)
axs[1].set_ylabel("Density (kg/m$^3$)")
axs[1].set_xlabel("time (ps)")

fig.tight_layout()
plt.show()