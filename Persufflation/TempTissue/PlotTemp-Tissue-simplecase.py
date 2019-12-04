#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 09:27:57 2019

@author: anujchaudhri
"""
from math import log
from math import exp
import numpy as np
import matplotlib.pyplot as plt

Ttv = 269.15
Tgi = 173.15

rg = 1500e-6
rt = 15000e-6
L = 50e-3

# Argon
rhog_Ar = 2.7875 #kg/m^3
cg_Ar = 0.52521e3 #J/Kg.K
Kg_Ar = 0.010926 #W/m.K
# Helium
rhog_He = 0.2778 #kg/m^3
cg_He = 5.1932e3 #J/Kg.K
Kg_He = 0.1069 #W/m.K

pi = 3.14159265359
Kt = 0.316 #tissue(CPA) thermal conductivity W/m.K

Nu = 4.0
d1 = log(rt/rg)
d2_Ar = 2.0*(Kt/Kg_Ar)*Nu
sigma_Ar = (2.0*pi)/(d1+d2_Ar-0.5)

d2_He = 2.0*(Kt/Kg_He)*Nu
sigma_He = (2.0*pi)/(d1+d2_He-0.5)

def TempTissue(sigma,z):
    Tt = Tgi + ((Ttv-Tgi)*(1.0-exp(sigma*z))/(1.0-exp(sigma*L)))
    return Tt

Nr = 1500
Nz = 1000
zlen = np.linspace(0,L,Nz,endpoint=True,dtype=np.double)
rad = rad = np.linspace(rg,rt,Nr,endpoint=True,dtype=np.double)


Tt_Ar = np.zeros(shape=(zlen.size,rad.size))
Tt_He = np.zeros(shape=(zlen.size,rad.size))

for i in range(zlen.size):
        for j in range(rad.size):
            Tt_Ar[i,j] = TempTissue(sigma_Ar,zlen[i])
            Tt_He[i,j] = TempTissue(sigma_He,zlen[i])

X, Y = np.meshgrid(rad,zlen)
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(5,5),dpi=300)
lvl = np.linspace(Tgi,Ttv,20)
cp_He = ax1.contourf(X,Y,Tt_He,levels=lvl,extend='both',cmap='RdYlBu_r')
cp_Ar = ax2.contourf(X,Y,Tt_Ar,levels=lvl,extend='both',cmap='RdYlBu_r')
fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8,
                    wspace=0.02, hspace=0.02)
cb_ax = fig.add_axes([0.83,0.1,0.02,0.8])
cbar = fig.colorbar(cp_He,cax=cb_ax)
plt.show()


