#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:29:18 2019

@author: anujchaudhri
"""

import numpy as np
import matplotlib.pyplot as plt

Ttv = 269.15
Tgi = 173.15

rg = 1500e-6
L = 100e-3

rhog = 0.2778 #kg/m^3
vdotg = 20e-3/60.0 #m^3/sec
mdot = rhog*vdotg #Kg/s
cg = 5.1932e3 #J/Kg.K
Kg = 0.1069 #W/m.K
pi = 3.14159265359
Ac = pi*rg*rg #m^2
UmbyAlpha = (mdot*cg)/(Kg*Ac)
P = 2.0*pi*rg #m
Kt = 0.316 #tissue(CPA) thermal conductivity W/m.K

def Tempgas(r,z):
    dTbyL = (Ttv-Tgi)/L
    heatflux = -Kt*dTbyL
    dTmbydz = heatflux*P/(mdot*cg)
    rbyrg = (3/16)+((1/16)*(r/rg)*(r/rg)*(r/rg)*(r/rg))-((1/4)*(r/rg)*(r/rg))
    R1 = dTbyL*z
    R2 = Tgi
    R3 = (2.0*UmbyAlpha*rg*rg)*dTmbydz*rbyrg
    Tgas = R1+R2-R3
    return Tgas

Nr = 1500
Nz = 1000
rad = np.linspace(0,rg,Nr,endpoint=True,dtype=np.double)
zlen = np.linspace(0,L,Nz,endpoint=True,dtype=np.double)

#print(rad)
#print(zlen)

Tg = np.zeros(shape=(zlen.size,rad.size))

for i in range(zlen.size):
    for j in range(rad.size):
        Tga = Tempgas(rad[j],zlen[i])
        Tg[i,j] = Tga
    
print(Tg)
X, Y = np.meshgrid(rad,zlen)
fig = plt.figure(figsize=(5,5),dpi=300)
lvl = np.linspace(Tgi,Ttv,20)
cp = plt.contourf(X,Y,Tg,levels=lvl,extend='both',cmap='RdYlBu_r')
cbar = plt.colorbar(cp)
plt.show()

