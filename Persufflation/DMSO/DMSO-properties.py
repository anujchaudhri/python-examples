#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 10:12:34 2019

@author: anujchaudhri
"""

import numpy as np
import matplotlib.pyplot as plt

# properties of 7.05 M DMSO

# Density

def density(T):
    den = 1058 - (0.41*T)
    return den

def specificheat(T):
    spheat = 2804 +(4.205*T)-(0.054*T*T)-(4.902*T*T*T)
    return spheat

def thermalcond(T):
    thcond = (0.356 + (7.42e-4*T) - (1.29e-6*T*T) - (6.87e-8*T*T*T) - 
              (2.95e-10*T*T*T*T))
    return thcond

Ti = -20
Tf = -135
N = 116
Temp = np.linspace(Ti,Tf,N,endpoint=True,dtype=np.double)

den = np.zeros(shape=(Temp.size))
spheat = np.zeros(shape=(Temp.size))
thcond = np.zeros(shape=(Temp.size))
alpha = np.zeros(shape=(Temp.size))
t = np.zeros(shape=(Temp.size))
dia = 100e-6 # metres, endothelial cell wall thickness
R0 = dia/2.0

for i in range(Temp.size):
    den[i] = density(Temp[i])
    spheat[i] = specificheat(Temp[i])
    thcond[i] = thermalcond(Temp[i])
    alpha[i] = thcond[i]/(den[i]*spheat[i])
    t[i] = (R0*R0)/alpha[i]
    
fig1 = plt.figure(figsize=(5,5),dpi=300)
plt.plot(Temp,t,'bo-')
plt.xlabel('Temperature (C)')
plt.ylabel('Thermal Diffusion Time (s)')
plt.title('7.05M DMSO, 100 microns length')
plt.show()
#fig2 = plt.figure(figsize=(4,4),dpi=100)
#plt.plot(Temp,t,'ro-')
#plt.show()