#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 11:40:57 2019

@author: anujchaudhri
"""

import numpy as np
import matplotlib.pyplot as plt
import BesselRoots as br

T1 = -4
T2 = -80
Nu = 4.0
a = 0.5
kHe = 0.1152 # at -80 C
kCPA = 0.310 # at -80 C
Bi = (Nu/(2.0*a))*(kHe/kCPA)
nroots = 40
lamroots = np.zeros(shape=(nroots))
start_point = 0.0
increment = 2.0
ctr = 0

def An(lam,r):
    numer = - 2.0*Bi*a*br.phi(lam,a)
    denom1 = lam*lam*br.phi(lam,1)*br.phi(lam,1)
    denom2 = lam*lam*a*a*br.phi(lam,a)*br.phi(lam,a)*(1.0+((Bi*Bi)/(lam*lam)))
    A = numer/(denom1-denom2)
    return A

def exptime(lam,t):
    Et = np.exp(-lam*lam*t)
    return Et

while (ctr < nroots):
    lroots = br.eigenroot(start_point,Bi,a)
    if round(lroots,7) in np.round(lamroots[:],7):
        start_point += increment
    else:
        lamroots[ctr] = lroots
        ctr += 1        

N = 100
Nt = 11
tend = 1e-1
r = np.linspace(a,1,N,endpoint=True,dtype=np.double)
t = np.linspace(0,tend,Nt,endpoint=True,dtype=np.double)
Temp = np.zeros(shape=(t.size,r.size))

for i in range(t.size):
    for j in range(r.size):
        for k in range(lamroots.size):
            term = (An(lamroots[k],r[j])*br.phi(lamroots[k],r[j])
                                        *exptime(lamroots[k],t[i]))
            Temp[i][j] += term

for i in range(t.size):        
    plt.plot(r,(Temp[i,:]*(T1-T2))+T2,'ro')
            
    