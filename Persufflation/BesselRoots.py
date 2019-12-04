#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 09:53:35 2019

@author: anujchaudhri
"""

import numpy as np
from scipy.special import j0,j1,y0,y1,jn,yn
from scipy import optimize
import matplotlib.pyplot as plt

def Jp(m,r):
    jp = 0.5*(jn(m-1,r)-jn(m+1,r))
    return jp

def Yp(m,r):
    yp = 0.5*(yn(m-1,r)-yn(m+1,r))
    return yp

def phi(lam,r):
    phir = (j0(lam*r)*y1(lam))-(j1(lam)*y0(lam*r))
    return phir
    
def dphi(lam,r):
    dphir = -lam*((j1(lam*r)*y1(lam))-(j1(lam)*y1(lam*r)))
    return dphir

def ddphi(lam,r):
    ddphir = -lam*((Jp(1,lam*r)*y1(lam))-(j1(lam)*Yp(1,lam*r)))
    return ddphir

def eigencond(lam,Bi,a):
    eigeneqn = dphi(lam,a) + (Bi*phi(lam,a))
    return eigeneqn

def deigencond(lam,Bi,a):
    deigeneqn = ddphi(lam,a) + (Bi*dphi(lam,a))
    return deigeneqn

def eigenroot(lamstart,Bi,a):
    dlam = 0.1
    lam0 = lamstart+1e-6
    eigen0 = eigencond(lam0,Bi,a)
    lam0 = lam0 + dlam
    eigen1 = eigencond(lam0,Bi,a)
    while(eigen0*eigen1 > 0.0):
        lam0 = lam0 + dlam
        eigen0 = eigen1
        eigen1 = eigencond(lam0,Bi,a)
    lamr = lam0 - ((eigen1*dlam)/(eigen1-eigen0))
    #lamroot = (optimize.root_scalar(eigencond,x0=lamr,fprime=deigencond,
    #                               method='newton'))
    #return lamroot.root
    return lamr

#lamroot1 = (optimize.root_scalar(eigencond,x0=5,fprime=deigencond,
#                                   method='newton'))
#lamroot2 = optimize.root_scalar(eigencond,bracket=[25,30],method='bisect')

#x = np.linspace(0,100,100)
#plt.plot(x,eigencond(x),'ro-')
#
#nroots = 15
#lamroots = np.zeros(shape=(nroots))
#start_point = 0.0
#increment = 2.0
#ctr = 0
#
#while (ctr < nroots):
#    lroots = eigenroot(start_point)
#    if round(lroots,7) in np.round(lamroots[:],7):
#        start_point += increment
#    else:
#        print(lroots)
#        lamroots[ctr] = lroots
#        ctr += 1        
#
#print(lamroots)
#print(lamroot1.root)
#print(lamroot2.root)