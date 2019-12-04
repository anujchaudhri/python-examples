#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:58:14 2019

@author: anujchaudhri
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

Vo = 5.24e-13 # m^3, initial cell volume Vc
Vso = 0.0 # m^3, initial solute volume in cell Vs
vbf = 0.25
Vb = vbf*Vo
# Vc = Vw + Vs + Vb, total volume of cell at time t
# Vo = Vwo + Vso + Vbo, total volume of cell at time 0
# Vbo = Vb, since cell solids volume does not change
Vwo = Vo - (Vso + Vb)# m^3, initial water volume in cell Vw

# function that returns dz/dt = [dx/dt dy/dt]
def model(z,t):
    Vsbar = 0.071e-3 # m^3/mol
    Mino = 0.3 # osmoles/Kg
    Lp = 1.0*(1e-6/60.0) # m/(sec*atm)
    A = 3.14e-8 # m^2
    R = 8.2057338e-5 # m^3*atm/(K*mol)
    T = 293 # K
    Mes = 2.0 # osmoles/Kg
    Men = 0.3 # osmoles/Kg
    Ps = 0.03*(1e-2/60.0) # m/sec
    
    A1 = (((Lp*A*R*T)/z[0])*((z[1]/Vsbar)+(Mino*Vwo))) - ((Lp*A*R*T)*(Mes+Men))
    A2 = ((Ps*A*Vsbar*Vsbar/2.0)*(Mes+(z[1]/(z[0]*Vsbar)))
            *(Mes+Men-(z[1]/(z[0]*Vsbar))-(Mino*Vwo/z[0])))
    A3 = (((Ps*Ps*A*Vsbar*Vsbar*Vsbar)/(2.0*Lp*R*T))
    *((Mes*Mes)-((z[1]/(z[0]*Vsbar))*(z[1]/(z[0]*Vsbar)))))
    A4 = -(Ps*A*z[1]/z[0]) + (Ps*A*Vsbar*Mes)
    
    dxdt = A1+A2-A3
    dydt = A4-A2+A3
    dzdt = [dxdt,dydt]
    return dzdt

# initial conditions
z0 = [Vwo,Vso]

# time points in secs
ty = 3600
N = 3600
t = np.linspace(0,ty,N,endpoint=True,dtype=np.int)
#print(t)

# solve ODE
z = odeint(model,z0,t)

# Note that z[:,0] is solution at all time points of dxdt i.e. Vw
# z[:,1] is solution at all time points of dydt i.e. Vs
# Vc = Vw+Vs+Vb and Vn = (Vw+Vs+Vb)/Vo or Vn = (z[:,0]+z[:,1]+Vb)/Vo

Vc = np.add(z[:,0],z[:,1])
Vbarray = np.full_like(Vc,Vb,dtype=np.double)
Vc = np.add(Vc,Vbarray)
Vn = Vc/Vo

#plt.axis([-50,ty,0.9999,1])
plt.plot(t,Vn,'g-') 
plt.show()


