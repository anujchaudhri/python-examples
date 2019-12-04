#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 12:21:38 2019

@author: anujchaudhri
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Vo = 5.24e-13 # m^3, initial cell volume Vc
Vso = 0.0 # m^3, initial solute volume in cell Vs
vbf = 0.25
Vb = vbf*Vo
# Vc = Vw + Vs + Vb, total volume of cell at time t
# Vo = Vwo + Vso + Vbo, total volume of cell at time 0
# Vbo = Vb, since cell solids volume does not change
Vwo = Vo - (Vso + Vb)# m^3, initial water volume in cell Vw

rhoWater = 997.0 # Kg/m^3
Vsbar = 0.071e-3 # m^3/mol
Mino = 0.3*rhoWater # osmoles/Kg * Kg/m^3
Lp = 1.0*(1e-6/60.0) # m/(sec*atm)
A = 3.14e-8 # m^2
R = 8.2057338e-5 # m^3*atm/(K*mol)
T = 293.0 # K
Mes = 2.0*rhoWater # osmoles/Kg * Kg/m^3
Men = 0.3*rhoWater # osmoles/Kg * Kg/m^3
Ps = 0.0*(1e-2/60.0) # m/sec

# function that returns dz/dt = [dx/dt dy/dt]
def model(t,z):
    dxdt = ((((Lp*A*R*T)/z[0])*((z[1]/Vsbar)+(Mino*Vwo))) 
                - ((Lp*A*R*T)*(Mes+Men)))
    dydt = -(Ps*A*z[1]/z[0]) + (Ps*A*Vsbar*Mes)
    dzdt = [dxdt,dydt]
    return dzdt

def jacobian(t,y):
    J11 = -(Lp*A*R*T/(y[0]*y[0]))*((y[1]/Vsbar)+(Mino*Vwo))
    J12 = (Lp*A*R*T)/(y[0]*Vsbar)
    J21 = (y[1]*Ps*A)/(y[0]*y[0])
    J22 = -(Ps*A)/y[0]
    Jacob = np.array([[J11,J12],[J21,J22]])
    return Jacob

# initial conditions
z0 = [Vwo,Vso]

# time points in secs
tf = 240
N = 240
t = np.linspace(1,tf,N,endpoint=True,dtype=np.int)
#print(t)

# solve ODE
sol = (solve_ivp(model,[1,tf],z0,method='BDF',t_eval=t,
                 vectorized=True,jac=jacobian))

# Note that z[:,0] is solution at all time points of dxdt i.e. Vw
# z[:,1] is solution at all time points of dydt i.e. Vs
# Vc = Vw+Vs+Vb and Vn = (Vw+Vs+Vb)/Vo or Vn = (z[:,0]+z[:,1]+Vb)/Vo

print("status: explanations")
print("-1: Integration step failed.")
print("0: The solver successfully reached the end of tspan.")
print("1: A termination event occurred. ")
print("status: ", sol.status)
#print(sol.y)
print("Number of function evaluations: ", sol.nfev)
print("Number of Jacobian evaluations: ", sol.njev)

Vc = np.add(sol.y[0,:],sol.y[1,:])
Vbarray = np.full_like(Vc,Vb,dtype=np.double)
Vc = np.add(Vc,Vbarray)
Vn = Vc/Vo

#plt.axis([-50,tf,0.9999,1])
plt.plot(t,Vn,'g-') 
plt.show()
