#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:58:17 2019

@author: anujchaudhri
"""


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

Vo = 5.24e5 # microm^3, initial cell volume Vc
Vso = 0.0 # microm^3, initial solute volume in cell Vs
vbf = 0.25
Vb = vbf*Vo
# Vc = Vw + Vs + Vb, total volume of cell at time t
# Vo = Vwo + Vso + Vbo, total volume of cell at time 0
# Vbo = Vb, since cell solids volume does not change
Vwo = Vo - (Vso + Vb)# microm^3, initial water volume in cell Vw

rhoWater = 997.0e-18 # Kg/microm^3
Vsbar = 0.071e15 # microm^3/mol
Mino = 0.3*rhoWater # osmoles/Kg * Kg/m^3
Lp = 1.0 # microm/(min*atm)
A = 3.14e4 # microm^2
R = 8.2057338e13 # microm^3*atm/(K*mol)
T = 293.0 # K
Mes = 2.0*rhoWater # osmoles/Kg * Kg/m^3
Men = 0.3*rhoWater # osmoles/Kg * Kg/m^3
Ps = 0.01e4 # microm/min
Psarray = np.array([0.0e4,0.001e4,0.01e4,0.03e4]) # microm/min


    
# function that returns dz/dt = [dx/dt dy/dt]
def model(t,z,Ps):
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
tf = 4
N = 50
t = np.linspace(0,tf,N,endpoint=True,dtype=np.double)
#print(t)
Vn = np.zeros(shape=(t.size,Psarray.size))

# solve ODE
for i in range(Psarray.size):
    sigma = 1.0-((Psarray[i]*Vsbar)/(R*T*Lp))
#    print("sigma: ",sigma)
    
    sol = solve_ivp(fun=lambda t, z: model(t, z, Psarray[i]),t_span=[0,tf]
                    ,y0=z0,method='BDF',t_eval=t,vectorized=True)
    
    # Note that z[:,0] is solution at all time points of dxdt i.e. Vw
    # z[:,1] is solution at all time points of dydt i.e. Vs
    # Vc = Vw+Vs+Vb and Vn = (Vw+Vs+Vb)/Vo or Vn = (z[:,0]+z[:,1]+Vb)/Vo
    
    #print("status: explanations")
    #print("-1: Integration step failed.")
    #print("0: The solver successfully reached the end of tspan.")
    #print("1: A termination event occurred. ")
    #print("status: ", sol.status)
    ##print(sol.y)
    #print("Number of function evaluations: ", sol.nfev)
    #print("Number of Jacobian evaluations: ", sol.njev)
    
    Vc = np.add(sol.y[0,:],sol.y[1,:])
    Vbarray = np.full_like(Vc,Vb,dtype=np.double)
    Vc = np.add(Vc,Vbarray)
    Vn[:,i] = Vc/Vo

print(Vn)

#plt.axis([-0.2,tf,0.0,1.2])
#plt.plot(t,Vn[:,0],'ro',t,Vn[:,1],'go',t,Vn[:,2],'bo',t,Vn[:,3],'ko') 
#plt.show()

