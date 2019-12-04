#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:00:42 2019

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
Mesarray = np.array([2.0,4.0,10.0])*rhoWater # osmoles/Kg * Kg/m^3
Men = 0.3*rhoWater # osmoles/Kg * Kg/m^3
Psarray = np.array([0.01e4]) # microm/min

# function that returns dz/dt = [dx/dt dy/dt]
def model2P(t,z,Ps,Mes):
    dxdt = ((((Lp*A*R*T)/z[0])*((z[1]/Vsbar)+(Mino*Vwo))) 
                - ((Lp*A*R*T)*(Mes+Men)))
    dydt = -(Ps*A*z[1]/z[0]) + (Ps*A*Vsbar*Mes)
    dzdt = [dxdt,dydt]
    return dzdt

# function that returns dz/dt = [dx/dt dy/dt]
def modelKK(t,z,Ps,Mes):
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


def jacobian2P(t,y,Ps):
    J11 = -(Lp*A*R*T/(y[0]*y[0]))*((y[1]/Vsbar)+(Mino*Vwo))
    J12 = (Lp*A*R*T)/(y[0]*Vsbar)
    J21 = (y[1]*Ps*A)/(y[0]*y[0])
    J22 = -(Ps*A)/y[0]
    Jacob = np.array([[J11,J12],[J21,J22]])
    return Jacob

# initial conditions
z0 = [Vwo,Vso]

# time points in mins
tf = 8
N = 50
t = np.linspace(0,tf,N,endpoint=True,dtype=np.double)
#print(t)
Vn2P = np.zeros(shape=(t.size,Mesarray.size))
VnKK = np.zeros(shape=(t.size,Mesarray.size))

# solve ODE
for i in range(Mesarray.size):
    
    sol2P = (solve_ivp(fun=lambda t, z: model2P(t, z, Psarray[0], Mesarray[i]),t_span=[0,tf],y0=z0,
                     method='BDF',t_eval=t,vectorized=True,
                     jac=lambda t, y: jacobian2P(t, y, Psarray[0])))
    solKK = solve_ivp(fun=lambda t, z: modelKK(t, z, Psarray[0], Mesarray[i]),t_span=[0,tf]
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
    
    Vc2P = np.add(sol2P.y[0,:],sol2P.y[1,:])
    Vbarray = np.full_like(Vc2P,Vb,dtype=np.double)
    Vc2P = np.add(Vc2P,Vbarray)
    Vn2P[:,i] = Vc2P/Vo
    
    VcKK = np.add(solKK.y[0,:],solKK.y[1,:])
    Vbarray = np.full_like(VcKK,Vb,dtype=np.double)
    VcKK = np.add(VcKK,Vbarray)
    VnKK[:,i] = VcKK/Vo


#print(Vn2P)
#print(VnKK)

plt.figure(figsize=(5,5),dpi=300)
plt.axis([-0.2,8.5,0.5,1.5])
plt.xlabel("time (min)")
plt.ylabel("Vn")

plt.plot(t,Vn2P[:,0],'r-',label='2P,M=2.0')
plt.plot(t,VnKK[:,0],'ro',label='KK,M=2.0')

plt.plot(t,Vn2P[:,1],'g-',label='2P,M=4.0')
plt.plot(t,VnKK[:,1],'go',label='KK,M=4.0')

plt.plot(t,Vn2P[:,2],'b-',label='2P,M=10.0')
plt.plot(t,VnKK[:,2],'bo',label='KK,M=10.0')

plt.legend()
plt.show()