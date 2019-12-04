#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:12:17 2019

@author: anujchaudhri
"""
import numpy as np

Vo = 5.24e-13 # m^3, initial cell volume Vc
Vso = 0.0 # m^3, initial solute volume in cell Vs
vbf = 0.25
Vb = vbf*Vo
# Vc = Vw + Vs + Vb, total volume of cell at time t
# Vo = Vwo + Vso + Vbo, total volume of cell at time 0
# Vbo = Vb, since cell solids volume does not change
Vwo = Vo - (Vso + Vb)# m^3, initial water volume in cell Vw
Vsbar = 0.071e-3 # m^3/mol
Mino = 0.3 # osmoles/Kg
Lp = 1.0*(1e-6/60.0) # m/(sec*atm)
A = 3.14e-8 # m^2
R = 8.2057338e-5 # m^3*atm/(K*mol)
T = 293 # K
Mes = 2.0 # osmoles/Kg
Men = 0.3 # osmoles/Kg
Ps = 0.03*(1e-2/60.0) # m/sec
sigma = 1.0 - ((Ps*Vsbar)/(R*T*Lp))
print(sigma)

x = (Mino*Vwo/Men)
y = (Vsbar*Mes*x)
print(x)
print(y)

J11 = -(Lp*A*R*T/(x*x))*((y/Vsbar)+(Mino*Vwo))
J12 = (Lp*A*R*T)/(x*Vsbar)
J21 = (y*Ps*A)/(x*x)
J22 = -(Ps*A)/x

J = np.matrix([[J11,J12],[J21,J22]])
print(J)

eigenvalues = np.linalg.eigvals(J)
print(eigenvalues)
print("Condition Number: ", abs(eigenvalues[1]/eigenvalues[0]))