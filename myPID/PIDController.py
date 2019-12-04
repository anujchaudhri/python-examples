#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

PID Controller
Anuj Chaudhri
Arigos Biomedical
Sep 2018

"""
from datetime import datetime
import random
from array import array
import matplotlib.pyplot as plt
from PID import PID

#tm = datetime.now()
#Time_now = tm.second
#print (Time_now)

time_now = 0.0
# sample_time is how quickly the PV is sampled by sensor
sample_time = 1.0
# controller_time is how frequently the contoller is called
controller_time = 2.0 
# max_time is maximum time for experiment
max_time = 10000.0

# current_value is current reading of sensor
current_value = 20.0
# fixed_value is the set point
fixed_value = 40.0
# bump_value and bump time are to simulate changes in fixed_value with time
bump_value = 0.0
bump_time = 0.0

# initialize PV_out to current_value initially in case the controller
# is not called right away
PV_out = current_value

_PIDPV = array('d')
_current = array('d')
_time = array('d')
_fixedPV = array('d')

# Initialize the P,I,D values 
_PID = PID(1.0,0.0,0.0,time_now)

_time.append(time_now)
_PIDPV.append(current_value)
_fixedPV.append(fixed_value)

while (time_now < max_time):
    time_now += sample_time
    if(time_now == bump_time):
        fixed_value = current_value + bump_value
        bump_time += bump_time
    if((time_now%controller_time) == 0):
        [P_out, I_out, D_out] = _PID.updatePID(current_value,fixed_value, \
                                                time_now)
        PV_out = P_out + I_out + D_out

#    Since we dont have a real sensor here, the current_value which
#    represents the sensor here is set to the output of the controller
#    PV_out
#    In real cases,the current_value will be the input of a sensor
#    (process variable, PV)
#    after the controller output (CO, here PV_out) is used to make some changes to the
#    process (final control element, FCE) and manipulated variable (MV)
#    i.e. say the speed of the servo motor (FCE) and flow rate (MV)
        
#    current_value = PV_out
        
#   Bump current_value by 20% every time to simulate change in PV
#   i.e. sensor reading by changing the manipulated variable 
#   like say the speed of a servo motor
    current_value = current_value + (0.2*current_value)
    print(PV_out)
    print(current_value)
    _PIDPV.append(current_value)
    _time.append(time_now)
    _fixedPV.append(fixed_value)

    
plt.plot(_time,_fixedPV,'b-',_time,_PIDPV,'ro-')
