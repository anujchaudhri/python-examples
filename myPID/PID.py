#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: anujchaudhri
"""

# Define a class for PID controls

class PID:
    
    def __init__(self, P_gain, I_gain, D_gain, now):
        
        self.previous_error = 0.0
        self.last_time = now      
        self.P_gain = P_gain
        self.I_gain = I_gain
        self.D_gain = D_gain      
        self.I_error = 0.0
        
    def updatePID(self,current_value,fixed_value,now):
        dt = (now - self.last_time)
#        dt = now
        error = fixed_value - current_value          
        P_error = error
        
        self.I_error += error    
        I_error = self.I_error
        
        D_error = (error - self.previous_error)
        
        P_out = P_error * self.P_gain
        I_out = I_error * self.I_gain * dt
        if(dt != 0.0):
            D_out = D_error * (self.D_gain / dt)
#        print(error)
        self.previous_error = error
        self.last_time = now
     
        return P_out,I_out,D_out
    