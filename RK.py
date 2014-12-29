# -*- coding: utf-8 -*-
"""
Created on Thu Dec 25 21:21:22 2014

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt\

"""
the runge-kutta method with 4order
"""
def integrate_RK4(F,x,y,xstop,h):
    
    def run_kut4(F,x,y,h):       
        K0=h*F(x,y)
        K1=h*F(x+h/2,y+K0/2)
        K2=h*F(x+h/2,y+K1/2)
        K3=h*F(x+h,y+K2)
        return (K0+2*K1+2*K2+K3)/6      
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    while x<xstop:
        h = min(h,xstop-x)
        y = y + run_kut4(F,x,y,h)
        x = x + h
        X.append(x)
        Y.append(y)
    return np.array(X),np.array(Y)

"""
the adaptive runge-kutta method with 4order working with cash-karp coefficients
"""
def integrate_RK5(F,x,y,xstop,h):
    
    def run_kut5(F,x,y,h):
        C = np.array([])
    
