# -*- coding: utf-8 -*-
"""
Created on Thu Dec 25 21:21:22 2014

@author: Ronald.Dai
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
def integrate_RK5(F,x,y,xstop,h,tol=1e-6):
    
    def run_kut5(F,x,y,h):
        C = np.array([37/378,0,250/621,125/594,0,512/1771])
        D = np.array([2825/27648,0,18575/48384,13525/55296,277/14336,1/4])
        n = len(Y)
        K = np.zeros((6,n))
        K[0] = h*F(x,y)
        K[1] = h*F(x+h/5,y+K[0]/5)
        K[2] = h*F(x+3*h/10,y+3*K[0]/40+9*K[1]/40)
        K[3] = h*F(x+3*h/5,y+3*K[0]/10-9*K[1]/10+6*K[2]/5)
        K[4] = h*F(x+h,y-11*K[0]/54+5*K[1]/2-70*K[2]/27+35*K[3]/27)
        K[5] = h*F(x+7*h/8,y+1631*K[0]/55296+175*K[1]/512+575*K[2]/13824+44275*K[3]/110592+253*K[4]/4096)
        
        E = np.zeros(n)
        dy = np.zeros(n)
        for i in range(6):
            dy = dy + C[i]*K[i]
            E = E + (C[i]-D[i])*K[i]
            e = np.sqrt(sum(E**2)/n)
        return dy,e
        
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    stopper = 0
    for i in range(10000):
        dy,e = run_kut5(F,x,y,h)
        if e <= tol:
            y = y+dy
            x = x+h
            X.append(x)
            Y.append(y)
            if stopper ==1: break
        if e != 0:
            hNext = 0.9*h*(tol/e)**0.2
        else hNext = h
        
        if (h>0) == ((x+hNext)>=xstop):
           hNext = xstop-x
           stopper = 1
        h = hNext
    return np.array(X),np.array(Y)
    
