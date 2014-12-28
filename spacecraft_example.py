"""
Created on Sun Dec 28 20:19:47 2014

@author: Administrator
"""

"""
A spacecraft is launched at the altutide H = 772km above the sea level with the 
speed v = 6700 m/s in the direction shown, the differential equations descrbing the 
motion of teh spacecraft are
       
       r(2) = r*theta(1)^2-GM/r^2
       theta(2) = -2*r(1)*theta(1)/r
G = 6.627x10^-11
M = 5.974210^24
R = 6378.14km

y=[y0,y1,y2,y3]=[r,r(1),theta,theta(1)]
y(1) = [y0(1),y1(1),y2(1),y3(1)]=[r(1),r(2),theta(1),theta(2)]=[y1,y0*y3^2-GM/y0^2,y3,-2*y1*y3/y0]

"""

import numpy as np
import run_kut4 as rk4
import printSoln as psl
import matplotlib.pyplot as plt

def F(x,y):
    F = np.zeros(4)
    F[0] = y[1]
    F[1] = y[0]*(y[3]**2)-3.986e14/(y[0]**2)
    F[2] = y[3]
    F[3] = -2*y[1]*y[3]/y[0]
    return F

x = 0
xstop = 1200
y = np.array([7.15014e6,0,0,0.937045e-3])
h = 50
freq = 2

X,Y = rk4.integrate(F,x,y,xstop,h)
psl.printSoln(X,Y,freq)
