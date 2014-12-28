"""
Created on Sun Dec 28 18:42:18 2014

@author: Administrator
"""
"""
H = h^j/j!
j=1: H1 =h/1
j=2: H2 = h^2/(1x2)
j=3: H3 = h^3/(1x2x3)=(h/3)*(h^2/(1x2))
"""
import numpy as np

def taylor(deriv,x,y,xstop,h):
    X = []
    Y = []
    X.append(x)
    Y.append(y)
    
    while x < xstop:
        h = min(h,xstop-x)
        D = deriv(x,y)
        H = 1.0
        for j in range(4):
            H = H*h/(j+1)
            y = y+D[j]*H
            x = x + h
            X.append(x)
            Y.append(y)
    return np.array(X), np.array(Y)
