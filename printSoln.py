"""
Created on Sun Dec 28 18:49:17 2014

@author: Administrator
"""
import matplotlib.pyplot as plt

def printSoln(X,Y,freq):
    
    def printHead(n):
        print "\n x",
        for i in range(n):
            print "y[",i,"]",
        print
    
    def printLine(x,y,n):
        print "%13.4e"% x,
        for i in range(n):
            print "%13.4e"% y[i],
        print
        
    m = len(Y)
    try: n = len(Y[0])
    except TypeError: n =1
    if freq == 0: freq = m
    printHead(n)
    for i in range(0,m,freq):
        printLine(X[i],Y[i],n)
    if i != m -1: 
        printLine(X[m-1],Y[m-1],n)
        plt.plot(X,Y,'*')
