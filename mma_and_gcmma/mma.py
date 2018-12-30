from sympy import *


x1=Symbol('x1')
x2=Symbol('x2')


def pqr(f, curval, preval, ppval, xmin, xmax, k):
    l = []
    u = []
    gamma = []
    for i in range(len(curval)):
        gamma.append((curval[i]-preval[i])*(preval[i]-ppval[i]))
        if gamma[i]<0:
            gamma[i]=0.7
        elif gamma[i]>0:
            gamma[i]=1.2
        else:
            gamma[i]=1
    
    if k<3 
    for i in range(len(curval)):
        l.append(curval[i]-0.5*(xmax[i]-xmin[i]))
        u.append(curval[i]+0.5*(xmax[i]-xmin[i]))
    





f=x1**2+x2**2
xmin = (-1,-2)
xmax = (1,2) 