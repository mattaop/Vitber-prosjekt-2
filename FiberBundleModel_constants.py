import numpy as np
from random import randint
N = 4
g = 9,81
#beta = 9
#alpha = 5.3e-9
#alfa = 1
b_k = 2

def B(r,E):
    return (E*np.pi*r**4/4)
    #plotSolution(x)

def alpha(g, Lambda, l_between_wire, B):
    return g*Lambda*(l_between_wire**3)/B

def kappa(E_wire,r_wire,l_wire):
    return E_wire*np.pi*r_wire**2/l_wire

def beta(E_wire,r_wire, l_wire, l_between_wire, B,systemSize):


    beta=kappa(E_wire,r_wire, l_wire)*(l_between_wire**3)/B
    return beta


