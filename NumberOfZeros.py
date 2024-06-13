import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
mp.mp.dps=15

# Choose T
T = 10**(9)
def A(T):
    return T/(2*mp.pi)*mp.log(T/(2*mp.pi)) - T/(2*mp.pi) + 7/8

def B(T):
    a1,a2,a3 = 0.122,0.278,2.510
    return a1*mp.log(T) + a2*mp.log(mp.log(T)) + a3

upper1 = A(T+1) - A(T) + B(T+1) + B(T)

print(upper1)
