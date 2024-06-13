import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
import random
import math

mp.mp.dps = 100
mp.pretty = False
TOL = 1.e-100
error=0.7

window=0

alphaR=np.array([error*(2*random.random()-1)+1 for ii in range(10001)])
alphaI=np.array([error*(2*random.random()-1) for ii in range(10001)])

"""
Line 27 is for Newton's method or Relaxed Newton's method by changing alphaC between 0 and 1
Line 28 is for Random Relaxed Newton's method
"""
def newton(z0, f, fprime, MAX_IT=10000):
    z = z0
    for i in range(MAX_IT):
        dz = f(z) / fprime(z)
        if abs(f(z)) < TOL:
            break
        alphaC=1
        #alphaC=alphaR[i]+alphaI[i]*1j
        z -= alphaC*dz
    return z

def plot_newton_fractal(f, fprime, n=120, domain=(-6, 6, -6, 6)):
    m = np.array([[0.01 for _ in range(-n, n+1)] for _ in range(-n, n+1)])
    v = np.array([0.5+0.1*np.random.rand(),window+0.1*np.random.rand()])

    stepstepx=0.5/n
    stepstepy=4.5/n
    count1 = 0
    count2 = 0

    for ix in range(-n, n+1):
        for iy in range(-n, n+1):
            x = v[0] + ix * stepstepx
            y = 10*(v[1] + iy*stepstepy)
            z = mp.mpc(x,y)
            r = newton(z, f, fprime)
            print("Iterate=", )
            print("The initial point=", z)
            print("The last point=", r)
            print("Value of zeta at the last point=", f(r))

            if abs(r - Root1) < 0.0001:
                m[ix,iy]=0
            elif abs(r - Root2) < 0.0001:
                m[ix,iy]=1
            elif abs(r - Root3) < 0.0001:
                m[ix, iy]=2
            elif abs(r - Root4) < 0.0001:
                m[ix, iy] = 3
            elif abs(r - Root5) < 0.0001:
                m[ix, iy] = 4
            elif abs(r - Root6) < 0.0001:
                m[ix, iy] = 5
            elif abs(r - Root7) < 0.0001:
                m[ix, iy] = 6
            elif abs(r - Root8) < 0.0001:
                m[ix, iy] = 7
            else:
                m[ix, iy] = 8


    def __arrow__(x, y, dx, dy, width, length):
        plt.arrow(x, y, dx, dy,
            color='k',
            clip_on=False,
            head_width=width,
            head_length=length
        )

    xlim = (-1, 2)
    ylim= (window-5,window+5)
    figsize = (15, 15)
    fig, ax = plt.subplots(figsize=figsize)
    plt.xlim(xlim)
    plt.ylim(ylim)

    spacing = 1
    minorLocator = plt.MultipleLocator(spacing)
    plt.gca().yaxis.set_minor_locator(minorLocator)
    plt.gca().xaxis.set_minor_locator(minorLocator)

    spacing_major = 1
    majorLocator = plt.MultipleLocator(spacing_major)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().yaxis.set_major_locator(majorLocator)

    plt.grid(which='both')
    plt.gca().set_aspect("equal")

    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    __arrow__(xlim[1], 0, 0.01, 0, 0.3, 0.2)  # x-axis arrow
    __arrow__(0, ylim[1], 0, 0.01, 0.2, 0.3)  # y-axis arrow

    for ix in range(-n, n+1):
        for iy in range(-n, n+1):
            if m[ix, iy] == 0:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="green")
            elif m[ix, iy] == 1:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="yellow")
            elif m[ix, iy] == 2:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="blue")
            elif m[ix, iy] == 3:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="red")
            elif m[ix, iy] == 4:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="pink")
            elif m[ix, iy] == 5:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="cyan")
            elif m[ix, iy] == 6:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="orange")
            elif m[ix, iy] == 7:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="purple")
            else:
                plt.plot(v[0] + ix * stepstepx, v[1] + iy * stepstepy, marker=".", color="black")
    plt.show()
    plt.close()

"""
The Riemann xi function
"""
def pol(z):
    tol=0.5*z*(z-1)*mp.pi**(-0.5*z)*mp.gamma(0.5*z)*mp.zeta(z)
    #tol = (z-Root1)*(z-Root2)*(z-Root3)*(z-Root4)*(z-Root5)*(z-Root6)*(z-Root7)*(z-Root8)
    return tol

def polDer(z):
    return mp.diff(pol, z)


Root1=0.5+14.13472514173*1j
Root2=0.5-14.13472514173*1j
Root3=0.5+21.02203963877*1j
Root4=0.5-21.02203963877*1j
Root5=0.5+25.01085758014*1j
Root6=0.5-25.01085758014*1j
Root7=0.5+30.42487612585*1j
Root8=0.5-30.42487612585*1j

f = lambda z: pol(z)
fprime = lambda z: polDer(z)

plot_newton_fractal(f, fprime, n=120)
