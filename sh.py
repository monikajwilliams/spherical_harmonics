#!/usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt

# => Coordinate transform <= #

def cart_to_sphere(
    x,
    y,
    z,
    ):

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z/r)
    theta[r==0.0] = 0.0
    phi = np.arctan2(y,x)

    return r,theta,phi

def sphere_to_cart(
    r,
    theta,
    phi,
    ):

    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)

    return x,y,z

# => Spherical transforms <= #

def associated_legendre(
    t,  # cos(theta) 
    tc, # square root of (1-t^2) is  (sin(theta))
    L,
    ):

    df = lambda N:np.prod([x for x in range(1,N+1,2)])
    P = {}
    for m in range(L+1):
        P[(m,m)] = df(2*m-1)*(tc**m)

    for m in range(L):
        P[(m+1,m)] = t*(2*m + 1)*P[(m,m)]
        for l in range(m+2,L+1):
            P[(l,m)] = ((2*l - 1)*t*P[(l-1,m)] - (l+m-1)*P[(l-2,m)])/(l-m)

    return P
        

def sh(
    theta,
    phi,
    L, # maximum angular momentum desired
    racah=True, # normalization convention
    ):

    P = associated_legendre(
        np.cos(theta),
        np.sin(theta),
        L,
        )

    Y = {}
    for l in range(L+1):    
        pref = np.sqrt((2.0*l + 1.0)/(4.0*np.pi)) if racah else 1.0 # prefactor
        Y[(l,0)] = pref*P[(l,0)] 
        for m in range(1,l+1):
            pref2 = np.sqrt(np.prod([1.0/x for x in range(l-m+1,l+m+1)]))
            Y[(l,m)] = np.sqrt(2.0)*pref*pref2*P[(l,m)]*np.cos(m*phi)
            Y[(l,-m)] = np.sqrt(2.0)*pref*pref2*P[(l,m)]*np.sin(m*phi)

    return Y
            
        
def plot_sh(
    fn_out,
    L = 20,
    lm = (10,10)
    ):

    theta = np.linspace(0.0,np.pi,100)
    phi = np.linspace(0.0,2.0*np.pi,100)

    tt,pp = np.meshgrid(theta,phi)

    Y = sh(tt,pp,L)

    print Y.keys()
    
    plt.clf()
    plt.contourf(pp.T,tt.T,Y[lm].T,cmap='seismic')
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\theta$')
    plt.colorbar()
    plt.savefig(fn_out)
    
#plot_sh('test.pdf', L=20,lm=(10,5)) 
