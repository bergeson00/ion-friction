# -*- coding: utf-8 -*-
"""
test004.py

Created by SDB on 6/25/2018

testing the use of function calls

"""

import sys
sys.path.append('C:/Users/sdb58/Documents/python/scratch')

# load the necessary packages
import numpy as np
import matplotlib.pyplot as plt
import sdb_functions as sdb
import myConstants as const

import time

from scipy.interpolate import interp1d
from fast_histogram import histogram1d

plt.close('all')

# experimental parameters
r0 = 3e-4       #the initial rms size of the plasma
Te = 100 * const.c * 100 * const.h / const.kb * 2 / 3
T1 = 2
m1 = 40 * const.mp
vexp = np.sqrt(const.kb * Te / m1) 
vth1 = np.sqrt(const.kb * T1 / m1)
tau_exp = np.sqrt(m1 * r0**2 / (const.kb*(Te + T1)))

# set some computational quantities
N1 = np.int(5e4 )
nbins = np.int(201)
xmax = 3e-3
tmax = 5e-6
vmax_factor = 5
dt = 0.5e-8
nsteps = np.int(np.floor( tmax / dt ))
frac_cut = 10e-3

# %% and initiate arrays, etc.

# ... like the arrays that evolve in time ...
rrms1 = np.zeros((nsteps,1))
vrms1 = np.zeros((nsteps,1))
Tet   = np.zeros((nsteps,1))
t     = np.zeros((nsteps,1))
vrad  = np.zeros((nsteps,1))
pkdens= np.zeros((nsteps,1))

# ... and the grid in space (same in all dimensions) ...
rmax  = np.sqrt(3) * xmax
xgrid = np.linspace(-xmax, xmax, nbins)
dx    = xgrid[1] - xgrid[0]

# ... and the velocity grid ...
vxmax = vmax_factor / np.sqrt(3) * vexp
vxgrid = np.linspace(-vxmax, vxmax, nbins)
dvx = vxgrid[1] - vxgrid[0]

# ... and the radial spatial grid ...
n1    = np.int( (nbins-1) / 2 )
rgrid = np.sqrt( 3 ) * xgrid[n1:nbins] 
dr    = rgrid[1] - rgrid[0]

# ... and the radial velocity grid
vrmax  = vxmax
vrgrid = np.sqrt(3) * vxgrid[n1:nbins]
dvr = vrgrid[1] - vrgrid[0]

# ... and the initial positions of the particles ...
x1, y1, z1 = np.random.randn(3, N1) * r0
r1 = np.sqrt(x1**2 + y1**2 + z1**2)

# ... and the velocity distribution ...
vx1, vy1, vz1 = np.random.randn(3, N1) * vth1
vr1 = np.sqrt( vx1**2 + vy1**2 + vz1**2 )

# %% Perform a few initial checks

# 1. calculate the density from the data and compare to the "right" answer
d1 = sdb.calc_dens(x1,y1,z1,xmax,dx,nbins)
d0 = N1 * (2*np.pi*r0**2)**(-1.5) * np.exp(-rgrid**2 / 2 / r0**2) 
sdb.densplot(rgrid,d0,d1)

# 2. calculate the velocity distribution function from the data and compare
#    to the analytic function. Note that multiplying this by 4 * pi * vr**2 and 
#    integrating will give N1
hvx1 = histogram1d(vx1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvy1 = histogram1d(vy1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvz1 = histogram1d(vz1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
fv1 = hvx1 * hvy1 * hvz1 * N1**(-2) * dvx**(-3)
fv1 = fv1[n1:nbins]
fv0 = N1 * (2*np.pi*vth1**2)**(-1.5) * np.exp( -vrgrid**2 / 2 / vth1**2)
plt.figure()
plt.plot(vrgrid, fv0, vrgrid, fv1,'bo')

# %% Move some particles
"""
Here are the steps, as recorded in ross_003.m
0. let the particles take an euler step
1. calculate the rms size in rrms1
2. calculate the density, the radial velocity profile, and the dlogndr
3. calculate the electron temperature now
4. calculate the friction (ca/ca, ca/yb, yb/ca, yb/yb)
5. calculate the acceleration by the ambipolar field
6. let the particles change their velocities because of the ambipolar field
7. calculate the rms velocity

"""
t1 = time.time()
for n in range(0,nsteps):

    # 0. let the particles take an euler step
    t[n] = n * dt
    x1 = x1 + vx1 * dt
    y1 = y1 + vy1 * dt
    z1 = z1 + vz1 * dt
    
    # 1. calculate the rms size of the plasma
    rrms1[n] = np.sqrt( np.mean(x1**2) + np.mean(y1**2) + np.mean(z1**2) ) / np.sqrt(3)
    
    # 2a. calculate the density from the data
    d1 = sdb.calc_dens(x1,y1,z1,xmax,dx,nbins)
    pkdens[n] = d1[0]    
    
    # 2b. calculate the radial velocity profile --> not yet
    
    # 2c. calcluate dlogndr
    dlogndr = sdb.calc_dlogndr(rgrid, d1, frac_cut)
    
    # 3. calculate the electron temperature right now
    v1squared = np.mean( vx1**2 + vy1**2 + vz1**2 )
    Te_now = Te + T1 - m1 * v1squared / (3*const.kb ) 
    if Te_now < 0:
        Te_now = Tet[n-1]
    Tet[n] = Te_now
    
    # 4. calculate the friction --> not yet, we need 2b first
    
    # 5. calculate ambipolar field: interpolate to get a1 at each particle position
    a1 = - const.kb * Te_now / m1 * dlogndr
    a1_int = interp1d(rgrid,a1, kind='linear')
    a1ambi = a1_int(r1)
   
    # change the x,y,z velocities of the particles
    vx1 = vx1 + a1ambi * dt * x1 / r1
    vy1 = vy1 + a1ambi * dt * y1 / r1
    vz1 = vz1 + a1ambi * dt * z1 / r1
    
    # calculate the new rms velocity
    vrms1[n] = np.sqrt( np.mean(vx1**2 + vy1**2 + vz1**2) * 0.3333 )
    vrad[n]  = np.mean(np.sqrt(((vx1*x1)**2 + (vy1*y1)**2 + (vz1*z1)**2)/(x1**2 + y1**2 + z1**2)))
    
    # let the user know the calculation is alive
    if np.mod(n,25)==0:
        print(' {0:4f}'.format(n/nsteps))
        
e1 = time.time() - t1

plt.figure()
r0t = r0 * np.sqrt(1 + vexp**2 * t**2 / r0**2)
plt.plot(t,rrms1, label='sim')
plt.plot(t,r0t,   label='theory')
plt.legend(loc='best')

plt.figure()
plt.plot(t,vrms1, label = 'sim')
plt.plot(t, np.sqrt(vth1**2 + vexp**2 * t**2 / (t**2 + r0**2/(vexp**2 + vth1**2))), label = 'theory')
plt.xlabel('time')
plt.ylabel('vrms')
plt.legend(loc='best')

plt.figure()
plt.plot(t,Tet, label = 'sim')
plt.plot(t, Te/(1+t**2/tau_exp**2), label = 'theory')
plt.xlabel('time')
plt.ylabel('Te(t) (K)')
plt.legend(loc='best')

plt.figure()
plt.plot(t,vrad, label = 'sim')
plt.ylabel('radial velocity (m/s)')
plt.legend(loc='best')

plt.figure()
plt.plot(t,pkdens, label='sim')
plt.plot(t,np.max(pkdens) / (1+(vexp*t/r0)**2)**(1.5), label = 'theory')
plt.legend(loc='best')

print('elapsed time = {0:4e}'.format(e1))
