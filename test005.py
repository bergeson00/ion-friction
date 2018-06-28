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
Te  = 100       # the initial electron energy in /cm

# 1. Particles of type 1 -- typically the SMALLER mass particle
m1  = 40 * const.mp
r01 = 3e-4      # the initial rms size of the plasma in meters
n01 = 1.8e10    # the initial peak density of the plasma in /cm^3
T1  = 2         # the initial ion temperature in K

# 2. Particles of type 2 -- typically the LARGER mass particle. Same Te!
m2  = 40 * const.mp
r02 = 3e-4      # the initial rms size of the plasma in meters
n02 = 8.8e10    # the initial peak density of the plasma in /cm^3
T2  = 2         # the initial ion temperature in K

# set some computational quantities
N1 = np.int(1e5 )
N2 = np.int(1e5 )
nbins = np.int(801)
xmax = 6e-3
tmax = 5e-6
vmax_factor = 5
dt = 2e-8
nsteps = np.int(np.floor( tmax / dt ))
frac_cut = 1e-3

# some derived quantities
Te = Te * const.c * 100 * const.h / const.kb * 2 / 3    # convert to K

N1real = n01 * 1e6 * (2 * np.pi)**(1.5) * r01**3        # the "real" number of ions in the plasma
w1 = N1real / N1                                        # the particle weight
vexp1    = np.sqrt(const.kb * Te / m1)                  # asymptotic expansion velocity
vth1     = np.sqrt(const.kb * T1 / m1)                  # thermal velocity
tau_exp1 = np.sqrt(m1 * r01**2 / (const.kb*(Te + T1)))  # characteristic expansion time

N2real = n02 * 1e6 * (2 * np.pi)**(1.5) * r02**3        # the "real" number of ions in the plasma
w2 = N2real / N2                                        # the particle weight
vexp2    = np.sqrt(const.kb * Te / m2)                  # asymptotic expansion velocity
vth2     = np.sqrt(const.kb * T2 / m2)                  # thermal velocity
tau_exp2 = np.sqrt(m2 * r02**2 / (const.kb*(Te + T2)))  # characteristic expansion time

# %% and initiate arrays, etc.

# ... like the arrays that evolve in time ...
Tet     = np.zeros((nsteps,1))
t       = np.zeros((nsteps,1))

pkdens1 = np.zeros((nsteps,1))
rrms1   = np.zeros((nsteps,1))
vrms1   = np.zeros((nsteps,1))
vrad1   = np.zeros((nsteps,1))

pkdens2 = np.zeros((nsteps,1))
rrms2   = np.zeros((nsteps,1))
vrms2   = np.zeros((nsteps,1))
vrad2   = np.zeros((nsteps,1))

# ... and the grid in space (same in all dimensions) ...
rmax  = np.sqrt(3) * xmax
xgrid = np.linspace(-xmax, xmax, nbins)
dx    = xgrid[1] - xgrid[0]

# ... and the velocity grid ...
vxmax = vmax_factor / np.sqrt(3) * vexp1
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
x1, y1, z1 = np.random.randn(3, N1) * r01
r1 = np.sqrt(x1**2 + y1**2 + z1**2)

x2, y2, z2 = np.random.randn(3, N2) * r02
r2 = np.sqrt(x2**2 + y2**2 + z2**2)

# ... and the velocity distribution ...
vx1, vy1, vz1 = np.random.randn(3, N1) * vth1
vr1 = np.sqrt( vx1**2 + vy1**2 + vz1**2 )

vx2, vy2, vz2 = np.random.randn(3, N2) * vth2
vr2 = np.sqrt( vx2**2 + vy2**2 + vz2**2 )

# %% Perform a few initial checks

# 1. calculate the density from the data and compare to the "right" answer
d1 = sdb.calc_dens(x1,y1,z1,xmax,dx,nbins,w1)
d2 = sdb.calc_dens(x2,y2,z2,xmax,dx,nbins,w2)
d01 = N1real * (2*np.pi*r01**2)**(-1.5) * np.exp(-rgrid**2 / 2 / r01**2) 
d02 = N2real * (2*np.pi*r02**2)**(-1.5) * np.exp(-rgrid**2 / 2 / r02**2) 
sdb.densplot(rgrid,d01,d1)
sdb.densplot(rgrid,d02,d2)

# 2. calculate the velocity distribution function from the data and compare
#    to the analytic function. Note that multiplying this by 4 * pi * vr**2 and 
#    integrating will give N1
hvx1 = histogram1d(vx1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvy1 = histogram1d(vy1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvz1 = histogram1d(vz1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
fv1  = hvx1 * hvy1 * hvz1 * N1**(-2) * dvx**(-3) * w1
fv1  = fv1[n1:nbins]
fv01 = N1real * (2*np.pi*vth1**2)**(-1.5) * np.exp( -vrgrid**2 / 2 / vth1**2)

hvx2 = histogram1d(vx2, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvy2 = histogram1d(vy2, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvz2 = histogram1d(vz2, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
fv2  = hvx2 * hvy2 * hvz2 * N2**(-2) * dvx**(-3) * w2
fv2  = fv2[n1:nbins]
fv02 = N2real * (2*np.pi*vth2**2)**(-1.5) * np.exp( -vrgrid**2 / 2 / vth2**2)

plt.figure()
plt.plot(vrgrid, fv01, label='m1, theory')
plt.plot(vrgrid, fv1, 'bo', label='m1, data')
plt.plot(vrgrid, fv02, label='m2, theory')
plt.plot(vrgrid, fv2, 'k.', label='m2, data')
plt.legend(loc='best')

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
    
    x2 = x2 + vx2 * dt
    y2 = y2 + vy2 * dt
    z2 = z2 + vz2 * dt

    # 1. calculate the rms size of the plasma
    rrms1[n] = np.sqrt( np.mean(x1**2) + np.mean(y1**2) + np.mean(z1**2) ) / np.sqrt(3)
    rrms2[n] = np.sqrt( np.mean(x2**2) + np.mean(y2**2) + np.mean(z2**2) ) / np.sqrt(3)
    
    # 2a. calculate the density from the data
    d1 = sdb.calc_dens(x1,y1,z1,xmax,dx,nbins,w1)
    pkdens1[n] = d1[0]    
    
    d2 = sdb.calc_dens(x2,y2,z2,xmax,dx,nbins,w2)
    pkdens2[n] = d2[0]    

    # 2b. calculate the radial velocity profile --> not yet
    
    # 2c. calcluate dlogndr
    dlogndr = sdb.calc_dlogndr(rgrid, d1 + d2, frac_cut)
    
    # 3. calculate the electron temperature right now
    v1squared = np.mean( vx1**2 + vy1**2 + vz1**2 )
    v2squared = np.mean( vx1**2 + vy1**2 + vz1**2 )
    Te_now = ( Te +
                ( N1real * ( T1 - m1 * v1squared / (3*const.kb ) ) + 
                  N2real * ( T2 - m2 * v2squared / (3*const.kb ) ) ) / 
                ( N1real + N2real )
                )
    Tet[n] = Te_now
    
    # 4. calculate the friction --> not yet, we need 2b first
    
    # 5. calculate ambipolar field: interpolate to get a1 at each particle position
    a = - const.kb * Te_now / m1 * dlogndr

    a_int = interp1d(rgrid,a, kind='linear')
    a1ambi = a_int(r1)
    a2ambi = a_int(r2)
   
    # change the x,y,z velocities of the particles
    vx1 = vx1 + a1ambi * dt * x1 / r1
    vy1 = vy1 + a1ambi * dt * y1 / r1
    vz1 = vz1 + a1ambi * dt * z1 / r1

    vx2 = vx2 + a2ambi * dt * x2 / r2
    vy2 = vy2 + a2ambi * dt * y2 / r2
    vz2 = vz2 + a2ambi * dt * z2 / r2
    
    # calculate the new rms and radial velocities
    vrms1[n] = np.sqrt( np.mean(vx1**2 + vy1**2 + vz1**2) * 0.3333 )
    vrad1[n] = np.mean( np.sqrt(((vx1*x1)**2 + (vy1*y1)**2 + (vz1*z1)**2)/(x1**2 + y1**2 + z1**2)))

    vrms2[n] = np.sqrt( np.mean(vx2**2 + vy2**2 + vz2**2) * 0.3333 )
    vrad2[n] = np.mean( np.sqrt(((vx2*x2)**2 + (vy2*y2)**2 + (vz2*z2)**2)/(x2**2 + y2**2 + z2**2)))
    
    # let the user know the calculation is alive, show a few plots if you want
    if np.mod(n,25)==0:
        print(' {0:2f}'.format(n/nsteps))
        #plt.figure()
        #plt.plot(rgrid,dlogndr)
        #plt.plot(rgrid,d1)
        
e1 = time.time() - t1

plt.figure()
r01t = r01 * np.sqrt(1 + vexp1**2 * t**2 / r01**2)
r02t = r02 * np.sqrt(1 + vexp2**2 * t**2 / r02**2)
plt.plot(t * 1e6,rrms1, label='sim1'   )
plt.plot(t * 1e6,r01t,  label='theory1')
plt.plot(t * 1e6,rrms2, label='sim2'   )
plt.plot(t * 1e6,r02t,  label='theory2')
plt.xlabel('time (us)')
plt.legend(loc='best')

plt.figure()
plt.plot(t * 1e6,vrms1, label = 'sim1')
plt.plot(t * 1e6, np.sqrt(vth1**2 + vexp1**2 * t**2 / (t**2 + r01**2/(vexp1**2 + vth1**2))), label = 'theory1')
plt.plot(t * 1e6,vrms2, label = 'sim2')
plt.plot(t * 1e6, np.sqrt(vth2**2 + vexp2**2 * t**2 / (t**2 + r02**2/(vexp2**2 + vth2**2))), label = 'theory2')
plt.xlabel('time (us)')
plt.ylabel('vrms')
plt.legend(loc='best')

plt.figure()
plt.plot(t * 1e6,Tet, label = 'sim')
plt.plot(t * 1e6, Te/(1+t**2/tau_exp1**2), label = 'theory')
plt.xlabel('time (us)')
plt.ylabel('Te(t) (K)')
plt.legend(loc='best')

plt.figure()
plt.plot(t * 1e6,vrad1, label = 'sim1')
plt.plot(t * 1e6,vrad2, label = 'sim2')
plt.xlabel('time (us)')
plt.ylabel('radial velocity (m/s)')
plt.legend(loc='best')

plt.figure()
plt.plot(t * 1e6,pkdens1, label='sim')
plt.plot(t * 1e6,np.max(pkdens1) / (1+(vexp1*t/r01)**2)**(1.5), label = 'theory1')
plt.plot(t * 1e6,pkdens2, label='sim')
plt.plot(t * 1e6,np.max(pkdens2) / (1+(vexp2*t/r02)**2)**(1.5), label = 'theory2')
plt.xlabel('time (us)')
plt.legend(loc='best')

print('elapsed time = {0:3e}'.format(e1))

