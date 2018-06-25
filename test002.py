# -*- coding: utf-8 -*-
"""
test002.py

Created by SDB on 6/25/2018

This script loads a 2E6 ions into an x,y,z Cartesian grid. It calculates
the radial position of the particles and then calculates the density. 

It also loads the ions with a thermal velocity into a Cartesian velocity
space, calculates and plots the velocity distribution function

"""

# load the necessary packages
import numpy as np
import matplotlib.pyplot as plt
from fast_histogram import histogram1d
plt.close('all')

# define some constants
mp = 1.67e-27
kb = 1.38e-23
c = 2.99792458e8
h = 6.63e-34

# experimental parameters
r0 = 3e-4       #the initial rms size of the plasma
Te = 100 * c * 100 * h / kb * 2 / 3
T1 = 2
m1 = 40 * mp
vexp = np.sqrt(kb * Te / m1) 
vth1 = np.sqrt(kb * T1 / m1)

# set some computational quantities
N1 = np.int(5e6 )
nbins = np.int(201)
xmax = 3e-3

# define the grid in x,y,z, (same in all dimensions)
ymax = xmax
zmax = xmax
rmax = np.sqrt(xmax**2 + + ymax**2 + zmax**2)
n1 = np.int( (nbins-1) / 2 )    # we'll need this for the radial grid
xgrid = np.linspace(-xmax, xmax, nbins)
dx = xgrid[1] - xgrid[0]

# fill the arrays wtih random numbers
x1, y1, z1 = np.random.randn(3, N1) * r0

# calculate the density from the data
hx1 = histogram1d(x1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
hy1 = histogram1d(y1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
hz1 = histogram1d(z1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
d1 = hx1 * hy1 * hz1 * N1**(-2) * dx**(-3)
d1 = d1[n1:nbins]

# calculate the known analytic density 
rgrid = np.sqrt( 3 ) * xgrid[n1:nbins] 
dr = rgrid[1] - rgrid[0]
d0 = N1 * (2*np.pi*r0**2)**(-1.5) * np.exp(-rgrid**2 / 2 / r0**2) 

# plot and compare the densities
fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.plot(rgrid, d0, rgrid, d1, 'ko')

ax2 = fig.add_subplot(212)
ax2.plot(rgrid,(d0-d1) / d0)

plt.show()

# ~~~ now the velocity distribution function ~~~ #

vxmax = 5 / np.sqrt(3) * vexp
vymax = vxmax
vzmax = vxmax
vrmax = np.sqrt(vxmax**2 + vymax**2 + vzmax**2)

vxgrid = np.linspace(-vxmax, vxmax, nbins)
dvx = vxgrid[1] - vxgrid[0]

vx1, vy1, vz1 = np.random.randn(3, N1) * vth1

# calculate the velocity distribution function from the data
hvx1 = histogram1d(vx1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvy1 = histogram1d(vy1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
hvz1 = histogram1d(vz1, range=[-vxmax - 0.5*dvx, vxmax + 0.5*dvx], bins = nbins)
fv1 = hvx1 * hvy1 * hvz1 * N1**(-2) * dvx**(-3)
fv1 = fv1[n1:nbins]

vrgrid = np.sqrt( 3 ) * vxgrid[n1:nbins]
dvr = vrgrid[1] - vrgrid[0]
fv0 = N1 * (2*np.pi*vth1**2)**(-1.5) * np.exp( -vrgrid**2 / 2 / vth1**2)

# plot the density
plt.figure()
plt.plot(vrgrid, fv0, vrgrid, fv1,'bo')
