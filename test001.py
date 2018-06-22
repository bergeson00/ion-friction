# -*- coding: utf-8 -*-
"""
test001.py

Created by SDB on 6/22/2018

This script loads a 20E6 ions into an x,y,z Cartesian grid. It calculates
the radial position of the particles and then calculates the density. 

This script compares two methods for calculating the density. The first
uses a function called histogram1d from the package fast_histogram. The
second uses the numpy function bincount. These outputs are then compared
to the analytical density. For reasons that are not entirely clear, probably
something to do with the binning of the two programs, the bins need to be
shifted. histogram1d systematically underestimates the "true" density, 
although it is possible that I am making a mistake in the binning. Note that 
the shift is different if you use a grid that is 0,1,2,3.. instead of
0.5 1.5 2.5 ...

It takes about 3 seconds to run

"""

import numpy as np
import matplotlib.pyplot as plt
import time
from fast_histogram import histogram1d

# initiate particles of type 1
N1 = 2e6 
xmax = 1
nbins = 101
ymax = xmax
zmax = xmax
rmax = np.sqrt(xmax**2 + + ymax**2 + zmax**2)

# the grid
dr = rmax / (nbins-1)
r = np.arange(0, rmax+dr, dr)
r = np.arange(0.5*dr, rmax+1.5*dr, dr)

# the location of the particles of type 1
t = time.time()
x1 = np.random.randn(np.int(N1))
y1 = np.random.randn(np.int(N1))
z1 = np.random.randn(np.int(N1))
r1 = np.sqrt(x1**2 + y1**2 + z1**2)
e0 = time.time() - t

# calculate the density
dens0 = N1 / ( 2 * np.pi )**1.5 * np.exp(-r**2 / 2 )
bvol = 4*np.pi * r**2 * dr + np.pi * dr**3 / 3

t = time.time()
dum = histogram1d(r1, range=[0, rmax], bins = nbins)
e1 = time.time()-t

r1_int   =  np.round( (r1-0.5*dr)  / dr ) 
r1_int = r1_int.astype(int)

r_int = r / dr
r_int = r1_int.astype(int)

t = time.time()
bcount = np.bincount(r1_int)
e2 = time.time()-t

n1 = np.size(bcount)
n2 = np.min((n1,nbins))

plt.close()
plt.semilogy(r,dum/ bvol, 'bo',r[0:n2],bcount[0:n2]/bvol[0:n2],'r.' , r , dens0)

print(" histogram1d: {0:4e}, bcount: {1:4e}, initialize arrays: {2:4e}".format(e1,e2,e0))
