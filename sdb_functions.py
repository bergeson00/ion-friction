# -*- coding: utf-8 -*-
"""
Created:    Wed Jun 27 08:50:59 2018
Author:     sdb58

Description: 
"""

from fast_histogram import histogram1d
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def calc_dens(x1,y1,z1,xmax,dx,nbins):
    """
    This function calculates the radial density of N1 number of particles
    from three cartesian arrays, x1, y1, z1. These arrays:
        1. Have the same length, N1
        2. Record the x, y, and z positions of the N1 particles.
        
    The function also requires a defined grid. The grid is expected to be
    zero-based, and contain nbins grid points. For example, the grid could
    run from -3 to 3 and contain the points -3, -2, -1, 0, 1, 2, 3. From
    the grid, this function requires the inputs xmax, dx, and nbins.
    
    This function returns the radial density in particles / m**3 at each 
    grid point. The density is defined on a new radial grid where the 
    grid spacing is sqrt(3) * the cartesian grid spacing.
    
    If the x1,y1,z1 arrays contain N1 pseudo-random normally-distributed
    points with an rms width (1 standard deviation) of r0, and if the
    radial grid is 
        rgrid = np.sqrt( 3 ) * xgrid[n1:nbins] 
    Then the theoretical density is equal to 
        d0 = N1 * (2*np.pi*r0**2)**(-1.5) * np.exp(-rgrid**2 / 2 / r0**2) 
    """
    N1 = np.size(x1)
    n1 = np.int( ( nbins - 1 ) / 2 )
    hx1 = histogram1d(x1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
    hy1 = histogram1d(y1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
    hz1 = histogram1d(z1, range=[-xmax - 0.5*dx, xmax + 0.5*dx], bins = nbins)
    d1 = hx1 * hy1 * hz1 * N1**(-2) * dx**(-3)
    d1 = d1[n1:nbins]
    return d1

def densplot(x,y1,y2):
    """
    This is a specialized density difference plot. 
        x  = radial grid
        y1 = analytic density
        y2 = calculated density from the data
    The top plot is the density.
    The bottom plot is the fractional difference in the density
    """
    fig = plt.figure()

    # make a big subplot and turn 
    ax = fig.add_subplot(111)

    ax1 = fig.add_subplot(211)
    ax1.plot(x, y1, label='analytic')
    ax1.plot(x, y2, 'r.', label='data')
    ax1.set_ylabel('density /m**3')

    ax2 = fig.add_subplot(212)
    ax2.plot(x,(y1-y2) / y1)
    ax2.set_ylim(-0.1, 0.1)
    ax2.set_ylabel('(analytic - data) / analytic')

    # Turn off axis lines and ticks of the big subplot and add a common x-label
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_xlabel('radius (m)')

def calc_dlogndr(rgrid,d1,frac_cut):
    """
    Calculate and smooth (dn/dr) / n. Here are some important things to remember
        1.  Below a certain density value, the density can be zero. So to avoid 
            a divide-by-zero problem, we can cut off the density at a certain
            fraction of the max. Typically we'll use frac_cut = 5e-3.
        2.  When the density falls below this value, we back up dn spaces and
            fit those 5 points to a line and then fill in the tail
        3.  Smooth the data
    """
    dn = 5
    dr = rgrid[1] - rgrid[0]
    dlogndr = np.zeros(np.shape(rgrid))
    b = np.argmin( (d1 - frac_cut*np.max(d1))**2 )
    dlogndr[0:b] = np.gradient(d1[0:b]) / d1[0:b]
    coef = np.polyfit( rgrid[b-dn:b], dlogndr[b-dn:b], 1 )
    dlogndr[b-dn:np.size(d1)] = np.polyval(coef, rgrid[b-dn:np.size(d1)])
    dlogndr = dlogndr / dr
    dlogndr[0] = 0
    dlogndr_s = savgol_filter(dlogndr, 11, 3)
    return dlogndr_s