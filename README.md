# ion-friction
This is the 1D1V PIC code used for ion friction calculation. Creation date: 4/23/2018

This code uses the acceleration derived from the Navier-Stokes equation to calculate
the velocity distribution in a dual-species ultracold neutral plasma. We add an
ion friction term to describe the interaction between the Ca and Yb ions. Initial
versions of the code used a constant term for the Coulomb logarithm. 

1. (6/22/18) I started writing a python calculation. As a first setp I loaded 2e6 particles onto a cartesian grid, calculated their radial position and then calculated the density using two different methods. It looks like the numpy.bincounts function is reasonably fast and closer to the analytic density than using histogram1d from the fast_histogram library. --> test001.py

2. (6/25/18) I changed the density calculation to reduce the noise near the origin. Added thermal velocity calculation. I changed the binning and grid so that histogram1d produces the proper "known" result for density and velocity. --> test002.py

3. (6/27/18) I streamlined the code a little, defining functions in sdb_functions.py and constants in myConstants.py. The file test004.py moves 5e5 particles 1000 steps in 156 seconds. The calculated rms velocity, rms size, electron temperature, and density are reasonably close to the analytic values. However, there are still systematic differences. They appear to be related to how the (Grad n) / n is calculated and depend on smoothing, etc. A reasonable compromise is no smoothing on the density, window = 11, poly = 2 on the golay filtering of dlogndr.

4. (6/28/18) I updated the smoothing on calc_dlogndr in sdb_functions.py. I added in the particles of mass 2 in the main function and ran afew checks to make sure it wasn't broken. --> test005.py. The next step is to think about friction.
