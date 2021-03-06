'''This is Katy's code (https://github.com/04301998/Lupe_project-/blob/master/ACC.py) to plot the black hole luminosity.
         I made some adjustements and added helpful comments '''

import pynbody
import pylab
import numpy as np
import matplotlib.pylab as plt
import readcol
import itertools as it
from itertools import tee
import pandas as pd
import warnings
import decimal
import statistics

# Loading files
files = readcol.readcol('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.orbit')
BHID= 89425759
i =np.where(files[:,0]==  BHID)
print (files[:,0][i])

#i= np.where(BHID== 89425759)
# Convertions:
# The following numbers are from the simulation
m_sol= 2.31e15
l_kpc  = 25000
timee = 38.78  # Time conversion: simulation units time to Gyr 
d_timee = 1.22386438e18 # Time conversion: simulation units to seconds
#t_square= d_time ^2
t_square = 1.49784401e36
m_g = 1.989e33 # Sun mass in gram
l_cm = 3.086e21 # kpc to cm

         
# delta energy
Denergy=( files[:,13][i]* m_sol*( l_kpc**2) *m_g *(l_cm**2))/t_square #units here are ergs
#delta time
Dtime = files[:,14][i]*d_timee
dEdt = Denergy/Dtime
Time=((files[:,1][i]))*timee
#print(Time)
#print(timee)


# Functions:

''' Create 2 parallel iterators (a,b) pointing to the first element of the original iterable.
The second iterator, b is moved 1 step forward (the next(b, None)) call).  a points to c0 and b points to c1.
Both a and b can traverse the original iterator independently - the izip function takes the two iterators and makes pairs 
of the returned elements, advancing both iterators at the same pace.'''

def pair(iterable):
    "c -> (c0,c1), (c1,c2), (c2, c3), ..." # This function creates ordered pairs 
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def float_range(start, stop, step):
  while start < stop:                 # Float Range function 
    yield float(start)
    start += decimal.Decimal(step)
    
#this is time interval from 0 to 13.8 Gyr with an interval of 0.01 which equals to 10 million years
intervals = pair(float_range(0,13.8,0.01))
# tmin,tmax for each time interval
centers = [(tmin+tmax)/2. for tmin, tmax in intervals]

def combining(Time,dEdt,intervals):          
                            # Calculate median valuea given intervals
    warnings.simplefilter("ignore")
    out = []

    for tmin, tmax in intervals:
        filter = (Time >= tmin) & (Time < tmax)
        out.append(np.mean(dEdt[filter]))
    return np.array(out)


b = len(intervals)
#print(centers)
#print(combining(Time, dEdt, intervals))


plt.title(" $\Delta$E/$\Delta$t vs Time")
plt.plot(centers, combining(Time, dEdt, intervals),'ro', label=('CptMarvel')) # FATSO
#plt.scatter(Time, dEdt)
plt.legend(loc = 'upper right')
plt.xlabel("Time(Gyrs)")
plt.ylabel("$\Delta$E/$\Delta$t(Erg/s)")
plt.yscale('log')
plt.ylim(10e35,10e38)
#plt.xlim(0,3)
#plt.show()
plt.savefig("CptMarvel.png")
