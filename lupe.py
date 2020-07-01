#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
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
Hfiles = readcol.readcol('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/supersample/highres/cptmarvel.test.orbit')
Ffiles =  readcol.readcol('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/supersample/fatso/cptmarvel.fatso.orbit')

# DATA
# Convertions
m_sol= 2.31e15
print m_sol.units
l_kpc  = 25000
m_g = 1.989e33 # Sun mass in gram
l_cm = 3.086e21 # kpc to cm
timee = 38.78
d_timee = 1.22386438e18
t_square = 1.49784401e36
         
# Highes Accretion
Denergy =( Ffiles[:,13]* m_sol*( l_kpc**2) *m_g *(l_cm**2))/t_square
Dtime = Ffiles[:,14]*d_timee
dEdt = Denergy/Dtime
Time =((Ffiles[:,1])-0.9999999)*timee
print(Time)
print(timee)


# Functions

def pair(iterable):
    "c -> (c0,c1), (c1,c2), (c2, c3), ..." # This function creates ordered pairs 
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def float_range(start, stop, step):
  while start < stop:                 # Float Range function 
    yield float(start)
    start += decimal.Decimal(step)

intervals = pair(float_range(0,2,0.25))
centers = [(tmin+tmax)/2. for tmin, tmax in intervals]

def combining(Time,dEdt,intervals):          
                            # Calculate median valuea given intervals
    warnings.simplefilter("ignore")
    out = []

    for tmin, tmax in intervals:
        mask = (Time >= tmin) & (Time < tmax)
        out.append(np.mean(dEdt[mask]))
    return np.array(out)

b = len(intervals)
print(centers)
print(combining(Time, dEdt, intervals))

plt.title(" $\Delta$E/$\Delta$t vs Time")
plt.plot(centers, combining(Time, dEdt, intervals),'ro', label=('FATSO Simulation')) # FATSO
#plt.scatter(Time, dEdt)
plt.legend(loc = 'upper right')
plt.xlabel("Time(Gyrs)")
plt.ylabel("$\Delta$E/$\Delta$t(Erg/s)")
plt.yscale('log')
plt.ylim(10e35,10e38)
#plt.xlim(0,3)
plt.show()
plt.savefig("A_FATSO.png")
