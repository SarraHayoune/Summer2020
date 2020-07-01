#import all the good stuff!
#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import pynbody
import pylab
import numpy as np
import matplotlib.pylab as plt
import astropy.units as u

#Separate stars from everything (Done)
#Find AVG velocity of stars

#Now I need a code that will load the simulation (s will stand for simulation)
s = pynbody.load("cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096.5.std")

#The following code will change the units to make it more appealing
s.physical_units()

#This code creates a sphere around the black holes position
radius = 0.5
sphere= pynbody.filt.Sphere(radius, cen= (-634.00464133,1258.07020815, 29.86851614))
print(sphere)

#This code tells us how many stars are in this section
num_of_stars = s.stars[0:]
in_sphere = num_of_stars[100]
total_stars = len(in_sphere)
#print(total_stars)

#Find the velocities of each star
velocity = in_sphere['vel']
#print(velocity)

#Now we need to find the velocity of these stars in x,y,z 
x = np.array([vel[0] for vel in velocity])
y = np.array([vel[1] for vel in velocity])
z = np.array([vel[2] for vel in velocity])

#Now we can find the average of these by dividing by the total
vel_answer = np.sqrt((x)**2 + (y)**2 + (z)**2)
#print(vel_answer)
print(s['vel'].units)

#Now divide by total number of stars
velocity = vel_answer.sum() / total_stars
print(velocity)

plt.figure(1)
plt.subplot(221)
plt.hist(x, color = 'r', bins = 100)
plt.xlabel("x")
plt.ylabel("frequency")
plt.axvline(x.mean(), color='k', linestyle='dashed', linewidth=1)
plt.legend()
plt.grid()

plt.subplot(222)
plt.hist(y, color = 'b', bins = 100)
plt.xlabel("y")
plt.ylabel("frequency")
plt.axvline(y.mean(), color='k', linestyle='dashed', linewidth=1)
plt.legend()
plt.grid()

plt.subplot(223)
plt.hist(z, color = 'g', bins = 100)
plt.xlabel("z")
plt.ylabel("frequency")
plt.axvline(z.mean(), color='k', linestyle='dashed', linewidth=1)
plt.legend()
plt.grid()
#plt.show()

#This is the function that found the black hole
def findBH(s):
    BH = s.stars[pynbody.filt.LowPass('tform', 0.0)]
    return BH
Found_one= findBH(s)


#Now I just ask the computer to print the mass of the object found in this funct
BH_mass = [Found_one['mass'], Found_one['r']]
print(BH_mass)


'''
x = np.array([vel[0] for vel in velocity])
y = np.array([vel[1] for vel in velocity])
z = np.array([vel[2] for vel in velocity])
i = x / total_stars
j = y / total_stars
k = z / total_stars
#Now we can find the average of these by dividing by the total
vel_answer = np.sqrt((i)**2 + (j)**2 + (k)**2)
print(vel_answer)
print(s['vel'].units)'''
