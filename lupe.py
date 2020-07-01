#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import pynbody
import pylab
import numpy as np
import matplotlib.pylab as plt
import readcol
'''
r = (G * BH_Mass) / (stars_vel**2)
G = 6.674e-11
'''

#Now I need a code that will load the snapshots(s will stand for )
Path = "/media/jillian//cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/supersample/highres/"
files = readcol.readcol('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/supersample/highres/files.list')
all_files = files[:,0]

#Tell where BH is Function
def findBH(s):
    BH = s.stars[pynbody.filt.LowPass('tform', 0.0)]
    return BH

#std stands for standard deviation (velocity dispersion)
def DispersionVelocity(s):
    velocity = s.stars['vel']
    x = np.std(velocity[0])
    y = np.std(velocity[1])
    z = np.std(velocity[2])
    #print(x_direct)
    #print(y_direct)
    #print(z_direct)
    dispersion_velocity = np.sqrt( (x)**2 + (y)**2 + (z)**2)
    #print("Dispersion velocity: ",dispersion_velocity)
    return dispersion_velocity

#Need my units to match up so the calculations go correctly
#Couldn't find a way to convert G, so i converted everything else and then converte back to KPC
def RadInfluence(s):
    G = 6.674e-11
    #G is in m**3 * kg**-1 * s**-2
    BH = findBH(s)
    BH_Mass = BH['mass'].in_units('kg')
    #Kg mtches kg in G
    stars_vel = DispersionVelocity(s) * 1e3
    r = (G * BH_Mass) / (stars_vel**2)
    return r * 3.24e-20*3

#Finally converted back to KPC (the conversion is * 3.24e-20)   
for i in all_files:
    s = pynbody.load(Path + i)
    s.physical_units()
    #Don't forget to center the galaxy with this
    pynbody.analysis.angmom.faceon(s)
    BH = findBH(s)
    BH_pos = BH['pos']
    BHx = BH_pos[:,0]
    BHy = BH_pos[:,1]
    BHz = BH_pos[:,2]
    BH_position = np.array([BHx[0], BHy[0], BHz[0]])
    pos_magnitude = np.sqrt((BHx)**2 + (BHy)**2 + (BHz)**2)
    Mass_Msol = BH['mass']
    MassBH = Mass_Msol[0]
    pos_magni = pos_magnitude[0]
    #print(BH_pos)
    #dispersion = DispersionVelocity(s)
    #print(dispersion)
    #The radius here is an array, we need the center to be an integer
    radius = RadInfluence(s)
    radius_influence = radius[0]
    #print(radius)
    #BH_pos is a three int array so it will be the center
    sphere = pynbody.filt.Sphere(radius_influence, cen = BH_position)
    #print(sphere)
    stars = s.stars[np.where(s.stars["tform"]>0)]
    in_sphere = stars[sphere]
    total_stars = len(in_sphere)
    #print("Total stars: ",total_stars)

    #This find their velocity
    velocity = in_sphere['vel']
    #Now we need to find the velocity of these stars in x,y,z
    x = np.array([vel[0] for vel in velocity])
    y = np.array([vel[1] for vel in velocity])
    z = np.array([vel[2] for vel in velocity])
    #Now we can find the average of these by dividing by the total
    vel_answer = np.sqrt((x)**2 + (y)**2 + (z)**2)
    #Now divide by total number of stars
    vel_stars_sphere= vel_answer.sum() / total_stars
    #print("Velocity of the stars in the sphere: ", vel_stars_sphere)

    # Velocity of the stars with respect of the BH
    stars_BH_X = x - BH['vel'][:,0]
    stars_BH_Y = y - BH['vel'][:,1]
    stars_BH_Z = z - BH['vel'][:,2]
    stars_xyz = np.sqrt((stars_BH_X)**2 + (stars_BH_Y)**2 +(stars_BH_Z)**2)
    stars_magnitude = np.sum(stars_xyz)/total_stars
    #print("Velocity of the stars with respect the BH: ", stars_magnitude)

    # VELOCITY OF THE GALAXY
    galaxy = s['vel']
    Galaxy_vx = galaxy[:,0]
    Galaxy_vy = galaxy[:,1]
    Galaxy_vz = galaxy[:,2]
    mass = s['mass']

    # AVERAGE VELOCITY
    Av_Vx = sum(Galaxy_vx*mass)/(sum(mass))
    #print(Av_Vx)
    Av_Vy = sum(Galaxy_vy*mass)/(sum(mass))
    #print(Av_Vy)
    Av_Vz = sum(Galaxy_vz*mass)/(sum(mass))
    #print(Av_Vz)
    
    # BLACK HOLE VELOCITY
    BH_x = Av_Vx - BH['vel'][:,0]
    #print(BH_x)
    BH_y= Av_Vy - BH['vel'][:,1]
    #print(BH_y)
    BH_z= Av_Vz - BH['vel'][:,2]
    #print(BH_z)

    BH_MAGNITUDE= np.sqrt((BH_x)**2 + (BH_y)**2 + (BH_z)**2)
    BH_VEL = BH_MAGNITUDE[0]
    #print("Velocity of the black hole: ",BH_MAGNITUDE)

pynbody.analysis.angmom.faceon(s)
# create an image using  the default bands (i, v, u)
pynbody.plot.stars.render(s,width= '5 kpc',plot=True,ret_im=True,filename='halo11Faceon.png')

# create an image using  the default bands (i, v, u)
pynbody.plot.stars.render(s,width= '5 kpc',plot=True,ret_im=True,filename='halo11Edgeon.png')
plt.show()
