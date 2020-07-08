
import pynbody
import pylab
import numpy as np
import matplotlib.pylab as plt
import astropy.units as u


#Now I need a code that will load the simulation (s will stand for simulation)
s = pynbody.load('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096')
h=s.halos()
#The following code will change the units to make it more appealing
s.physical_units()
def findBH(s):
    BH = s.stars[pynbody.filt.LowPass('tform', 0.0)]
    return BH

pynbody.analysis.angmom.faceon(h[5])
BH = findBH(s)
BH_pos = BH['pos']
BHx = BH_pos[:,0]
BHy = BH_pos[:,1]
BHz = BH_pos[:,2]
BH_position = np.array([BHx[0], BHy[0], BHz[0]])

radius= 0.5 #kpc
sphere = pynbody.filt.Sphere(radius, cen =(BH_position))
#sphere=pynbody.analysis.halo.shrink_sphere_center(s, r=None, shrink_factor=0.7, min_particles=100)
stars = s.stars[0:]
in_sphere = stars[sphere]
total_stars = len(in_sphere)
print("Total stars: ",total_stars)




