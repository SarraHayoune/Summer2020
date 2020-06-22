import pynbody 
import numpy as np
from numpy import _NoValue
import pandas as pd
import matplotlib.pylab as plt
import readcol
import BH_functions as BHF


files = readcol.readcol('/media/jillian/cptmarvel/files.list')
files = files[:,0]

# function to find black hole
def findBH(s):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    #BHfilter = np.where((s.stars['iord']==60354630) | ( s.stars['iord']==60352986))
    BH = s.stars[BHfilter]
    return BH

#function to find the halos that the galaxy is in
def findBHhalos(s):
    BH = findBH(s)
    BHhalos = BH['amiga.grp']
    return BHhalos


def getz(s):
    return s.properties['z']
def gettime(s):
    return pynbody.analysis.cosmology.age(s)


for j in range (2,3):
    # loading the snapshotS
    s =pynbody.load('/media/jillian/cptmarvel/'+files[j])
    #turns the numbers into useful units
    s.physical_units()
    #halos= [5]
    #loading the halos aka galaxies
    h = s.halos()
    BH = findBH(s)
    BHhalos = findBHhalos(s)
    currenthalo= np.argsort(BHhalos)

     #print "black holes are one black hole now"
            #put the galaxy you care about in the center of the simulation
    pynbody.analysis.angmom.sideon(h[5])
   
    BHposition=BH['pos']

            #putting the x-values into a column

    BHx1= BHposition[[0],0]
            #print "x postion", BHx
                 
            #putting the y-values into a column
    BHy1= BHposition[[0],1]
            #print "y position", BHy
        #putting the z-values into a column
    BHz1= BHposition[[0],2]
            #print "z position", BHz
   
    BHF.render(s,width= '5 kpc',plot= True, ret_im= True)
    plt.plot(BHx1, BHy1,'+')        
    plt.title(str(getz(s)))
    plt.savefig(filename='h5'+',z='+str(getz(s))+'.png')
    
            #plt.plot(BHx2, BHy2,'+')
   # plt.show()


    
