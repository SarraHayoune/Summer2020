import pynbody 
import numpy as np
from numpy import _NoValue
import pandas as pd
import matplotlib.pylab as plt
import readcol
import BH_functions as BHF


#files = readcol.readcol('/media/jillian/cptmarvel/files.list')
#files = files[:,0]
s =pynbody.load('/media/jillian/cptmarvel/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096')
# function to find black hole
def findBH(s):
    BHfilter = pynbody.filt.LowPass('tform',0.0)
    BH = s.stars[BHfilter]
    return BH

#function to find the halos that the galaxy is in
def findBHhalos(s):
    BH = findBH(s)
    BHhalos = BH['amiga.grp']
    return BHhalos

#for file in files:
    
    # loading the snapshotS
    

h = s.halos()

pynbody.analysis.halo.center(h[5],mode='hyb')

BHhalos = findBHhalos(s)
s.physical_units()
BH = findBH(s)
currenthalo = np.argsort(BHhalos)
print BHhalos[currenthalo]
    
    

    
ID =[currenthalo,BH['iord']]

print ID 
