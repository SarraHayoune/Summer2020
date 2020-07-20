#this code is to track merger between dwarf galaxies that host supermassive black holes"

import pynbody 
import numpy as np
from numpy import _NoValue
import pandas as pd
import matplotlib.pylab as plt
import readcol
import BH_functions as BHF
from decimal import Decimal

files = readcol.readcol('/media/jillian/cptmarvel/files.list')
files = files[:,0]


# function to find black hole
def findBH(s):
    #BHfilter = pynbody.filt.LowPass('tform',0.0)
    #BHfilter = np.where((s.stars['iord']==243778457))# | ( s.stars['iord']==243771992))#storm
    ##BHfilter = np.where((s.stars['iord']==60353246) | ( s.stars['iord']==60355455))#ruth
    BHfilter = np.where((s.stars['iord']==89425759))
                        #( s.stars['iord']==101863739))#sandra(h148)
    BH = s.stars[BHfilter]
    return BH

#function to find the halos that the galaxy is in
def findBHhalos(s):
    BH = findBH(s)
    BHhalos = BH['amiga.grp']
    return BHhalos

#using the function the halos


def getz(s):
    return s.properties['z']
def gettime(s):
    return pynbody.analysis.cosmology.age(s)
 

for j in range (20,30):
    
    # loading the snapshotS
    s =pynbody.load('/media/jillian/cptmarvel/'+files[j])
    # convert the units 
    s.physical_units()
    #  load any available halo
    h = s.halos()
    BH = findBH(s)
    BHhalos = findBHhalos(s)
    #sorting the halos, indexes/indecis are like an exact address
    currenthalo = np.argsort(BHhalos)
    #print BHhalos[currenthalo]

    for i in currenthalo:
    
        #which halo are we on?
        currenthalo = BHhalos[i]
        print 'current halo: ', currenthalo
        bhid= BH['iord']
        print bhid[i]
        #print i
        
                # if there are two black holes
        if len(BH)== 2:
  
            if  BHhalos[0]  ==  BHhalos[1]:
                print "black holes are in the same halo"
                #put the galaxy you care about in the center of the simulation
                pynbody.analysis.angmom.faceon(h[currenthalo])
                BHposition=BH['pos']
                #putting the x-values into a column
                BHx1= BHposition[[0],0]
                BHx2= BHposition[[1],0]
                #print "x postion", BHx
                #putting the y-values into a column
                BHy1= BHposition[[0],1]
                BHy2= BHposition[[1],1]
               # print "y position", BHy
                #putting the z-values into a column
                BHz1= BHposition[[0],2]
                BHz2= BHposition[[1],2]
                #print "z position", BHz               
                BHF.render(s,width= '25 kpc',plot= True, ret_im= True)
                #plt.title(str(round(getz(s),2)))
                IDbh = BH['iord'][i]
                plt.plot(BHx1, BHy1,'r+',)
                plt.plot(BHx2, BHy2,'c+',)
               # plt.axis('off')
                p='z = '+str(round(getz(s),2))
                leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False)           
                for text in leg.get_texts():
                    plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
        

            else:
                print "black holes are in different halos"
                #put the galaxy you care about in the center of the simulation
                pynbody.analysis.angmom.faceon(h[currenthalo])
                #this is the position of black hole
                BHposition=BH['pos']
                #putting the x-values into a column
                BHx1= BHposition[[0],0]
                BHx2= BHposition[[1],0]
               # BHx3= BHposition[[2],0]
                #BHx4= BHposition[[3],0]
               # print "x postion", BHx
                #putting the y-values into a column
                BHy1= BHposition[[0],1]
                BHy2= BHposition[[1],1]
                #BHy3= BHposition[[2],1]
                #BHy4= BHposition[[3],1]
                #print "y positon", BHy
                #putting the z-values into a column
                BHz1= BHposition[[0],2]
                BHz2= BHposition[[1],2]
                #BHz3= BHposition[[2],2]
                #BHz4= BHposition[[3],2]
               # print "z  positon", BHz                
                # create an image using  the default bands (i, v, u)
                IDbh= BH['iord'][i]
               # IDbh2=BH['iord'][j]
                BHF.render(s,width= '80 kpc',plot=True,ret_im=True)
                #plt.title(str(round(getz(s),2)))
                plt.plot(BHx1, BHy1,'r+')
                plt.plot(BHx2, BHy2, 'c+')
               # plt.plot(BHx3, BHy3,'b+')
               # plt.plot(BHx4, BHy4, 'w+')
                #plt.axis('off')          
                p='z = '+str(round(getz(s),2))
                leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False)
                for text in leg.get_texts():
                    plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                
                plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),3))+'.png')
               
            
        else:
              
            print "one black hole "
            #put the galaxy you care about in the center of the simulation
            pynbody.analysis.angmom.faceon(h[currenthalo])
            BHposition=BH['pos']
            #putting the x-values into a column
            BHx1= BHposition[[i],0]
            #print "x postion", BHx1                 
            #putting the y-values into a column
            BHy1= BHposition[[i],1]
            #print "y position", BHy1
            #putting the z-values into a column
            BHz1= BHposition[[i],2]
            #print "z position", BHz1            
            IDbh= BH['iord'][i]
            BHF.render(s,width= '25 kpc',plot= True, ret_im= True) 
            plt.plot(BHx1, BHy1,'r+')
            #plt.axis('off')          
            p='z = '+str(round(getz(s),2))
            leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False)           
            for text in leg.get_texts():
               plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')           
            plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
           
            
                   
        
       




