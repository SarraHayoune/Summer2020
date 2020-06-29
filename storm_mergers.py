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
    BHfilter = np.where((s.stars['iord']=243778457)| ( s.stars['iord']==243771992))

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

#f =  open("cptmarvel.dat", "w+") 

for j in range (35,40):


    
    # loading the snapshotS
    s =pynbody.load('/media/jillian/cptmarvel/'+files[j])
    # convert the units 
    s.physical_units()
    #  load any available halo
    h = s.halos()
    BH = findBH(s)
    print BH
    print len(BH)
    BHhalos = findBHhalos(s)
    #sorting the halos, indexes/indecis are like an exact address
    currenthalo = np.argsort(BHhalos)
    print BHhalos[currenthalo]

    for i in currenthalo:
    
        #which halo are we on?
        currenthalo = BHhalos[i]
        print 'current halo: ', currenthalo
        print i
        
        
        # if there are two black holes
        if len(BH)== 2:
  
            if  BHhalos[0]  ==  BHhalos[1]:
                print "black holes are in the same halo"
                #put the galaxy you care about in the center of the simulation
                pynbody.analysis.angmom.sideon(h[currenthalo])
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
                BHF.render(s,width= '10 kpc',plot= True, ret_im= True)
                #plt.title(str(round(getz(s),2)))
                IDbh = BH['iord'][i]
                plt.plot(BHx1, BHy1,'+')
                plt.plot(BHx2, BHy2,'+')
               # plt.axis('off')
                p='z = '+str(round(getz(s),2))
                leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False,fancybox=True)           
                for text in leg.get_texts():
                    plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                distance1 =((BHx1**2)+(BHy1**2)+(BHz1**2))**(.5)
                print "the distance is:", distance1
                distance2 =((BHx2**2)+(BHy2**2)+(BHz2**2))**(.5)
                print "the distance is:", distance2

            else:
                print "black holes are in different halos"
                #put the galaxy you care about in the center of the simulation
                pynbody.analysis.angmom.sideon(h[currenthalo])
                #this is the position of black hole
                BHposition=BH['pos']
                #putting the x-values into a column
                BHx= BHposition[[i],0]
                # BHx2= BHposition[[1],0]
               # print "x postion", BHx
                #putting the y-values into a column
                BHy= BHposition[[i],1]
                # BHy2= BHposition[[i],1]
                print "y positon", BHy
                #putting the z-values into a column
                BHz= BHposition[[i],2]
                #BHz2= BHposition[[i],2]
               # print "z  positon", BHz                
                # create an image using  the default bands (i, v, u)
                IDbh= BH['iord'][i]
                BHF.render(s,width= '10 kpc',plot=True,ret_im=True)
                #plt.title(str(round(getz(s),2)))
                plt.plot(BHx, BHy,'+')
                plt.plot(BHx, BHy, '+')
                #plt.axis('off')          
                p='z = '+str(round(getz(s),2))
                leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False,fancybox=True)
                for text in leg.get_texts():
                    plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
                #the .5 is the square root , this is the distance formula
                distance =((BHx**2)+(BHy**2)+(BHz**2))**(.5)
                print "the distance is:", distance
                
        else:
              
            print " one black hole "
            #put the galaxy you care about in the center of the simulation
            pynbody.analysis.angmom.sideon(h[currenthalo])
            BHposition=BH['pos']
            #putting the x-values into a column
            BHx1= BHposition[[0],0]
            #print "x postion", BHx1                 
            #putting the y-values into a column
            BHy1= BHposition[[0],1]
            #print "y position", BHy1
            #putting the z-values into a column
            BHz1= BHposition[[0],2]
            #print "z position", BHz1            
            IDbh= BH['iord'][i]
            BHF.render(s,width= '10 kpc',plot= True, ret_im= True) 
            plt.plot(BHx1, BHy1,'r+')
            #plt.axis('off')          
            p='z = '+str(round(getz(s),2))
            leg= plt.legend([p],loc= 'upper left',fontsize='medium', frameon= False)           
            for text in leg.get_texts():
               plt.setp(text, color = 'w') #filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')           
            plt.savefig(filename='ID= '+str(IDbh)+',z= '+str(round(getz(s),2))+'.png')
            distance =((BHx1**2)+(BHy1**2)+(BHz1**2))**(.5)
            print "the distance is:", distance
            
                   
        
        #starmass = h[currenthalo].s['mass'].sum()
        #gasmass = h[currenthalo].g['mass'].sum()
        #virialmass = starmass+gasmass+h[currenthalo].d['mass'].sum()
        
        #data = [currenthalo, BH['iord'][i], gettime(s),getz(s), BH['mass'][i], BH['r'][i], starmass, gasmass, virialmass] 
        
        #data= str(data)
        #data= data[1:-1]
        #f.write(data+'\n')
    
        #print data

        
    
#f.close()





