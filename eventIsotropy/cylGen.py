#
# Created by C. Cesarotti (ccesarotti@g.harvard.edu) 05/2019
#
# Code to generate cylindrically or ring-like quasi isotropic
# events.
#
import sys
import warnings
import numpy as np
import math
import random

######################################
#
def eng(px,py,pz):
    return np.sqrt(px**2+py**2+pz**2)
#
#
def eta(px, py, pz):
    etaV = np.arctanh(pz/np.sqrt(px**2+py**2+pz**2))
    return etaV
#
#
def pT(px, py):
    return np.sqrt(px**2+py**2)
#
######################################
#
## Return the points of a uniformly grided cylinder in eta - phi space. Phi goes from 0 to 2 pi,                                                    
## Eta goes from -etaMax to etaMax                                                                                                                  
## piSeg is how many slices you want in the phi direction, must be integer
## Returns points on the cylinder with roughly equal gridding on the z axis 
def cylinderGen(piSeg, etaMax):

    flag = False
    # Check that piSeg is a positive integer                                                                                                                                                
    if float(piSeg).is_integer():
        if etaMax>0:
            flag=True
        
    if flag:
        # First, calculate the fraction of the points that is along the phi direction                                                                                                                       
        etaSeg = int(math.floor(etaMax*piSeg/np.pi))
        phiVals = [2*np.pi*(j+0.5)/piSeg for j in range(piSeg)]
        etaVals = [-1.0*etaMax + 2.0*etaMax*(j+0.5)/(etaSeg) for j in range(etaSeg)]

        # Build set of points in (phi, eta)                                                                                                                                             
        cylPoints=[]
        for j in range(len(phiVals)):
            points = [[etaVals[i], phiVals[j]] for i in range(len(etaVals))]
            cylPoints = cylPoints+points
        
        # RETURNS 
        # Array of points of the cylinder configuration in (phi, eta) space. 
        return np.array(cylPoints)

    else:
        raise Exception('Error: first argument must be a positive integer, second argument must be positive')
###############################
# Randomly offsets angular position
def cylinderGenShift(piSeg, etaMax):
    ## Return the points of a uniformly grided cylinder in eta - phi space. Phi goes from 0 to 2 pi,                                                              
    ## Eta goes from -etaMax to etaMax                                                                                                                            
    ## piSeg is how many slices you want in the phi direction.                                                                                                                
    ## Returns points on the cylinder with roughly equal gridding on the z axis                                                                     

    flag = False
    # Check that piSeg is a positive integer                                                                                            
                                                                                                                                        
    if float(piSeg).is_integer():
        if etaMax>0:
            flag=True

    if flag:
        ###
        randShift = random.uniform(0,2*np.pi/piSeg)
        etaSeg = int(math.floor(etaMax*piSeg/np.pi))
        phiVals = [2*np.pi*(j)/piSeg+randShift for j in range(piSeg)]
        etaVals = [-1.0*etaMax + 2.0*etaMax*(j+0.5)/(etaSeg) for j in range(etaSeg)]

        # Build set of points in (phi, eta)                                                                                                                        
        cylPoints=[]
        for j in range(len(phiVals)):
            points = [[etaVals[i], phiVals[j]] for i in range(len(etaVals))]
            cylPoints = cylPoints+points
        # RETURNS                                                                                                                                                                     # Array of points of the cylinder configuration in (phi, eta) space.                                                                                                   
        return np.array(cylPoints)


#################################
#
## Returns an array of the pT values for all the particles in the passed array                                                                                      
#
def pTFromVec(vecArray):
    if len(vecArray.shape) != 2 or vecArray.shape[1] != 3:
        raise Exception('cylGen Error: invalid format. Enter array of 3 vectors')
    pTvals = [pT(vec[0], vec[1]) for vec in vecArray]
    return np.array(pTvals)

## Returns an array of the eta values for all the particles in the passed array                                                                                                      
#
def etaFromVec(vecArray):
    if len(vecArray.shape) != 2 or vecArray.shape[1] != 3:
        raise Exception('cylGen Error: invalid format. Enter array of 3 vectors')
    etaVals = [eta(vec[0], vec[1], vec[2]) for vec in vecArray]
    return np.array(etaVals)
#

#### Defines just a ring of particles
def ringGen(piSeg):
   ## Return the points of a uniformly grided ring
    flag = False
    # Check that piSeg is a positive integer                                                                                                                                                    
    if float(piSeg).is_integer():
        flag=True

    if flag:
        # First, calculate the fraction of the points that is along the phi direction                                                                                                                                                                                                                                                  
        phiVals = [2*np.pi*(j+0.5)/piSeg for j in range(piSeg)]

        return np.array(phiVals)

    else:
        raise Exception('Error: first argument must be a positive integer')
#################################                  

#### Defines just a ring of particles                                                                                                                                          
## piSeg is an integer
def ringGenShift(piSeg):
   ## Return the points of a uniformly grided ring                                                                                                                             
    flag = False
    # Check that piSeg is a positive integer                   
    if float(piSeg).is_integer():
        flag=True

    if flag:
        # First, calculate the fraction of the points that is along the phi direction                                            
        randShift = random.uniform(0,2*np.pi/piSeg)
        #randShift = 0. # Don't need random shift for collider events, already random. Just for testing. 
        phiVals = [2*np.pi*j/piSeg+randShift for j in range(piSeg)]
        return np.array(phiVals)

    else:
        raise Exception('Error: first argument must be a positive integer')
#################################                               
