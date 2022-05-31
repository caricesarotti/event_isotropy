#
# Created by C. Cesarotti (ccesarotti@g.harvard.edu) 04/2019
# Last updated: 04/16/20
# 
# Contains code to calculate distance matrices for spherical, 
# ring like, and cylindrical geometries. 
# Also include code to calculate the EMD between two particle
# physics collider events. 
#
# Dependencies: POT library
#
##########################################
import sys
import time
import warnings
import numpy as np
from numpy import linalg as LA
import random
import ot
from ot.lp import emd2

#########################################                                                                                                                                            
# PROCESSING FUNCTIONS

# Calculate phi from 3 momentum                                                                                                                                                                   
def phi(vec): # Should only return values between 0 and 2 pi                                                                                                                    
#
    if len(vec) !=3:
        raise Exception('emdVar Error: invalid phi argument')
    px, py  = vec[0], vec[1]
    if px==0:
        if py>0:
            return np.pi/2.
        else:
            return 3*np.pi/2.
    phiV = np.arctan(py/px)
    if px<0:
        phiV=np.pi+phiV
    else:
        if py<0:
            phiV=2*np.pi+phiV
    return phiV

# Calculate eta from 3 momentum                                                                                                                                                                     
def eta(vec):
    if len(vec) !=3:
        raise Exception('emdVar Error: invalid eta argument')
    px, py, pz = vec[0], vec[1], vec[2]
    etaV = np.arctanh(pz/np.sqrt(px**2+py**2+pz**2))
    return etaV

# Define to handle periodicity of phi   
# Processes array of phi values                                                                                                                                                      
def preproc(X):
    return np.array([wrapCheck(x) for x in X])

# Ensures that all phi values are between 0 and 2pi                                                                                                                        
def wrapCheck(x):
    if x<0:
        return x+2*np.pi
    if x>2*np.pi:
        return wrapCheck(x-2*np.pi)
    return x


##########################################################################
# Define distance metrics.                                                                                                                                                 


#######################################
## CYLINDRICAL GEOMETRY     
####################################### 

# Calculates euclidean distance where the first column is eta, the second is phi, eg elements in both X and Y are (y,phi)
# This is the distance on the cylinder, beta = 2 measure                                                                                                                           
# ym is max rapidity, needed for correct normalization
def _cdist_phi_y(X,Y, ym):
    # define ym as the maximum rapidity cut on the quasi-isotropic event
    # Make sure the phi values are in range                                                                                                                                          
    phi1 = preproc(X[:,1])
    phi2 = preproc(Y[:,1])
    # Trick to account for phi distance periodicity
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:])) 
    y_d = X[:,0,np.newaxis] - Y[:,0]
    norm = 12.0/(np.pi*np.pi+16*ym*ym)
    dist = norm*(phi_d**2 + y_d**2)
    return dist

# Distance on cylinder, beta = 1 metric 
# first column is eta, the second is phi, eg elements in both X and Y are (y,phi) 
def _cdist_phi_y_sqrt(X,Y):
    # NOTE: THIS IS NOT NORMALIZED!! DOES NOT RUN FROM 0 TO 1
    # Make sure the phi values are in range                                                                   
    phi1 = preproc(X[:,1])
    phi2 = preproc(Y[:,1])
    # Trick to account for phi distance 'wrap around'                                                            
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:]))
    y_d = X[:,0,np.newaxis] - Y[:,0]
    dist = phi_d**2 + y_d**2
    return np.sqrt(dist)

#######################################                                           
## RING LIKE GEOMETRY
####################################### 
# Calculates distance on ring, phi metric
# X, Y are arrays of phi
def _cdist_phi(X,Y):
    phi1 = preproc(X)
    phi2 = preproc(Y)
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:])) # GIVES MATRIX OF DIFFERENCE OF PHI VALUES  
    return (4/np.pi)*phi_d

# Calculates distance on ring, cos phi measure
# X, Y are arrays of phi                                                                                                                                                     
def _cdist_phicos(X,Y):
    phi1 = preproc(X)
    phi2 = preproc(Y)
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:])) # GIVES MATRIX OF DIFFERENCE OF PHI VALUES       
    return (np.pi/(np.pi-2))*(1-np.cos(phi_d))

#######################################
## SPHERICAL GEOMETRY
#######################################
# Calculates the distance on the sphere, angluar distance                                         
# X, Y are arrays of 3 momenta of the particles in the event
def _cdist_sphere(X,Y):
    theta_d=np.array([[np.arccos(np.around(np.dot(X[i],Y[j]), decimals=5)/np.around(LA.norm(Y[j])*LA.norm(X[i]),decimals=5)) for j in range(len(Y))] for i in range(len(X))])
    return theta_d

# Calculates disntace on sphere, cos distance
# X, Y are arrays of 3 momenta of the particles in the event
def _cdist_cos(X,Y):
    cos_d=np.array([[2*(1-np.around(np.dot(X[i],Y[j]), decimals=5)/np.around(LA.norm(Y[j])*LA.norm(X[i]),decimals=5)) for j in range(len(Y))] for i in range(len(X))])
    return cos_d

# Distance on sphere, sqrt cos distance
# X, Y are arrays of 3 momenta of the particles in the event
def _cdist_sqrt_cos(X,Y):
    cos_d=np.array([[(3./2.)*np.sqrt(1-np.around(np.dot(X[i],Y[j]), decimals=5)/np.around(LA.norm(Y[j])*LA.norm(X[i]),decimals=5)) for j in range(len(Y))] for i in range(len(X))])
    return cos_d


                     
######################################
# EMD CALCULATION                                                                                                                                                                    

## ev0 is the ARRAY of eng. measure values from one event (e.g. pT or energy), ev1 is the ARRAY of the values from                                                                                
## the other event. M is the distance MATRIX between the two events. Can be calculated by user with different 
## distance measures or with the distance measures defined above. When calculating event isotropy, one of the 
## events must be the quasi-uniform event.                                                                                                                  

def emd_Calc(ev0,ev1,M,maxIter=1000000):
    # NORMALIZE IF NOT NORMALIZED                                                                                                                                                    
    ev0norm = ev0[:]/ev0[:].sum()
    ev1norm = ev1[:]/ev1[:].sum()

    cost, log = emd2(ev0norm, ev1norm, M, numItermax=100000000,log=True)
    # Should only return 0 when two events are identical. If returning 0 otherwise, problems in config
    if cost==0:
        print(log['warning'])
    # returns the normalized EMD between events (e.g. multiply by event pT, eng, etc. to get dimensional value)                                                                      
    return cost
