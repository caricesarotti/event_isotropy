########################################
# emdVar.py
# Code to calcualte the EMD from POT
# written by Cari Cesarotti
# Needs ot library
#######################################
import sys
import time
import warnings
import numpy as np
from numpy import linalg as LA
import random
import ot
from ot.lp import emd2

#########################################                                                                                                                                         
# Functions to calculate some variables
   
# Calculate phi from 3 mom
def phi(vec): # Should only return values between 0 and 2 pi                                                                                                                    
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

# Calculate eta from 3 mom                                                                                                                                                           
def eta(vec):
    if len(vec) !=3:
        raise Exception('emdVar Error: invalid eta argument')
    px, py, pz = vec[0], vec[1], vec[2]
    etaV = np.arctanh(pz/np.sqrt(px**2+py**2+pz**2))
    return etaV

# Define to handle periodicity of phi, processes array of phi values
def preproc(X):
    return np.array([wrapCheck(x) for x in X])

# Ensures that all phi values are between 0 and 2pi                                                                                                                        
def wrapCheck(x):
    if x<0:
        return x+2*np.pi
    if x>2*np.pi:
        return wrapCheck(x-2*np.pi)
    return x

#############################################################
# Define distance metrics.                                   
                                                                                                              
# Calculates euclidean distance squared where the first column is eta, the second is phi         
def _cdist_phi_y(X,Y, ym):
    # define ym as the maximum rapidity cut on the quasi-isotropic event
    # Make sure the phi values are in range                                                                                                                                          
    phi1 = preproc(X[:,1])
    phi2 = preproc(Y[:,1])
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:]))
    y_d = X[:,0,np.newaxis] - Y[:,0]
    norm = 12.0/(np.pi*np.pi+16*ym*ym)
    dist = norm*(phi_d**2 + y_d**2)
    return dist

# Calculates just phi distance
def _cdist_phi(X,Y):
    phi1 = preproc(X)
    phi2 = preproc(Y)
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:]))
    return (4/np.pi)*phi_d

# Calculates just phi distance                                                                                                                                                                  
def _cdist_phicos(X,Y):
    phi1 = preproc(X)
    phi2 = preproc(Y)
    phi_d =np.pi -np.abs(np.pi-np.abs(phi1[:,np.newaxis] - phi2[:]))
    return (np.pi/(np.pi-2))*(1-np.cos(phi_d))

# Calculates the distance on the sphere                                         
def _cdist_sphere(X,Y):
    theta_d=np.array([[np.arccos(np.around(np.dot(X[i],Y[j]), decimals=5)/np.around(LA.norm(Y[j])*LA.norm(X[i]),decimals=5)) for j in range(len(Y))] for i in range(len(X))])
    return theta_d

# We include a scaling factor of pi/2 such that the range is also between 0 to pi
def _cdist_cos(X,Y):
    cos_d=np.array([[2*(1-np.around(np.dot(X[i],Y[j]), decimals=5)/np.around(LA.norm(Y[j])*LA.norm(X[i]),decimals=5)) for j in range(len(Y))] for i in range(len(X))])
    return cos_d
                     
######################################
# EMD CALCULATION                                                                                                                                                                    

## ev0 is the array of values from one event (e.g. pT or energy), ev1 is the array of the values from                                                                                
## the other event. M is the distance matrix between the two events.                                                                                                                 
def emd_Calc(ev0,ev1,M,maxIter=1000000):
    # NORMALIZE IF NOT NORMALIZED                                                                                                                                                    
    ev0norm = ev0[:]/ev0[:].sum()
    ev1norm = ev1[:]/ev1[:].sum()

    cost, log = emd2(ev0norm, ev1norm, M, numItermax=100000000,log=True)
    if cost==0:
        print(log['warning'])
    # returns the normalized EMD between events (e.g. multiply by event pT, eng, etc. to get dimensional value)                                                                      
    return cost
