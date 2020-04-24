##
## Created by C. Cesarotti (ccesarotti@g.harvard.edu) 03/2019     
## Code to generate spherical quasi isotropic events with the healpix library
##
##   

import sys
import warnings
import numpy as np
import astropy_healpix
import astropy_healpix.healpy as hp
import random

##################
## Use to generate random directions with normalized vectors                                                                                                                                                                                                  
def sample_spherical(ndim=3):
    vec = np.random.randn(ndim)
    vec /= np.linalg.norm(vec, axis=0)
    return list(vec)

# DEFINED FOR MASSLESS PARTICLES
def eng(px,py,pz):
    return np.sqrt(px**2+py**2+pz**2)

def eta(px, py, pz):
    etaV = np.arctanh(pz/np.sqrt(px**2+py**2+pz**2))
    return etaV

def pT(px, py): 
    return np.sqrt(px**2+py**2)

##################

## Function generates the point on a evenly tiled sphere using HEALPIX from Astropy_HEALPIX
## Returns arrays of the points (3 vectors), the 3 momentum of the particles
def sphericalGen(nVal, etaMax=100):
    
    # Only allowed values of number of points is 12*(2**2i) for i an integer. 
    
    # etaMax must be positive
    if etaMax<0:
        raise Exception('spherGen Error: Invalid eta value.')
    # nVal must be positive integer
    if not (float(nVal).is_integer()) or nVal < 0:
        raise Exception('spherGen Error: Invalid number value')
        
    nside=2**nVal 
    numPix = 12*(nside**2)
    spherPoint = []
    ## USE HEALPIX TO ACCESS POINTS
    for j in range(numPix):
        vecOnSpher = hp.pix2vec(nside, j)
        point = np.array([vecOnSpher[0], vecOnSpher[1], vecOnSpher[2]])
        if abs(eta(vecOnSpher[0], vecOnSpher[1], vecOnSpher[2])) < etaMax:
            spherPoint.append(point)
    return np.array(spherPoint)

def sphericalThetaGen(nVal):
    # Returns only the theta information of the event
    nside=2**nVal
    numPix=12*(nside**2)
    spherPoint=[]
    # Use HEALPix to generate points
    for j in range(numPix):
        theta, phi = hp.pix2ang(nside, j)
        if phi>np.pi:
            point = np.array([np.sin(theta), 0,np.cos(theta)])
            spherPoint.append(point)
        else:
            point = np.array([-np.sin(theta), 0,np.cos(theta)])
            spherPoint.append(point)
    return np.array(spherPoint)
    

## Returns an array of the energy values (assuming massless particles) in the passed array
def engFromVec(vecArray):
    if len(vecArray.shape) != 2 or vecArray.shape[1] != 3:
        raise Exception('spherGen Error: invalid format. Enter array of 3 vectors')
    engVals = [eng(vec[0], vec[1], vec[2]) for vec in vecArray]
    return np.array(engVals)

## Returns an array of the pT values for all the particles in the passed array
def pTFromVec(vecArray):
    if len(vecArray.shape) != 2 or vecArray.shape[1] != 3:
        raise Exception('spherGen Error: invalid format. Enter array of 3 vectors')
    pTvals = [pT(vec[0], vec[1]) for vec in vecArray]
    return np.array(pTvals)

## Returns an array of the eta values for all the particles in the passed array
def etaFromVec(vecArray):
    if len(vecArray.shape) != 2 or vecArray.shape[1] != 3:
        raise Exception('spherGen Error: invalid format. Enter array of 3 vectors')
    etaVals = [eta(vec[0], vec[1], vec[2]) for vec in vecArray]
    return np.array(etaVals)

## Returns a sphere + dijet event                                                                                                                                                                        
def sphereAndDijet(nVal, djFrac):
    # Only allowed values of number of points is 12*(2**2i) for i an integer.                                                                                                                            
    # nVal must be positive integer                                                                                                                                                                      
    if not (float(nVal).is_integer()) or nVal < 0:
        raise Exception('spherGen Error: Invalid number value')

    nside=2**nVal
    numPix = 12*(nside**2)
    spherPoint = []
    ## USE HEALPIX TO ACCESS POINTS                                                                                                                                                                      
    djPoint1 = np.array([1,0,0])
    djPoint2 = np.array([-1,0,0])
    spherPoint.append(djPoint1)
    spherPoint.append(djPoint2)
    vec = sample_spherical()
    ux, uy, uz = vec[0], vec[1], vec[2]
    x=random.uniform(0, 2*np.pi)
    rotMat = np.array([[ux*ux*(1-np.cos(x)) + np.cos(x), ux*uy*(1-np.cos(x))-uz*np.sin(x), ux*uz*(1-np.cos(x))+uy*np.sin(x)], [ux*uy*(1-np.cos(x))+uz*np.sin(x), uy*uy*(1-np.cos(x))+np.cos(x), uy*uz*(1-np.cos(x))-ux*np.sin(x)], [ux*uz*(1-np.cos(x))-uy*np.sin(x), uy*uz*(1-np.cos(x))+ux*np.sin(x), uz*uz*(1-np.cos(x))+np.cos(x)]])
    for j in range(numPix):
        vecOnSpher = hp.pix2vec(nside, j)
        point = np.array([vecOnSpher[0], vecOnSpher[1], vecOnSpher[2]])
        pointRot = rotMat.dot(np.array(point))
        spherPoint.append(pointRot)
    return np.array(spherPoint)

