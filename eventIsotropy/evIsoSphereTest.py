##
## Created by C. Cesarotti (ccesarotti@g.harvard.edu) 04/2019                                                                                                         
## Last updated: 04/24/20                                                                                                                                                                            
##                                          
## This file allows users to upload lists of 4-momenta and calculate its EMD from spherical geometry
## Event must be formatted as <event> .... </event>
## Particle information per line is E px py pz
##
## Calculates spherical event isotropy of single event
######################
import sys
import time
import warnings
import numpy as np
import matplotlib.pylab as plt

from spherGen import sphericalGen, engFromVec
from emdVar import _cdist_cos, emd_Calc

from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
rc('text', usetex=True)
from prettytable import PrettyTable

############################################
## Specify input file
if len(sys.argv)<2:
    print 'Error: user did not specify input file'
    sys.exit(2)

## Generate spherical sample
sphereSample = np.array([sphericalGen(i) for i in range(5)])
sphereEng = np.array([engFromVec(sphereSample[j]) for j in range(5)])

## Choose sphere n points
spherePoints1 = sphereSample[2]
sphereEng1 = sphereEng[2]

fileName=sys.argv[1]

momenta=[]
engL=[]
file = open(fileName, 'r') 

if file.closed:
    print('Error: could not open file')
    sys.exit(2)

nextline=file.readline()
while nextline[0:7]=="<event>":
    nextline=file.readline()
    while nextline[0:8]!="</event>":
        particle = [float(n) for n in nextline.split()]
            eng, px, py, pz = particle[0], particle[1], particle[2], particle[3]   
            if eng > 1e-05:
                momenta.append(np.array([px, py, pz]))
                engL.append(eng)
        nextline = file.readline()
    nextline=file.readline()
        
file.close()

## Calculate the \semd values 
M = _cdist_cos(spherePoints1,np.array(momenta)) # Calculates distance with 1 - \cos metric
emdval = emd_Calc(sphereEng1, np.array(engL), M) # Computes EMD
print(emdval)



