## Created by C. Cesarotti (ccesarotti@g.harvard.edu) 04/2019                                                                                                   
## Last updated: 04/24/20                                            
## 
## Calculates the event isotropy for ring configurations with
## random orientations.
##
#############################
#
import sys
import time
import warnings
import numpy as np
import matplotlib.pylab as plt
import random

from cylGen import ringGen, ringGenShift
from emdVar import _cdist_phicos, emd_Calc

from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
rc('text', usetex=True)
from prettytable import PrettyTable

##############################

# First, generate rings and pT lists
nList=[4,8,16,32,64]

ringSample = np.array([ringGen(nList[i]) for i in range(5)]) # THESE ARE THE PHI VALUES OF THE PARTICLES
ringPtSample=np.array([np.full(len(ringSample[i]), 1.) for i in range(5)])  # THE UNORMALIZED WEIGHT: ALL OF EQUAL PT. NORMALIZATION IN EMD CALC

for i in range(5):
    ringPoints1 = ringSample[i]
    ringPT1 = ringPtSample[i]
    for j in range(5):
        emdSpec=[]
        # SET THE SECOND EVENT WITH j 
        ringPT2 = ringPtSample[j]
        for num in range(1000):
            ringPoints2 = ringGenShift(nList[j]) # The shift just randomly orients the ring, doesn't change particle spacing
            M = _cdist_phicos(ringPoints1,ringPoints2)  # M is distance in phi according to 1 - cos phi metric
            emdval = emd_Calc(ringPT1, ringPT2, M)
            emdSpec.append(emdval)
        f= open("emdRingtoRing"+str(i)+"_"+str(j)+".dat","w+")
        for emdVal in emdSpec:
            f.write(str(emdVal)+ ' ')
        f.close()






