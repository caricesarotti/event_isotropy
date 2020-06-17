## Created by C. Cesarotti (ccesarotti@g.harvard.edu) 06/2019
## Last updated: 04/24/20
##
##
## Calculates the event isotropy of a (trasnverse) dijet configuration
## Should be close to zero for n large
##
## Note that event isotropy of a dijet with a cylinder occurs when the
## phi values are opposite but the y values are equal
#############################
import sys
import time
import warnings
import numpy as np
import matplotlib.pylab as plt
import random

from eventIsotropy.cylGen import cylinderGen
from eventIsotropy.emdVar import _cdist_phi_y, emd_Calc

from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
rc('text', usetex=True)
from prettytable import PrettyTable

##############################
# Generate cylinders and pT lists

yMax=2

nList=[4,8,16,32,64]

cylSample = np.array([cylinderGen(nList[i],yMax) for i in range(5)])
cylPtSample=np.array([np.full(len(cylSample[i]), 1.) for i in range(5)])

numPart = np.array([len(cylPtSample[i]) for i in range(5)])

## Calculate 2 particle perfect jet samples in random directions
 # Number of events is 1000
cylPoints1 = cylSample[2]
cylPT1 = cylPtSample[2]

NPENCIL = 100000
pencilPoint= []
for num in range(NPENCIL):
    yVal = yMax
    phiVal = np.pi*random.random()
    penEvent=np.array([[yVal,phiVal],[yVal,np.pi+phiVal]])
    pencilPoint.append(penEvent)
numPencil=2
pencilPt = np.full(numPencil, np.float(1./numPencil))


for i in range(5):
    # SET THE FIRST EVENT WITH i
    cylPoints1 = cylSample[i]
    cylPT1 = cylPtSample[i]
    # CALC EMD FOR ALL PENCIL-LIKE
    emdSpec=[]
    for event in pencilPoint:
        M = _cdist_phi_y(cylPoints1,event, yMax)
        emdval = emd_Calc(cylPT1, pencilPt, M)
        emdSpec.append(emdval)
    filename="emdSpec"+str(len(cylPoints1))+"_CylJetMax.dat"
    f= open(filename,"w+")
    for emdVal in emdSpec:
        f.write(str(emdVal)+ ' ')
    f.close()
