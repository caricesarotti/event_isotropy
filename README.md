# Event Isotropy
### A Robust Measure of Event Isotropy at Colliders (https://arxiv.org/abs/2004.06125)
### C. Cesarotti and J. Thaler
Repository reated December 2019 by Cari Cesarotti (ccesarotti@g.harvard.edu)
$`a^2+b^2=c^2`$
## Event Isotropy Code
To run this code, one neads the Python Optimal Transport Library (https://pot.readthedocs.io/en/stable/)

 - emdVar.py
This is the file with the functions to calcuate the distances between sets for several distance measures as well as the event isotropy. 

### Geometries
For the spherical case, one can calculate the $`1-\cos\theta`$ measure (_cdist_cos(X,Y)) or the $\sqrt{1 - \cos \theta}$ measure (_cdist_sqrt_cos(X,Y)). To calculate, pass these functions arrays X, Y of the 3 momenta of each set. 

For the cylindrical case, one can calculate the squared Euclidean distance in $y-\phi$ space (_cdist_phi_y(X,Y,ymax)) or the unnormalized Euclidean distance in $y-\phi$ (_cdist_phi_y_sqrt(X,Y)) space. Pass these functions arrays of the position in ($y, \phi$) space and the maximum value of $y$ ymax. 

For the ring case, one can calculate the distance in $\phi$ (_cdist_phi(X,Y)) and $1-\cos\phi$ (_cdist_phicos(X,Y)). Pass the function the arrays of $\phi$ values X, Y.

### Event Isotropy Calculation
To calcuate the event isotropy, use the function emd_Calc(ev0,ev1,M) where ev0, ev1 are the energy weights of the event and the uniform event, and M is the distance matrix between them as computed by one of the previous functions. 
Note, this function will also accept user defined distance matrices of the correct dimension.

- spherGen.py

Generates spherical samples and some related quantites. To generate a spherical quasi-uniform event with $n = 12 \times 2^{2i}$ particles, use sphericalGen(i) for $i\in\mathbb{Z}$.
 

## Examples
To run these examples, one neads the Python Optimal Transport Library (https://pot.readthedocs.io/en/stable/) and HEALPix (https://astropy-healpix.readthedocs.io/en/latest/)

emdVar.py

This file includes the calculation of distance between data sets for several ground metrics, as well as the computation of the actual EMD. 



