# Event Isotropy
### A Robust Measure of Event Isotropy at Colliders (https://arxiv.org/abs/2004.06125)
### C. Cesarotti and J. Thaler
Repository created December 2019 by Cari Cesarotti (ccesarotti@g.harvard.edu)

## Installation

### From this repository

In your Python environment from the top level of this repository run

```
python -m pip install .
```

### From GitHub

In your Python environment run

```
python -m pip install "git+https://github.com/caricesarotti/event_isotropy.git#egg=eventIsotropy"
```

## Event Isotropy Code

### `emdVar.py`
This is the module with the functions to calculate the distances between sets for several distance measures as well as the event isotropy.

- Geometries
For the spherical case, one can calculate the <img src="https://render.githubusercontent.com/render/math?math=1-\cos\theta"> measure (`_cdist_cos(X,Y)`) or the <img src="https://render.githubusercontent.com/render/math?math=\sqrt{1-\cos\theta}"> measure (`_cdist_sqrt_cos(X,Y)`). To calculate, pass these functions arrays X, Y of the 3 momenta of each set.

For the cylindrical case, one can calculate the squared Euclidean distance in <img src="https://render.githubusercontent.com/render/math?math=y-\phi"> space (`_cdist_phi_y(X,Y,ymax)`) or the unnormalized Euclidean distance in <img src="https://render.githubusercontent.com/render/math?math=y-\phi"> (`_cdist_phi_y_sqrt(X,Y)`) space. Pass these functions arrays of the position in <img src="https://render.githubusercontent.com/render/math?math=(y,\phi)"> space and the maximum value of `y` ymax.

For the ring case, one can calculate the distance in <img src="https://render.githubusercontent.com/render/math?math=\phi"> (`_cdist_phi(X,Y)`) and <img src="https://render.githubusercontent.com/render/math?math=1-\cos\phi"> (`_cdist_phicos(X,Y)`). Pass the function the arrays of <img src="https://render.githubusercontent.com/render/math?math=\phi"> values X, Y.

- Event Isotropy Calculation
To calculate the event isotropy, use the function `emd_Calc(ev0,ev1,M)` where ev0, ev1 are the energy weights of the event and the uniform event, and M is the distance matrix between them as computed by one of the previous functions.
Note, this function will also accept user defined distance matrices of the correct dimension.

### `spherGen.py`

Generates spherical samples and some related quantities. To generate a spherical quasi-uniform event with <img src="https://render.githubusercontent.com/render/math?math=n=12\times2^{2i}"> particles for <img src="https://render.githubusercontent.com/render/math?math=i\in\mathbb{Z}">, use `sphericalGen(i)`

### `cylGen.py`

Generates cylindrical samples and ring-like samples, as well as related quatities. To generate a cylinder with uniform tiling in the <img src="https://render.githubusercontent.com/render/math?math=y-\phi"> plane, use `cylinderGen(piSeg, etaMax)` where piSeg is an integer, the number of slices in <img src="https://render.githubusercontent.com/render/math?math=\phi">.

For a ring-like sample, use `ringGen(piSeg)` where piSeg is an integer, the number of slices in <img src="https://render.githubusercontent.com/render/math?math=\phi">.

## Examples

There are three example programs in the `examples` directory.
The examples have external dependencies, but can be installed through the `test` extra with

```
python -m pip install --upgrade -e .[test]
```

### `evIsoSphere.py` file i

Calculates the spherical event isotropy of an uploaded event. User specifies `file` name of event and spherical index `i`. The file should be saved as the four momenta of the particles. `i` is the sphere index, not particle number, and corresponds to the tiling of the quasi uniform comparison sphere. The number of particles in the sphere will be <img src="https://render.githubusercontent.com/render/math?math=n=12\times2^{2i}">.

### `evIsoCyl.py`

Calculates the cylindrical event isotropy of very coarse to very finely gridded cylinders for random extremal transverse dijet configurations. In the limit of large `n`, we expect to see the event isotropy drop to zero. No user input needed.

### `evIsoRing.py`

Calculates the ring isotropy between two rings at random orientations and with potentially different particle number. Generates 1000 random configurations to average event isotropy for given configuration. No user input needed.
