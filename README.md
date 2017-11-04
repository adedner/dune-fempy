Python Bindings for the DUNE-FEM module
=======================================

dune-fempy aims to provide Python bindings for the [dune-fem][femlink] discretization
module. Here an example for solving Poisson's equation:

```python
import math
from ufl import *
import dune.ufl
import dune.create as create

grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)

# set up a diffusion reaction model using UFL
uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
# provide an exact solution that will be used to add suitable forcing and dirichlet b.c.
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
model = create.model("elliptic", grid, a==0, exact=exact, dirichlet={ 1:exact } )

# set up a space and a conforming finite element scheme and solve the PDE
space  = dune.create.space("Lagrange", grid, dimrange=1, order=1)
scheme = create.scheme("h1", space, model, "cg")
uh,info = scheme.solve()

# make 'exact' into a grid function for output and uh into an UFL coefficient for error computation
exact_gf = create.function("ufl", grid, "exact", 5, exact)
uh_coeff = dune.ufl.GridCoefficient(uh)
# now define a grid function representing the pointwise error
l2error_gf = create.function("ufl", grid, "error", 5, as_vector([(exact[0]-uh_coeff[0])**2]) )

error = math.sqrt( l2error_gf.integrate() )
print("size:",grid.size(0),"L2-error:",error)
grid.writeVTK("laplace", pointdata=[ uh, l2error_gf, exact_gf ])
```

See the file COPYING for full copying permissions.

Dependencies
------------

dune-fempy depends on the following DUNE modules:
- dune-common 2.7+
- dune-geometry 2.7+
- dune-grid 2.7+
- dune-corepy 2.7+
- dune-fem 2.7+

In addition to the dependencies of these DUNE modules, the following software
packages are required:

- a C++14 compatible C++ compiler (e.g., GCC 4.9+, clang 3.5+)
- Python 2.7+

We strongly recommend installing the following Pyhton packages to make full
use of dune-fempy:

- numpy
- mpi4py
- ufl (2016.1.0+)


Quick Dive-In Using Docker
--------------------------

Users who simply want to use the functionality of `dune-fempy` as-is, e.g.,
for experimenting, can do so using Docker and a web browser.

To run the docker container with the latest release just type
```
docker run --rm -v dune:/dune -p 127.0.0.1:8888:8888 registry.dune-project.org/dune-fem/dune-fempy
```
at the command prompt and connect to `localhost:8888` using your favorite web
browser.
Log into Jupyter using the password `dune` and have fun.
The development version can be accessed using the `unstable` tag, i.e.,
```
docker run --rm -v dune:/dune -p 127.0.0.1:8888:8888 registry.dune-project.org/dune-fem/dune-fempy:unstable
```

For your convenience, the demo notebooks are provided in the folder
`dune-fempy`.

With the above command, your notebooks will be kept persistent in a Docker
volume called `dune`.
To remove this volume, simply type
```
docker volume rm dune
```
at the command prompt.


Installation
------------

To download/install `dune-fempy` on your system
follow the instructions in [dune-python][corepy].
After building `dune-fempy` add the `python` folder in the build directory
to your `PYTHON_PATH` or (having set the cmake `PYTHON_INSTALL_LOCATION` flag)
use the `setup-dunepy.py` script or call `make python_install1` in the
build directory of dune-fempy to install the `dune-fem` part of the Dune
python package.

[corepy]: https://gitlab.dune-project.org/staging/dune-python
[femlink]: https://gitlab.dune-project.org/dune-fem/dune-fem
