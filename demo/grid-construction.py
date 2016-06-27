import numpy
from mpi4py import MPI
import dune.common
import dune.fem
import dune.fem.grid

from scipy.spatial import Delaunay
# import matplotlib.pyplot as plt

n_angles = 36
n_radii = 8

radii = numpy.linspace(1.0 / n_radii, 1.0, n_radii)
angles = numpy.linspace(0, 2*numpy.pi, n_angles, endpoint=False)
angles = numpy.repeat(angles[..., numpy.newaxis], n_radii, axis=1)

x = numpy.append(0, (radii*numpy.cos(angles)).flatten())
y = numpy.append(0, (radii*numpy.sin(angles)).flatten())

points = numpy.stack((x,y), axis=-1)
triangles = Delaunay(points).simplices

# plt.triplot(points[:,0], points[:,1], triangles.copy())
# plt.plot(points[:,0], points[:,1], 'o')
# plt.show()

m_alugrid = dune.fem.grid.get("ALUSimplexGrid", dimgrid=2, refinement="conforming")

print("GridFactory 1")
alugrid = m_alugrid.LeafGrid(m_alugrid.reader( {'vertex':points, 'simplex':triangles} ))
output = alugrid.vtkWriter()
output.write("grid_construction000")

print("GridFactory 2")
alugrid = dune.fem.grid.leafGrid({'vertex':points, 'simplex':triangles}, "ALUSimplexGrid", dimgrid=2, refinement="conforming")

print("from file 1")
alugrid = dune.fem.grid.leafGrid("../data/unitcube-2d.dgf", "ALUSimplexGrid", dimgrid=2, refinement="conforming")
output = alugrid.vtkWriter()
output.write("grid_construction001")

print("from file 2")
alugrid = dune.fem.leafGrid( (dune.fem.reader.dgf,"../data/unitcube-2d.dgf"), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
output = alugrid.vtkWriter()
output.write("grid_construction002")

print("string 2 dgf")
dgf = """
INTERVAL
0  0
1  1
16 16
#
"""
alugrid = dune.fem.leafGrid( dune.fem.string2dgf(dgf), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
output = alugrid.vtkWriter()
output.write("grid_construction003")

print("cartsesian domain")
alugrid = dune.fem.leafGrid( dune.fem.cartesianDomain([0,0],[1,1],[16,16]), "ALUSimplexGrid", dimgrid=2, refinement="conforming")
