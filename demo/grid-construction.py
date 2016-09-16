import numpy
from mpi4py import MPI
import dune.common
import dune.grid
import dune.fem

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

alugrid = dune.grid.create("ALUConform", {'vertex':points, 'simplex':triangles}, dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction000")

print("GridFactory 2")
alugrid = dune.grid.create("ALUConform", {'vertex':points, 'simplex':triangles}, dimgrid=2)

print("from file 1")
alugrid = dune.grid.create("ALUConform", "../data/unitcube-2d.dgf", dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction001")

print("from file 2")
alugrid = dune.grid.create("ALUConform", (dune.fem.reader.dgf,"../data/unitcube-2d.dgf"), dimgrid=2)
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
alugrid = dune.grid.create("ALUConform", dune.fem.string2dgf(dgf), dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction003")

print("cartsesian domain")
alugrid = dune.grid.create("ALUConform", dune.fem.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
