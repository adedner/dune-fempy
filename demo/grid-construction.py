from __future__ import print_function

import numpy

import dune.common
import dune.fem

import dune.create as create

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

alugrid = create.grid("ALUConform", {'vertices':points, 'simplex':triangles}, dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction000")

print("GridFactory 2")
alugrid = create.grid("ALUConform", {'vertices':points, 'simplex':triangles}, dimgrid=2)

print("from file 1")
alugrid = create.grid("ALUConform", "../data/unitcube-2d.dgf", dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction001")

print("from file 2")
alugrid = create.grid("ALUConform", (dune.common.reader.dgf,"../data/unitcube-2d.dgf"), dimgrid=2)
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
alugrid = create.grid("ALUConform", dune.grid.string2dgf(dgf), dimgrid=2)
output = alugrid.vtkWriter()
output.write("grid_construction003")

print("cartsesian domain")
alugrid = create.grid("ALUConform", dune.grid.cartesianDomain([0,0],[1,1],[16,16]), dimgrid=2)
