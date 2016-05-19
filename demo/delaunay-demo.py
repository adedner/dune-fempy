import numpy
import dune.fem.grid

from scipy.spatial import Delaunay

n_angles = 36
n_radii = 8

radii = numpy.linspace(1.0 / n_radii, 1.0, n_radii)
angles = numpy.linspace(0, 2*numpy.pi, n_angles, endpoint=False)
angles = numpy.repeat(angles[..., numpy.newaxis], n_radii, axis=1)

x = numpy.append(0, (radii*numpy.cos(angles)).flatten())
y = numpy.append(0, (radii*numpy.sin(angles)).flatten())

points = numpy.stack((x,y), axis=-1)
triangles = Delaunay(points).simplices

m_alugrid = dune.fem.grid.get("ALUSimplexGrid", dimgrid=2, refinement="conforming")
alugrid = m_alugrid.LeafGrid(m_alugrid.makeSimplexGrid(points, triangles))

output = alugrid.vtkWriter()
output.write("delaunay-demo")
