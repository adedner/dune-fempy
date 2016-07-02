from __future__ import absolute_import
import ufl

class Space(ufl.VectorElement):
    def __init__(self, dimDomain, dimRange, dimWorld=None):
        if not dimWorld:
            dimWorld = dimDomain
        dimWorld = int(dimWorld)
        if dimDomain == 1:
            cell = ufl.Cell("interval", dimWorld)
        elif dimDomain == 2:
            cell = ufl.Cell("triangle", dimWorld)
        elif dimDomain == 3:
            cell = ufl.Cell("tetrahedron", dimWorld)
        else:
            raise NotImplementedError('dune.ufl.Space not implemented for dimension' + str(dimDomain) + '.')
        ufl.VectorElement.__init__(self, "Lagrange", cell, 1, dimRange)
