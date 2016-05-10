from __future__ import print_function

import dune.fem.grid as grid
import dune.fem.space as space

yaspGrid = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)
lagrangeSpace = space.create("Lagrange", yaspGrid)
