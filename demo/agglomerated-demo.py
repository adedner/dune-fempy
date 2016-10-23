"""Solve the Laplace equation
"""
from __future__ import print_function

import math
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d, cKDTree

from ufl import *

import dune.ufl
import dune.fem
import dune.fem.function as gf

import dune.create as create

dune.fem.parameter.append("../data/parameter")

def plot(grid, solution):
    try:
        from matplotlib import pyplot
        from numpy import amin, amax, linspace

        triangulation = grid.triangulation(4)
        data = solution.pointData(4)

        levels = linspace(amin(data[:,0]), amax(data[:,0]), 256)

        pyplot.gca().set_aspect('equal')
        pyplot.triplot(grid.triangulation(), antialiased=True, linewidth=0.2, color='black')
        pyplot.tricontourf(triangulation, data[:,0], cmap=pyplot.cm.rainbow, levels=levels)
        pyplot.show()
    except ImportError:
        pass


def compute(grid,NX,NY):
    voronoi_points = np.random.rand(NX*NY, 2)
    voronoi_kdtree = cKDTree(voronoi_points)
    vor = Voronoi(voronoi_points)
    voronoi_plot_2d(vor).savefig("agglomerate_voronoi"+str(NX*NY)+".pdf", bbox_inches='tight')
    ind = set()
    def agglomerate(en):
        p = en.geometry.center
        if 0: # Cartesian
            x = p[0]
            y = p[1]
            nx = int(x * NX)
            ny = int(y * NY)
            index = nx*NY+ny
        else:  # Voronoi
            test_point_dist, test_point_regions = voronoi_kdtree.query([p], k=1)
            index = test_point_regions[0]
        ind.add(index)
        return index
    spc  = create.space("AgglomeratedDG", grid, agglomerate, dimrange=1, order=3, storage="istl")
    assert len(ind)==NX*NY, "missing or too many indices provided by agglomoration object. Should be "+str(NX*NY)+" was "+str(len(ind))

    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())
    exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
    a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
    # a = a + 20./(u[0]*u[0]+1.) * v[0] * dx
    model = create.model("elliptic", grid, a==0, exact=exact, dirichlet={ 1:exact } )

    scheme = create.scheme("dggalerkin", spc, model, 10)
    df = scheme.solve(name="solution")
    df_interpol = spc.interpolate( create.function("ufl",grid,"exact",4,exact), "exact" )

    df_coeff = dune.ufl.GridCoefficient(df)
    l2error_gf = create.function("ufl", grid, "error", 5,
            as_vector([(exact[0]-df_coeff[0])**2]) )
    error = math.sqrt( l2error_gf.integrate() )
    print("size:",spc.size,"L2-error:", error)
    # grid.hierarchicalGrid.globalRefine(1) <- this leads to not responding program possibly due to disabled RP
    grid.writeVTK("agglomerate"+str(NX*NY),
        celldata=[ df, df_interpol, create.function("local",grid,"cells",1,lambda en,x: [agglomerate(en)]) ])

grid = create.view("adaptive", grid="ALUCube", constructor=dune.grid.cartesianDomain([0, 0], [1, 1], [100, 100]), dimgrid=2)
compute(grid,3,3)
compute(grid,6,6)
compute(grid,12,12)
