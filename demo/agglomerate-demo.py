"""Solve the Laplace equation
"""
from __future__ import print_function

import math
import numpy
from scipy.spatial import Voronoi, voronoi_plot_2d, cKDTree

from ufl import *

import dune.ufl
import dune.fem
import dune.fem.function as gf

import dune.create as create

from voronoi import triangulated_voronoi

dune.fem.parameter.append("../data/parameter")

bounding_box = numpy.array([0., 1., 0., 1.]) # [x_min, x_max, y_min, y_max]

def compute(N,cube=False):
    numpy.random.seed(1234)
    voronoi_points = numpy.random.rand(N, 2)
    voronoi_kdtree = cKDTree(voronoi_points)
    vor = Voronoi(voronoi_points)
    voronoi_plot_2d(vor).savefig("agglomerate_voronoi"+str(N)+".pdf", bbox_inches='tight')
    ind = set()
    def agglomerate(en):
        p = en.geometry.center
        test_point_dist, test_point_regions = voronoi_kdtree.query([p], k=1)
        index = test_point_regions[0]
        ind.add(index)
        return index

    points, triangles = triangulated_voronoi(voronoi_points, bounding_box)
    if cube:
        grid = create.view("ALUCube", constructor=dune.grid.cartesianDomain([0, 0], [1, 1], [100, 100]), dimgrid=2)
    else:
        grid = create.grid("ALUSimplex", {'vertex':points, 'simplex':triangles}, dimgrid=2)

    spc    = create.space("AgglomeratedDG", grid, agglomerate, dimrange=1, order=1, storage="istl")
    spcVEM = create.space("AgglomeratedVEM", grid, agglomerate, dimrange=1, order=1, storage="istl")
    assert len(ind)==N, "missing or too many indices provided by agglomoration object. Should be "+str(N)+" was "+str(len(ind))

    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1, field="double")
    x = SpatialCoordinate(uflSpace.cell())
    exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )

    df          = spcVEM.interpolate( create.function("ufl",grid,"exact",4,exact), "exactvem" )
    df_interpol = spc.interpolate( create.function("ufl",grid,"exact",4,exact), "exactdg" )

    df_coeff = dune.ufl.GridCoefficient(df)
    l2error_gf = create.function("ufl", grid, "error", 5,
            as_vector([(exact[0]-df_coeff[0])**2]) )
    error = math.sqrt( l2error_gf.integrate() )
    df_interpol_coeff = dune.ufl.GridCoefficient(df_interpol)
    l2error_gf = create.function("ufl", grid, "error", 5,
            as_vector([(exact[0]-df_interpol_coeff[0])**2]) )
    error_dg = math.sqrt( l2error_gf.integrate() )

    print("dg  size:",spc.size,   "L2-error:", error_dg)
    print("vem size:",spcVEM.size,"L2-error:", error)
    grid.writeVTK("agglomerate"+str(N),
        pointdata=[ df, df_interpol ],
        celldata= [ create.function("local",grid,"cells",1,lambda en,x: [agglomerate(en)]) ])

compute(9,False)
compute(36,False)
compute(144,False)
