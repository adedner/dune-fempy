"""Solve the Laplace equation
"""
from __future__ import print_function

import math
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

def compute():
    grid = create.view("adaptive", grid="ALUConform", constructor=dune.grid.cartesianDomain([0, 0], [1, 1], [8, 8]), dimgrid=2)
    # spc  = create.space("DGONB", grid, dimrange=1, order=1)
    spc  = create.space("AgglomeratedDG", grid, lambda en: 0, dimrange=1, order=1)
    expr_global = lambda x: [-(x[1] - 0.5)*math.sin(x[0]*12)]*dimRange

    gf = create.function("discrete", spc, "expr_global", 1, expr_global)
    grid.writeVTK("test", pointdata=[gf])

compute()
