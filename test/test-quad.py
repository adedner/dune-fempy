from __future__ import print_function, unicode_literals
import time, math, numpy
import dune
from dune.generator import algorithm

from dune.grid import cartesianDomain, yaspGrid
domain = cartesianDomain([0, 0], [1, 0.25], [60, 16])
yaspView = yaspGrid(domain)

@dune.grid.gridFunction(yaspView)
def function(x):
    return numpy.cos(2.*numpy.pi/(0.3+abs(x[0]*x[1])))
def interpolate(grid):
    mapper = grid.mapper({dune.geometry.vertex: 1})
    data = numpy.zeros(mapper.size)
    for v in grid.vertices:
        data[mapper.index(v)] = function(v.geometry.center)
    return mapper, data
mapper, data = interpolate(yaspView)
@dune.grid.gridFunction(yaspView)
def p12dEvaluate(element,x):
    indices = mapper(element)
    bary = 1-x[0]-x[1], x[0], x[1]
    return sum( bary[i] * data[indices[i]] for i in range(3) )

@dune.grid.gridFunction(yaspView)
def error(element,x):
    return numpy.abs(p12dEvaluate(element,x)-function(element,x))

rules = dune.geometry.quadratureRules(5)

start = time.time()
l2norm2 = 0
for e in yaspView.elements:
    hatxs, hatws = rules(e.type).get()
    weights = hatws * e.geometry.integrationElement(hatxs)
    l2norm2 += numpy.sum(error(e, hatxs)**2 * weights, axis=-1)
print("Python:",math.sqrt(l2norm2))
print("time used:", round(time.time()-start,2))

#algo = algorithm.load('l2norm2', 'test_quad.hh', yaspView, rules, error)
algo = algorithm.load('l2norm2FemQuad', 'test_quad.hh', yaspView, rules, error)
start = time.time()
l2norm2 = algo(yaspView,rules,error)
print("C++:",math.sqrt(l2norm2))
print("time used:", round(time.time()-start,2))

try:
    import dune.geometry.quadpy as quadpy
    rules = quadpy.rules({dune.geometry.quadrilateral: ("C2 7-2","Stroud")})

    start = time.time()
    l2norm2 = 0
    for e in yaspView.elements:
        hatxs, hatws = rules(e.type).get()
        weights = hatws * e.geometry.integrationElement(hatxs)
        l2norm2 += numpy.sum(error(e, hatxs)**2 * weights, axis=-1)
    print("Python:",math.sqrt(l2norm2))
    print("time used:", round(time.time()-start,2))

    start = time.time()
    l2norm2 = algo(yaspView,rules,error)
    print("C++:",math.sqrt(l2norm2))
    print("time used:", round(time.time()-start,2))
except ImportError:
    pass
