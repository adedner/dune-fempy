.. _usageintro:

################################
Usage introduction
################################

As explained in the introduction, Dune-fempy provides a python interface for solving PDEs using Dune and Dune-Fem. Here we will explain how this can be used to set up the various parts of a numerical problem. For a complete example of how this works, see :ref:`usageexample`.

.. contents::

################################
Setting up a computational grid
################################

Grids can be set up easily in python using the database found in python/database/grid. This allows the user to specify the grid file they want to use, and the Dune type that the grid takes. An example of this in python is the following

.. code-block:: python

  grid2d = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)

An explanation of how leafGrid works is given in the following docstring

.. autofunction:: dune.fem.grid.leafGrid()

This relies on get(), which does the following

.. autofunction:: dune.fem.grid.get()

###############################################
Setting up a mathematical model using UFL
###############################################

In Dune-fempy, the **model** refers to the part of the problem that contains the PDE itself, the dimensions of the space and the boundary conditions. UFL is used to express the model, which is used to generate a Dune model file.

The following docstrings give an in-depth explanation of a file that performs the model generation.

femufl
------

.. automodule:: dune.models.femufl
   :members: generate
   :undoc-members:

Once a model has been generated using the above, it can be used in python code by calling

.. code-block:: python

  Model = model.makeAndImport(grid2d)
  pyModel = Model.get()

here *grid2d* must be previously defined in the code.

################################
Setting up a numerical scheme
################################

In Dune-fempy, the **scheme** contains information about the method used to solve the PDE. Schemes can be set up in a similar way to grids using the database found in python/database/scheme. An example of this in python is the following

.. code-block:: python

  scheme = scheme.scheme("FemScheme", grid, model, "solution", solver="fem")

Here *grid* and *model* must both be previously defined, as shown above. An explanation of how scheme works is given in the following docstring

.. autofunction:: dune.fem.scheme.scheme()

Here as before, the function get() creates the C++ scheme class using the information given to it.

.. _usageexample:

################################
A full example
################################

Here we give a complete example for a problem that uses all the above methods. Other such examples can be found in the demo directory.

.. code-block:: python

  ############################################################################
  # initialize
  ############################################################################
  from __future__ import print_function
  import sys,os,math, getopt
  sys.path.append("../python")
  import timeit

  import sympy
  import ufl

  import dune.models.femufl as duneuflmodel
  import dune.fem.gridfunction as gf
  import dune.fem.grid as grid
  import dune.fem.scheme as scheme

  grid2d = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2)
  grid2d.globalRefine(3)

  ############################################################################
  # build a transport model -eps Laplace(u) + velocity.grad(u) + gamma u = f
  # with unknown velocity field
  ############################################################################
  model    = duneuflmodel.DuneUFLModel(2,1,'Transport')
  vecmodel = duneuflmodel.DuneUFLModel(2,2,'Velocity')
  exact = [sympy.cos(2*math.pi*model.x0)*sympy.cos(2*math.pi*model.x1)]
  u = model.trialFunction()
  v = model.testFunction()
  x = model.spatialCoordinate()
  velo = vecmodel.coefficient('velocity')
  diff = model.coefficient('diffusion')
  a = ( (ufl.dot(velo,ufl.grad(u[0]))+0.01*u[0])*v[0] + diff[0]*ufl.inner(ufl.grad(u[0]),ufl.grad(v[0])) ) * ufl.dx(0) 
  L = 10./(1.+(x[0]*x[0]+x[1]*x[1])**4 )  *  v[0]*ufl.dx(0)
  model.setCoefficient("diffusion",[1+0.1*model.x0*model.x0])

  start_time = timeit.default_timer()
  model.generate(a,L,exact)
  Model = model.makeAndImport(grid2d)
  print("Building TransportModel took ", timeit.default_timer() - start_time, "s")
  m = Model.get()

  s = scheme.scheme("FemScheme", grid2d, m, "transport", polorder=1)

  ############################################################################
  # Now use different way to solve this problem with given velocity field
  ############################################################################

  print("define velocity field in global coordinates (using different approaches)")
  expr1 = gf.MathExpression(["-(x1-0.5)","x0-1./2."])    
  expr2 = gf.SympyExpression(["-(x1-0.5)","x0-1./2."])   
  def expr_func(x,r):
      r[0] = -(x[1]-0.5)
      r[1] = (x[0]-0.5)
  expr3 = gf.FuncExpression(2,expr_func)
  velocityGlobal = grid2d.getGlobal("global_velocity",expr1)

  m.setvelocity(velocityGlobal)
  s.solve()
  out_gl = grid2d.vtkOutput()
  out_gl.add(s.solution())
  print('solution: ', s.solution())
  out_gl.add(velocityGlobal)
  out_gl.write( "testmodel_global" )


