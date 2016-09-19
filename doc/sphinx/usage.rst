.. _usage:

.. raw:: html

  <h1> Usage guide </h1>

################################
Introduction
################################

In the :ref:`tutorial <tutorial>`, we looked at the basic functions available to Dune-Fempy in the context of the Laplace equation. Here we will explain in more detail how we set up the various parts of a numerical problem on the python side, and the tools we have at our disposal.

Behind all of the interface methods we use, the philosophy is that they are set up in a very similar way to the Dune-Fem structure. This should hopefully make the underlying C++ code transparent and easier to understand for a user of the python code, and vice versa. For more information about the C++ code, see the :ref:`advanced topics <advanced>` section.

.. contents::

################################
Setting up a computational grid
################################

In Dune-Fempy the **grid** (somewhat self-explanatorily) refers to the grid used in the numerical method. It contains information about the mesh used, the dimension of the domain, and the Dune type that the grid takes. Note that the `leafGrid` function always takes the following parameters

* Mesh (string): Either a dgf file or a description of the mesh (we saw the latter in the tutorial).
* Identifier (string): The grid type that Dune uses.

To explain further, the identifier looks into the grid database, which we will talk more about in :ref:`this section <database>`. For now, we will just say that the list of possible grids is located in python/database/grid, and that optional arguments may or may not be needed depending on what identifier we use. For instance a 1D grid on a unit cube could simply be defined using a dgf file as

.. code-block:: python

  onedgrid = dune.fem.leafGrid("../data/unitcube-1d.dgf", "OneDGrid")

On the other hand a grid with additional parameters might look like this

.. code-block:: python

  holegrid = dune.fem.leafGrid("../data/hole2.dgf", "ALUSimplexGrid", dimgrid=2, dimworld=1, refinement="conforming")

Refinements
"""""""""""

Once the grid has been created, it can be refined using

.. code-block:: python

  grid.hierarchicalGrid.globalRefine(2)

Here the inital mesh will be globally refined twice.

It is also possible to do adaptively refine the mesh (see e.g. afem.py in /demo).

.. _gridfunctions:

################################
Setting up functions on the grid
################################

Once we have set up the grid, it is possible to define **grid functions** on it using the method `function`. There are four different kinds of grid functions one can set up this way.

1. `globalExpr`: A function using global coordinates defined on the python side.
2. `localExpr`: A function using local coordinates defined on the python side.
3. `ufl`: A function defined using a UFL expression.
4. `code`: A function defined using C++ code that is inputted.

Let us first look at an example of the method itself, `function`. It always needs the following

* Name (string): The internal name python uses for the function.

It also takes the following keyword arguments

* order=(int): The order of the approximation (optional).
* type=(object): The type of function used (from the above list), followed by the object used to define it.

.. code-block:: python

    func = grid.function("global_velocity", order=0, globalExpr=somefunction)




###############################################
Setting up a space
###############################################

In Dune-Fempy the **space** refers to the function space used in our finite element method. The space can be set up in python in an identical way to the grid as follows.

.. code-block:: python

  space = dune.fem.create.space("Lagrange", grid, dimrange=1, polorder=2)

Interpolation
"""""""""""""

Once the space has been set up, we can use the method **interpolate** to create functions that can be accessed on the python side. This approach has the advantage that we can have a single solution vector to use for multiple schemes. We do this as follows

.. code-block:: python

  u = spc.interpolate(lambda x: [x[0]])
  scheme.solve(target = u)

In the first line we interpolate ``u`` over our space using a spatial coordinate ``x[0]`` (previously defined). This ``u`` can then be passed into our ``solve`` method by specifying a ``target``.

###############################################
Setting up a mathematical model using UFL
###############################################

In Dune-Fempy, the **model** refers to the part of the problem that contains the weak form of the PDE and its boundary conditions. UFL is used to express the PDE, and from this we can generate a Dune model file. The module generation is done in the file python/dune/models/elliptic.hh.

We have already seen an example of UFL usage in the tutorial, so here let's consider a more complex example, a model for mean curvature flow.

.. code-block:: python

  dt        = 0.0025
  theta     = 0.5

  uflSpace = dune.ufl.Space((surface.dimGrid, surface.dimWorld), surface.dimWorld)
  u = TrialFunction(uflSpace)
  v = TestFunction(uflSpace)
  u_n = Coefficient(uflSpace)

  a_im = (dt * theta * inner(grad(u), grad(v)) + inner(u, v)) * dx
  a_ex = (-dt * (1-theta) * inner(grad(u), grad(v)) + inner(u, v)) * dx
  lhsModel = dune.models.elliptic.importModel(surface, a_im == 0).get()
  rhsModel = dune.models.elliptic.importModel(surface, a_ex == 0).get()

As we can see, it is not very difficult to set up time-dependent problems since we can make separate models for the explicit and implicit parts of the equation. Constants such as `dt` and `theta` can be simply defined on the python side and put directly into the bilinear forms, and everything else can be acquired from the UFL side. `Coefficient` is a special variable that be set to different functions that we will talk more about below.

UFL code in Dune-Fempy is mostly identical to that in the original module, so the `documentation <http://fenicsproject.org/documentation/ufl/1.0-beta2/ufl.html>`_ is a useful resource.

Boundary conditions
"""""""""""""""""""

Boundary conditions can also be added to the model using UFL. Any *natural* boundary conditions (e.g. Neumann or Robin) can be added to the weak form directly by using a surface integral ds (instead of dx). On the other hand, *essential* boundary conditions can be added optionally using the **dirichlet** argument as follows.

.. code-block:: python

  g1 = [cos(x[0]), sin(x[0])]
  g2 = [x[1], 3]
  model = dune.models.elliptic.importModel(grid, a == b, dirichlet = {1:[g1], 2:[g2]}).get()

Here ``1:[g1]`` tells us that the function ``g1`` is set on the boundary assigned to ``1`` in the mesh file, and similarly ``2:[g2]`` sets boundary ``2`` to ``g2``. Multiple Dirichlet boundary conditions can be individually assigned to different boundaries in this way.

.. _coefficients:

Coefficients and constants
""""""""""""""""""""""""""

Suppose we want to have a scalar or vector in our model that can be set using either the solution from another scheme, or a function that we define ourselves in python. Suppose also that we might want to have a constant in our model that can be set to different values. We can do these things using the **Coefficient** and **Constant** variables as shown below (example taken from demo/heat.py).

First we set up an initial function via python called `initial` and set up two interpolated solutions with it.

.. code-block:: python

    initial = lambda x: [ math.atan( (10.*x[0]*(1-x[0])*x[1]*(1-x[1]))**2 ) ]
    solution = spc.interpolate(initial, name="u")
    old_solution = spc.interpolate(initial, name="u_n")

Next we set up a model using UFL in the usual way, only we also define a *Coefficient* `u_n` and a *Constant* `tau`. 

.. code-block:: python

    uflSpace = dune.ufl.Space((grid.dimGrid, grid.dimWorld), 1)
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    u_n = Coefficient(uflSpace)
    tau = Constant(triangle)
    a = (inner(u - u_n, v) + tau * inner(grad(theta*u + (1-theta)*u_n), grad(v))) * dx

Now when we create the model object, we add in an additional keyword argument which sets up the coefficient to our initial data. 

.. code-block:: python

    model = dune.fem.create.ellipticModel(grid, a == 0)(coefficients={u_n: old_solution})

Note that for each coefficient we defined in UFL, we must set its value using the dictionary format (`coefficients={coef1: value1, coef2: value2, coef3: value3}`). We can also set any constants in the same format (`constants={...}`). 

We also set up the scheme. 

.. code-block:: python

    scheme = dune.fem.create.scheme("FemScheme", spc, model, "scheme")

Finally before we solve the model, we must set our constant. Here we can do this at each step in a for loop, using `setConstant`.

.. code-block:: python

    deltaT = 0.01
    steps = int(1 / deltaT)
    for n in range(1, steps+1):
        model.setConstant(tau, [deltaT])
        old_solution.assign(solution)
        scheme.solve(target=solution)
        grid.writeVTK("heat", pointdata=[solution], number=n)

Note that we have to pass in `deltaT` using square brackets, since the method takes a list object. 

The advantage of using this method, is that it lets us change the model slightly without having to recompile it completely.

Note that we can also use coefficients and constants with :ref:`grid functions <gridfunctions>`. For ordinary grid functions the procedure is the same as above, but with grid functions created using C++ code, there are a couple of extra things to consider. See the following example.

.. code-block:: python
    
    func1 = """
    double s = sin(xGlobal[0]);
    double c = cos(@const:fac*xGlobal[1]);
    value[ 0 ] = s*s;
    value[ 1 ] = s*c;
    value[ 2 ] = c*c;
    """
    func2 = """
    double cx = cos(xGlobal[0]);
    double cy = cos(@const:fac*xGlobal[1]);
    double sx = sin(xGlobal[0]);
    double sy = sin(@const:fac*xGlobal[1]);
    value[ 0 ][ 0 ] = cx*sx*@gf:test[1];
    value[ 0 ][ 1 ] = 0;
    value[ 1 ][ 0 ] = cx*cy*@gf:test[0];
    value[ 1 ][ 1 ] = -@const:fac*sx*sy;
    value[ 2 ][ 0 ] = 0;
    value[ 2 ][ 1 ] = -2.*@const:fac*cy*sy;
    """
    code = { 'eval': func1, 'jac': func2 }

    coeffFunc = grid.function("global_velocity", order=1, globalExpr=lambda x: [1,2])
    func = grid.function("code", 3, code=code, coefficients={"test": coeffFunc}, constants={"fac": 1} )
    func.setConstant("fac", [factor])

Note that we have to add placeholders directly into the code that tell the compiler to replace them with coefficients/constants. Simply put,

1. For a coefficient, write `@gf:name`.
2. For a constant, write `@const:name`.

Then, to declare you are using this variable in your function, you would add `coefficients={"name", somefunction}` or `constants={"name"}`. `setConstant`. Constants can then be set using `func.setConstant("name", someConstant)`.

.. _dunemodel:

Stand-alone Dune model generation
"""""""""""""""""""""""""""""""""

It is possible to just create a C++ model file using UFL code for use within the Dune-Fem-Howto framework without using any of the other python interface tools. The advantage of this is to forgo the complicated process of manually writing a model file with functions for the source, flux, linSource, linFlux and so on. This can be done quite easily in the following way.

1. Create a UFL model file in a similar way to above. For examples of exactly what is required, see the *models* folder for reference.
2. Run the generateModel script in the build-cmake/demos directory. For example, to generate a model file for the transport equation example, you would run.

  .. code-block:: bash

    python generateModel.pyc ../../models/transport.py

  Optionally you can add -m or -t to the call to make a python module, or test it with a FEM scheme.
3. Use the generated model file in conjuction with your own Dune code to make a method. The file is outputted to *build-cmake/python/dune/generated* using the name given in the UFL file (e.g. *TransportModel.hh* in this case).

################################
Setting up a numerical scheme
################################

In Dune-Fempy, the **scheme** contains information about the method used to solve the PDE. Just as before, schemes can be set up in a similar way to grids and spaces using the database found in python/database/scheme. A simple example of this in python is the following.

.. code-block:: python

  scheme = dune.fem.create.scheme("FemScheme", space, model, "scheme")

Here `"FemScheme"` is the identifier for the default FEM scheme in Dune, and `space` and `model` are the Fempy objects defined above.

.. _usageexample:

################################
A full example
################################

Here we give a complete example for a problem that uses all the above methods (an example for Mean Curvature Flow found in demo/mcf.py). Other examples can be found in the demo directory.

.. literalinclude:: ../../demo/mcf.py
   :language: python
