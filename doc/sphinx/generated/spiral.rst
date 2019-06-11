
Spiral Wave
===========

This demonstrates the simulation of spiral waves in an excitable media.
It consists of system of reaction diffusion equations with two
components. Both the model parameters and the approach for discretizing
the system are taken from
http://www.scholarpedia.org/article/Barkley\_model.

We use the *Barkley model* in its simplest form:

.. raw:: latex

   \begin{align*}
     \frac{\partial u}{\partial_t}
          &= \frac{1}{\varepsilon}f(u,v) + \Delta u \\
     \frac{\partial v}{\partial_t} &= h(u,v)
   \end{align*}

where

.. raw:: latex

   \begin{gather}
     f(u,v(=u\Big(1-u\Big)\Big(u-\frac{v+b}{a}\Big)
   \end{gather}

The function :math:`h` can take different forms, e.g., in its simplest
form

.. raw:: latex

   \begin{gather}
     h(u,v) = u - v~.
   \end{gather}

Finally, :math:`\varepsilon,a,b` for more details on how to chose these
parameters check the web page provided above.

We employ a carefully constructed linear time stepping scheme for this
model: let :math:`u^n,v^n` be given functions approximating the solution
at a time :math:`t^n`. To compute approximations :math:`u^{m+1},v^{m+1}`
at a later time :math:`t^{n+1}=t^n+\tau` we first split up the non
linear function :math:`f` as follows:

.. raw:: latex

   \begin{align*}
     f(u,v) = f_I(u,u,v) + f_E(u,v)
   \end{align*}

where using :math:`u^*(V):=\frac{V+b}{a}`:

.. raw:: latex

   \begin{align*}
     f_I(u,U,V) &= \begin{cases}
       u\;(1-U)\;(\;U-U^*(V)\;) & U < U^*(V) \\
       -u\;U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
     \end{cases} \\
   \text{and} \\
       f_E(U,V) &= \begin{cases}
       0 & U < U^*(V) \\
       U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
     \end{cases} \\
   \end{align*}

Thus :math:`f_I(u,U,V) = -m(U,V)u` with

.. raw:: latex

   \begin{align*}
     m(U,V) &= \begin{cases}
       (U-1)\;(\;U-U^*(V)\;) & U < U^*(V) \\
       U\;(\;U-U^*(V)\;)    & U \geq U^*(V)
     \end{cases}
   \end{align*}

Note that :math:`u,v` are assumed to take values only between zero and
one so that therefore :math:`m(u^n,v^n) > 0`. Therefore, the following
time discrete version of the Barkley model has a linear, positive
definite elliptic operator on its left hand side:

.. raw:: latex

   \begin{align*}
     -\tau\Delta u^{n+1} +
      (1+\frac{\tau}{\varepsilon} m(u^n,v^n))\; u^{n+1}
          &= u^n + \frac{\tau}{\varepsilon} f_E(u^n,v^n) \\
     v^{n+1} &= v^n + \tau h(u^n,v^n)
   \end{align*}

Which can now be solved using a finite element discretization for
:math:`u^n,v^n`.

Note that by taking the slow reaction :math:`h(u,v)` explicitly, the
equation for :math:`v^{n+1}` is purely algebraic. We will therefore
construct a scalar model for computing :math:`u^{n+1}` only and compute
:math:`v^{{n+1}}` be using the interpolation method on the space applied
to :math:`v^n + \tau h(u^n,v^n)`.

Let's get started by importing some standard python packages, ufl, and
some part of the dune-fempy package:

.. code:: ipython3

    import math
    import ufl
    import dune.ufl
    import dune.grid
    import dune.fem

In our attempt we will discretize the model as a 2x2 system. Here are
some possible model parameters and initial conditions (we even have two
sets of model parameters to choose from):

.. code:: ipython3

    dimRange   = 1
    dt         = 0.25
    linearSpiral = True
    
    if linearSpiral:
        spiral_a   = 0.75
        spiral_b   = 0.02
        spiral_eps = 0.02
        spiral_D   = 1./100
        def spiral_h(u,v): return u - v
    else:
        spiral_a   = 0.75
        spiral_b   = 0.0006
        spiral_eps = 0.08
        def spiral_h(u,v): return u**3 - v
    
    initial_u = lambda x: [1   if x[1]>1.25 else 0]
    initial_v = lambda x: [0.5 if x[0]<1.25 else 0]

Now we set up the reference domain, the Lagrange finite element space
(second order), and discrete functions for :math:`(u^n,v^n(`,
:math:`(u^{n+1},v^{n+1})`:

.. code:: ipython3

    # domain = dune.grid.structuredGrid([0,0],[3.5,3.5],[40,40])
    gridView = dune.grid.structuredGrid([0,0],[2.5,2.5],[30,30])
    space    = dune.fem.space.lagrange( gridView, dimRange=dimRange, order=1 )
    
    uh   = space.interpolate( initial_u, name="u" )
    uh_n = uh.copy()
    vh   = space.interpolate( initial_v, name="v" )
    vh_n = vh.copy()

We define the model in two steps: - first we define the standard parts,
not involving :math:`f_E,f_I`: - then we add the missing parts with the
required *if* statement directly using C++ code

.. code:: ipython3

    u   = ufl.TrialFunction(space)
    phi = ufl.TestFunction(space)
    
    # right hand sie (time derivative part + explicit forcing in v)
    a_ex = ufl.inner(uh_n, phi) * ufl.dx
    # left hand side (heat equation in first variable + backward Euler in time)
    a_im = (dt * spiral_D * ufl.inner(ufl.grad(u), ufl.grad(phi)) +
            ufl.inner(u,phi)) * ufl.dx
    
    ustar = (vh_n[0]+spiral_b)/spiral_a
    a_ex += ufl.conditional(uh_n[0]<ustar, dt/spiral_eps* u[0]*(1-uh_n[0])*(uh_n[0]-ustar),
                                         dt/spiral_eps*uh_n[0]*(1-u[0]) *(uh_n[0]-ustar) ) * phi[0] * ufl.dx
    
    equation   = a_im == a_ex
    ode_update = ufl.as_vector([ vh_n[0] + dt*spiral_h(uh_n[0], vh_n[0]) ])

The model is now completely implemented and can be created, together
with the corresponding scheme:

.. code:: ipython3

    solverParameters =\
           {"newton.tolerance": 1e-3,
            "newton.verbose": False,
            "newton.linear.tolerance": 1e-5,
            "newton.linear.preconditioning.method": "ilu",
            "newton.linear.verbose": False}
    scheme = dune.fem.scheme.galerkin( equation, space, solver="cg", parameters=solverParameters)

To show the solution we make use of the *animate* module of
*matplotlib*. Here is the ``stepping`` functions:

.. code:: ipython3

    def init():
        data = uh.pointData(1)
        C = plt.tricontourf(triangulation, data[:,0], cmap=plt.cm.rainbow, levels=levels)
        return C.collections
    def animate(count):
        global t,stepsize,nextstep
        nextstep += stepsize
        # print(count,t,stepsize,nextstep)
        while t < nextstep:
            uh_n.assign(uh)
            vh_n.assign(vh)
            info = scheme.solve(target=uh)
            vh.interpolate( ode_update )
            # print("Computing solution a t = " + str(t + dt), "iterations: " + str(info["linear_iterations"]) )
            t += dt
        data = uh.pointData(1)
        C = plt.tricontourf(triangulation, data[:,0], cmap=plt.cm.rainbow, levels=levels)
        # gridView.writeVTK("spiral", pointdata=[uh], number=count)
        return C.collections

And generate the movie:

.. code:: ipython3

    import matplotlib.pyplot as plt
    from matplotlib import animation, rc
    
    from numpy import linspace
    fig, ax = plt.subplots()
    ax.set_xlim(( 0, 2.5))
    ax.set_ylim(( 0, 2.5))
    triangulation = gridView.triangulation(1)
    levels = linspace(-0.1, 1.1, 256)
    ax.set_aspect('equal');
    t        = 0.
    stepsize = 0.5
    nextstep = 0.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=20, interval=100, blit=True)
    
    try:
        movie = anim.to_html5_video()
        from IPython.display import HTML, display
        display( HTML(movie) )
    except: # ffmpeg probably missing
        anim.save("spiral.html")
        try:
            from IPython.display import IFrame
            IFrame(src='./nice.html', width=700, height=600)
        except:
            pass



.. raw:: html

    <video width="432" height="288" controls autoplay loop>
      <source type="video/mp4" src="data:video/mp4;base64,AAAAHGZ0eXBNNFYgAAACAGlzb21pc28yYXZjMQAAAAhmcmVlAABOi21kYXQAAAKuBgX//6rcRem9
    5tlIt5Ys2CDZI+7veDI2NCAtIGNvcmUgMTQ4IHIyNjQzIDVjNjU3MDQgLSBILjI2NC9NUEVHLTQg
    QVZDIGNvZGVjIC0gQ29weWxlZnQgMjAwMy0yMDE1IC0gaHR0cDovL3d3dy52aWRlb2xhbi5vcmcv
    eDI2NC5odG1sIC0gb3B0aW9uczogY2FiYWM9MSByZWY9MyBkZWJsb2NrPTE6MDowIGFuYWx5c2U9
    MHgzOjB4MTEzIG1lPWhleCBzdWJtZT03IHBzeT0xIHBzeV9yZD0xLjAwOjAuMDAgbWl4ZWRfcmVm
    PTEgbWVfcmFuZ2U9MTYgY2hyb21hX21lPTEgdHJlbGxpcz0xIDh4OGRjdD0xIGNxbT0wIGRlYWR6
    b25lPTIxLDExIGZhc3RfcHNraXA9MSBjaHJvbWFfcXBfb2Zmc2V0PS0yIHRocmVhZHM9NiBsb29r
    YWhlYWRfdGhyZWFkcz0xIHNsaWNlZF90aHJlYWRzPTAgbnI9MCBkZWNpbWF0ZT0xIGludGVybGFj
    ZWQ9MCBibHVyYXlfY29tcGF0PTAgY29uc3RyYWluZWRfaW50cmE9MCBiZnJhbWVzPTMgYl9weXJh
    bWlkPTIgYl9hZGFwdD0xIGJfYmlhcz0wIGRpcmVjdD0xIHdlaWdodGI9MSBvcGVuX2dvcD0wIHdl
    aWdodHA9MiBrZXlpbnQ9MjUwIGtleWludF9taW49MTAgc2NlbmVjdXQ9NDAgaW50cmFfcmVmcmVz
    aD0wIHJjX2xvb2thaGVhZD00MCByYz1jcmYgbWJ0cmVlPTEgY3JmPTIzLjAgcWNvbXA9MC42MCBx
    cG1pbj0wIHFwbWF4PTY5IHFwc3RlcD00IGlwX3JhdGlvPTEuNDAgYXE9MToxLjAwAIAAAAl1ZYiE
    AD///vdonwKbWkN6gOSVxSXbT4H/q2dwfI/pAwAAAwAArqxz6Lgn8I5StQGOF/GUYt8JbmbKZ4bG
    2Itibainkk8GrGj6YFVJ03/FZ4WgY43LmhPeiNvFDlgnksHNJpfzgAAAVw8w3vT3ZepAn+udIgf8
    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
    AB1kYXRhAAAAAQAAAABMYXZmNTYuNDAuMTAx
    ">
      Your browser does not support the video tag.
    </video>



.. image:: spiral_files/spiral_13_1.png


... if the movie is not showing you might have to rerun the notebook ...
