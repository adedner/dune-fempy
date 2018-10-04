from __future__ import absolute_import, division, print_function, unicode_literals

import dune.fem.scheme

def create(velocitySpace, pressureSpace, model, name, viscosity, timeStep, **kwargs):
    stokesScheme = dune.fem.scheme.create( "Stokes", ( velocitySpace, pressureSpace ), model, name+"Stokes",\
                   viscosity, timeStep, storage = "Istl" )
    burgersScheme = dune.fem.scheme.create( "Burgers", ( velocitySpace, pressureSpace ), model, name+"Burgers",\
                    viscosity, timeStep, storage = "Istl" )
    return stokesScheme, burgersScheme
