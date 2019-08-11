import math
import ufl
from ufl import div, grad, inner, dot, dx, dS, jump, avg
import dune.ufl
import dune.grid
import dune.fem
import dune.alugrid

from dune.grid import cartesianDomain, gridFunction

# <markdowncell>
# In our attempt we will discretize the model as a 2x2 system. Here are some possible model parameters and initial conditions (we even have two sets of model parameters to choose from):
# <codecell>

maxLevel     = 13
startLevel   = 3
endTime      = 20.
saveInterval = 0.5
maxTol       = 1e-5

# <markdowncell>
# Now we set up the reference domain, the Lagrange finite element space (second order), and discrete functions for $(u^n,v^n($, $(u^{n+1},v^{n+1})$:
# <codecell>

domain   = dune.grid.cartesianDomain([0,0],[2.5,2.5],[5,5])
#baseView = dune.alugrid.aluConformGrid(domain)
#baseView = dune.grid.ugGrid(domain)
baseView = dune.grid.ugGrid("2dgrid.dgf", dimgrid=2, dimworld=2)

print("BaseView created")
gridView = dune.fem.view.adaptiveLeafGridView( baseView )
print("Adaptive view created")
baseView.hierarchicalGrid.globalRefine(startLevel)
