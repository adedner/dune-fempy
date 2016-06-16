"""Can be used to make and test a model using a modelfile.py file as input.
"""

import getopt, math, os, sys
from mpi4py import MPI

def main(argv):
   inputname = ''
   try:
      opts, args = getopt.getopt(argv,"htm",["help","test","make","dgf=","solver="])
   except getopt.GetoptError:
      print('usage: test.py -m modelname modelfile.py')
      sys.exit(2)

   test = 0
   make = 0
   dgf = "../data/unitcube-2d.dgf"
   solver = "fem"
   for opt, arg in opts:
      if opt == '-h':
         print('usage: generatemodel.py modelfile.py')
         sys.exit(2)
      elif opt in ("-t", "--test"):
          test = 1
      elif opt in ("-m", "--make"):
          make = 1
      elif opt in ("--dgf"):
          dgf = arg
      elif opt in ("--solver"):
          solver = arg

   if not len(args) == 1:
      print ( 'usage: generatemodel.py modelfile.py')
      sys.exit(2)
   inputname = args[0]
   (path,name) = os.path.split(inputname)
   if name.endswith('.py'):
      name = name[:-3]
   if "." in name:
      print("error: remove the extension from ", args[0])
      sys.exit(2)

   sys.dont_write_bytecode = True
   # read in the ufl to model translater
   print('input script used:' , name , 'from ' , path)
   sys.path.append(path)
   Model = __import__(name)
   model = Model.model

   if make:
       import dune.fem.grid as grid
       grid2d = grid.leafGrid(dgf, "YaspGrid", dimgrid=2)
       model.make(grid2d)
       grid1d = grid.get("OneDGrid")
       model.make(grid1d)
   elif test:
       import dune.fem.grid as grid
       import dune.fem.scheme as scheme
       import dune.fem.space as space
       # set up a 2d grid
       grid2d = grid.get("YaspGrid", dimgrid=2)
       Model = model.makeAndImport(grid2d)
       m = Model.get()
       g = grid.leafGrid(dgf, grid2d)
       print('get space')
       dimR = m.getDimRange()
       sp = space.create("Lagrange", g, dimrange=dimR )

       print('get scheme')
       try:
         femSchemeModule = scheme.get("FemScheme", sp, g, dimR, polorder=1, solver=solver)
       except Exception as exception:
          print('could not compile an extension module')
          print(exception)
          # try default fem solvers
          femSchemeModule = scheme.get("FemScheme", sp, g, dimR, polorder=1)

       s = femSchemeModule.Scheme( g, m.wrap(), "solution" )
       s1 = femSchemeModule.Scheme( g, m.wrap(), "solution" )
       s.solve(True)
       s1.solve(True)
       error = s.error()
       est = 1 # s.estimate()
       print( "difference between the errors of two schemes which are the same: ", s.error()-s1.error() )
       s1 = 0
       print( 0, error, -1, est, -1)
       sol = s.solution()
       for n in range(1, 4):
          g.globalRefine(1)
          olderror = error
          oldest = est
          sol = "hallo"
          s.solve()
          error = s.error()
          est = 1 # s.estimate()
          print(n, error, math.log(error/olderror)/math.log(0.5),\
                   est, math.log(est/oldest)/math.log(0.5))
          w = "hallo"  # this is to check memory management
          sol = s.solution()
       g = "not a grid anymore" # let's check memory management a second time

       # can we do the whole thing twice?
       g = grid.leafGrid(dgf, "YaspGrid", dimgrid=2)
       m = Model.get()
       s = scheme.scheme( "FemScheme", g, m, "solution", solver="fem" )
       s.solve()
       print("second scheme: ", s.error())

       # can we do the whole thing again?
       s = scheme.scheme( "FemScheme", g, m, "solution", solver="fem", polorder=2 )
       # these are not needed anymore
       m = 0
       g = 0
       # but this still needs to work
       s.solve()
       print("second scheme with second order: ", s.error())

       print( "finalize" )

if __name__ == "__main__":
   main(sys.argv[1:])
   print("finalize II")
