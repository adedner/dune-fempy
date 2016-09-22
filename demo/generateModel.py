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
   ModelDescription = __import__(name)
   model = ModelDescription.model

   if make:
       import dune.grid as grid
       grid2d = dune.grid.create("Yasp", dgf, dimgrid=2)
       model.make(grid2d)
       # grid1d = dune.grid.create("OneD", dgf, dimgrid=2)
       # grid1d = grid.get("OneDGrid")
       # model.make(grid1d)
   elif test:
       import dune.grid as grid
       import dune.fem as fem
       import dune.fem.scheme as scheme
       import dune.fem.space as space
       # set up a 2d grid
       grid2d = grid.get("ALUSimplexGrid", dimgrid=2)
       Model = model.makeAndImport(grid2d)
       m = Model.get()
       if hasattr(ModelDescription,"DGF"):
          print("using: ",ModelDescription.DGF)
          g = dune.grid.create("ALUSimplex", grid.string2dgf(ModelDescription.DGF), dimgrid=2)
       else:
          g = dune.grid.create("ALUSimplex", dgf, dimgrid=2)
       print('get space')
       dimR = m.dimRange
       sp = space.create( "Lagrange", g, dimrange=dimR )

       print('get scheme')
       try:
         femSchemeModule = scheme.get( "FemScheme", sp, solver=solver )
       except Exception as exception:
          print('could not compile an extension module')
          print(exception)
          # try default fem solvers
          femSchemeModule = scheme.get( "FemScheme", sp)

       s   = femSchemeModule.Scheme( sp, m, "solution" )
       s1  = femSchemeModule.Scheme( sp, m, "solution" )
       vtk = g.vtkWriter()
       solution = s.solve()
       solution.addToVTKWriter(vtk, vtk.PointData)
       solution1 = s1.solve()
       error = s.error(solution)
       est = 1 # s.estimate()
       print( "difference between the errors of two schemes which are the same: ",\
               s.error(solution)-s1.error(solution1) )
       s1 = 0
       print( 0, error, -1, est, -1)
       vtk.write("generateModel"+str(0));
       for n in range(1, 4):
          g.hierarchicalGrid.globalRefine(1)
          olderror = error
          oldest = est
          solution.clear()
          s.solve( None, solution )
          # s.solve( target=solution )
          error = s.error(solution)
          est = 1 # s.estimate()
          print(n, error, math.log(error/olderror)/math.log(0.5),\
                   est, math.log(est/oldest)/math.log(0.5))
          w = "hallo"  # this is to check memory management
          vtk.write("generateModel"+str(n));
       g = "not a grid anymore" # let's check memory management a second time

       # can we do the whole thing twice?
       g = dune.grid.create("ALUSimplex", dgf, dimgrid=2)
       m = Model.get()
       s = scheme.create( "FemScheme", sp, m, "solution", solver="fem" )
       print("second scheme: ", s.error(s.solve()))

       # can we do the whole thing again?
       # sp = space.create("Lagrange", g, dimrange=dimR, order=1 )
       sp = space.create( "Lagrange", g, dimrange=dimR )
       s = scheme.create( "FemScheme", sp, m, "solution", solver="fem" )
       # these are not needed anymore
       m = 0
       g = 0
       sp = 0
       # but this still needs to work
       solution = s.solve()
       print("second scheme with second order: ", s.error(solution))

       print( "finalize" )

if __name__ == "__main__":
   main(sys.argv[1:])
   print("finalize II")
