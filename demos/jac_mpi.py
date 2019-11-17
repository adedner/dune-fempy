# <markdowncell>
# # Dense Vectors and the Geometry Classes
# A quick survey of Dune-Common and Dune-Geometry
# The core module Dune-Common provides some classes for dense linear algebra.
# The `FieldVector` and `FieldMatrix` classes are heavily used in the grid
# geometry realizations. The conceptional basis for these geometries is
# provided by Dune-Geometry, providing, for example, reference elements and
# quadrature rules.
# <codecell>

import numpy as np
import pylab
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

global_size = 500
dx = 1./(global_size)
print("dx = ", dx)

# DoF per rank
N = int(global_size / comm.size)

nstart = rank * N
a = nstart * dx + 0.5*dx
b = a + (N-1) * dx

h = np.linspace(a, b, N) # 'N' linearly spaced numbers
u = np.sin( h * np.pi )
b = u * np.pi * np.pi * dx * dx

xL = 0.
xR = 0.

L  = -1.
D  =  2.
U  = -1.

x = np.zeros_like(b)

parallel = True

converged = False

links = []
if comm.size > 1:
    # left side
    if rank > 0:
        links.append( [rank-1,0] )
    else:
        links.append( [-1,-1] )

    # right side
    if rank < comm.size-1:
        links.append( [rank+1, N-1] )
    else:
        links.append( [-1,-1] )

def synchronize( vec ):
    global links

    for link in links:
        if link[0] != -1:
            data = {'a': vec[ link[1] ] }
            #print("Rank ",rank," sending to ",link, " v=",value)
            req = comm.isend(data, dest=link[0] )
            req.wait()

    values = [ 0, 0]
    c = 0
    for link in links:
        if link[0] != -1:
            req2 = comm.irecv(source=link[0])
            data = req2.wait()
            values[ c ] = data["a"]
        c += 1
    return values

def globalSum( localErr ):
    err = np.array( [ localErr ] )
    #print("Rank ", rank, " err = ", localErr )
    errG = np.zeros_like( err )
    comm.Allreduce( err, errG, op=MPI.SUM)
    #print("Rank ", rank, " global err = ", err[0], errG[ 0 ] )
    return errG[ 0 ]


def jacobi(b, x, x_new):
    global L, D, U
    global xR
    global xL

    N = len(x)

    omega = 2./3.
    x_new[ 0   ] = omega * (b[ 0 ]   - L * xL - U * x[1] ) / D + (1.0 -omega) * x[0]
    for i in range(1,N-1):
        x_new[ i ] = omega * (b[ i ] - L * x[ i-1 ] - U * x[ i+1 ] ) / D + (1.0 -omega) * x[i]
    x_new[ N-1 ] = omega * (b[ N-1 ] - L * x[ N-2] - U * xR   ) / D + (1.0 -omega) * x[N-1]

### end Jacobi ###


def solve(b):
    global xL
    global xR

    converged = False

    N = len(b)

    x = np.zeros_like( b )
    x_new = np.zeros_like( x )

    for it_count in range(global_size*global_size):
        if converged:
            print("Iteration %s: converged = %s",it_count, converged)
            break

        # compute multi-threaded iteration
        jacobi( b, x, x_new )

        # update solution and compute error
        x_k = x - x_new
        err = globalSum( np.dot(x_k, x_k) )

        # update x with values of x_new
        np.copyto( x, x_new )

        # if | x - x_new | < 1e-9  stop iterating
        if err < 1e-10:
            converged = True

        if rank == 0 and it_count % 1000 == 0:
            print("Iteration %s: converged = %s, err = %s",it_count, converged, err )

        if comm.size > 1:
            values = synchronize( x )
            xL = values[ 0 ]
            xR = values[ 1 ]

    comm.barrier()
    return x

#### end of solve ##############
def gatherSolution( x, root=0 ):
    # compose plot
    recvbuff = None
    if rank == root:
        recvbuff = np.empty([comm.size, N], dtype='d')
    comm.Gather(x, recvbuff, root=root)

    X = None
    # unpack buffer
    if rank == root:
        X = np.empty( comm.size*N, dtype='d')
        for i in range(comm.size):
            for j in range(N):
                X[ i*N + j ] = recvbuff[ i ][ j ]
    return X
### end gatherSolution ###

# measure CPU time
start = time.time()

x = solve( b )

if rank == 0:
    print("time loop: %s sec",time.time()-start,flush=True)

# collect solution to rank 0
x = gatherSolution( x )

if rank == 0:
    h = np.linspace(0., 1., global_size) # 'N' linearly spaced numbers
    u = np.sin( h * np.pi )
    pylab.plot(h, x, 'co') # sin(x)/x
    pylab.plot(h, u, color='red') # sin(x)/x
    pylab.show() # show the plot
