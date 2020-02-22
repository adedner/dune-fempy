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
import multiprocessing
import time

global_size = 500
size = int(global_size)
dx = 1./size

h = np.linspace(0, 1, size) # 'size' linearly spaced numbers
u = np.sin( h * np.pi )
b = u * np.pi * np.pi * dx * dx

# matrix  3-point stencil
L  = -1.
D  =  2.
U  = -1.

nThreads = 4

# measure CPU time
start = time.time()

barrier  = multiprocessing.Barrier( nThreads )

def iterate(start, end, x, x_new ):
    global size

    omega   = 2./3.
    omega_1 = 1. - omega
    for i in range(start,end):
        s1 = 0.
        if i > 0:
            s1 = L * x[i-1]
        if i < size-1:
            s1 += U * x[i+1]
        x_new[ i ] = omega * ( b[ i ] - s1 ) / D + omega_1 * x[ i ]

def solve(num, X):
    global size
    global nThreads
    global barrier

    x     = np.zeros_like( b )
    x_old = np.zeros_like( b )
    x_new = np.zeros_like( b )

    # update local copy of X
    for i in range(size):
        x[ i ] = X[ i ]

    size_thread = int( size / nThreads )
    thStart = num * size_thread
    thEnd   = thStart + size_thread

    converged = False

    for it_count in range(size*size):
        if converged:
            break

        # store last iterate
        np.copyto( x_old, x )

        # compute multi-threaded iteration
        iterate( thStart, thEnd, x, x_new )

        # update global array
        for i in range(thStart, thEnd):
            X[ i ] = x_new[ i ]

        # wait all process to arrive here
        barrier.wait()

        # update local copy of X
        for i in range(size):
            x[ i ] = X[ i ]

        x_k = x - x_old
        err = np.dot(x_k, x_k)

        if err < 1e-9:
            converged = True

        # if num == 0 and it_count % 1000 == 0:
        #     print("Iteration %s: converged = %s, err = %s",it_count, converged, err)

    ## end for loop ##


#### end of solve ##############

X = multiprocessing.Array('d', size)

#x = np.zeros_like(b)
#x_new  = np.zeros_like(b)
for i in range(size):
    X[ i ] = 0.

thread_list = []
for t in range(nThreads):
    thread = multiprocessing.Process(target=solve, args=(t,X))
    thread.start()
    thread_list.append(thread)

for thread in thread_list:
    thread.join()

print("time loop: %s sec",time.time()-start,flush=True)
pylab.plot(h, X, 'co') # sin(x)/x
pylab.plot(h, u, color='red') # sin(x)/x
pylab.show() # show the plot

#print("Solution:")
#print(x)
#error = np.dot(A, x) - b
#print("Error:")
#print(error)
