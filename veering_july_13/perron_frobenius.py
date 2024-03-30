
import random
import math
import copy
try:
    from cdecimal import *
except ImportError:
    from decimal import *

# to install cdecimal put cdecimal.pyd in SnapPy directory
# also, take numbers.pyd from python26 and put it into Snappy directory
  
#getcontext().prec = 200

def vector_norm(vec):
    return sum(x*x for x in vec).sqrt()

""" 
from wikipedia:
5. The only eigenvectors whose components are all positive are those associated with the eigenvalue r.
this should guarantee if algorithm converges then it converges to correct eigenvector

We use the trick that matrices M and M+I have the same eigenvectors, but
the real part of the eigenvalues of M+I have increased by 1 
"""

def mat_mult(m1, m2):
    return [[sum(m1[i][k]*m2[k][j] for k in range(len(m1[i]))) 
            for j in range(len(m1))] for i in range(len(m1))]

def perron_eigenvector(mat_o, tol=Decimal('10e-30'), max_iters = 10000):
    getcontext().prec = 100
    mat = copy.deepcopy(mat_o)
    
    mat = [map(Decimal, row) for row in mat]
    
    for i in range(len(mat)):
        mat[i][i] += Decimal(1)
        
    for i in range(4):
        mat = mat_mult(mat,mat)
		
    prev_vec = [Decimal("%.15g" % random.random()) for i in range(len(mat))]
    next_vec = prev_vec
    
    found = False
    for iter in xrange(max_iters):
        nv_norm_recip = Decimal(1) / vector_norm(next_vec)
        next_vec, prev_vec = [sum(row[i]*next_vec[i]*nv_norm_recip
                                 for i in range(len(mat))) for row in mat], next_vec
        if vector_norm([next_vec[i] - prev_vec[i] for i in range(len(prev_vec))]) < tol:
            found = True
            break

    if not found:
        print 'failed to converge after:', max_iters, 'iterations', next_vec
    # reduce the eigenvalue by 1
    vn = vector_norm(next_vec)
    scale = ((vn.ln()/Decimal(16)).exp() - Decimal(1))/vn
    next_vec = [x * scale for x in next_vec]
    return next_vec

def normalize(vec):
    norm = vector_norm(vec)
    return [x/norm for x in vec]
