#import SnapPea
import os
import pprint

""" Henry Segerman sent me this code.
    If change_of_basis_matrix returns (a,c,b,d) then the change of basis matrix is given by
    [(a,b), (c,d)].
    
    Example usage: (new curve, old meridian, old longitude)
    change_of_basis_matrix([0, 8, 2, 0, 8, 2, 4, 16, 4, 2, 2, 4],
                           [1, -1,  1,  0, -2,  0,  1, -2,  0,  0, -1, -1],
                           [ 0,  2,  0,  0,  2,  0, -2,  0,  0,  0,  2,  0])      
"""

def eea(u, v):
    """extended euclidean algorithm: uu1 + vu2 = u3 = gcd(u,v)"""
    u1 = 1
    u2 = 0
    u3 = u
    v1 = 0
    v2 = 1
    v3 = v
    while v3 != 0:
        q = u3 / v3
        t1 = u1 - q * v1
        t2 = u2 - q * v2
        t3 = u3 - q * v3
        u1 = v1
        u2 = v2
        u3 = v3
        v1 = t1
        v2 = t2
        v3 = t3
    if u3 < 0:
        u1, u2 = -u1, -u2
    return u1, u2, u3

def convert_to_NZ(vec):
    """converts a vector in the form of a gluing equation from snappy into the kind of vector needed for the Neumann-Zagier stuff"""
    out1 = []
    out2 = []
    for i in range(len(vec)/3):
        out1.append(vec[3*i] - vec[3*i+2])
        out2.append(-vec[3*i+1] + vec[3*i+2])
    #print vec, '->', out1, out2 
    return out1,out2

def vec_add(v1, v2):
    out = []
    for j in range(len(v1)):
        out.append(v1[j] + v2[j])
    return out

def vec_scalar_mult(a, v):
    out = []
    for j in range(len(v)):
        out.append(a * v[j])
    return out

def dot(v1, v2):
    out = 0
    for i in range(len(v1)):
        out += v1[i]*v2[i]
    return out

def isect(v1, v2):
    """returns signed intersection number between v1 and v2"""
    a1,b1 = convert_to_NZ(v1)
    a2,b2 = convert_to_NZ(v2)
    temp = (-dot(b1,a2) + dot(a1,b2))
    if temp % 2 == 1:
        print 'warning, isect not divisible by 2'
    return (-dot(b1,a2) + dot(a1,b2))/2

def change_of_basis_matrix(new_curve, c1, c2):
    """makes new_curve direction the merid"""
    if isect(new_curve, c1) == 0: #do nothing
        out = 1,0,0,1
    elif isect(new_curve, c2) == 0: #swap so that ladder is first
        out = 0,1,-1,0
    else:
        i1 = isect(new_curve, c1)
        i2 = isect(new_curve, c2)
        # then new_curve is in same direction as i2 * c1 - i1 * c2 (or possibly -this, let's fix if it's wrong)
        u1,u2,u3 = eea(i2,-i1)
        # then u1 * i2 + u2 * -i1 = 1
        out = i2,-u2,-i1,u1
         #so other basis element is given by -u2 * c1 + u1 * c2
        
    basis_elt_2 = vec_add(vec_scalar_mult(out[1], c1), vec_scalar_mult(out[3], c2))
        
    if isect(new_curve, basis_elt_2) > 0:
        return out
    else:
        return -out[0], -out[1], -out[2], -out[3]
   
   
