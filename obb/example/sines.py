# Example code for oBB
#from obb import obb

#import from the file project
from bb import obb

from numpy import sin, cos, diag, ones, zeros

# Input Settings
# Algorithm (T1, T2_individual, T2_synchronised)
alg = 'T2_synchronised'

# Model type (q - norm quadratic, g/Hz/lbH/E0/Ediag - min eig. quadratic,
# c - norm cubic, gc - gershgorin cubic)
mod = 'c'

# Tolerance
tol = 1e-2

# Tolerance type (r - relative, a - absolute)
toltype = 'r'

# Visualisation !!Requires matplotlib!! (0 - off, 1 - on)
vis = 1

# QP solver (cvxopt, quadprog)
qpsolver = 'cvxopt'

# Tensor diagonal function
def diagt(v):
    T = zeros((D,D,D))
    for i in range(0,D):
        T[i,i,i] = v[i]
    return T

# Set up sum of sines test function
# Dimension
D = 2
# Constraints
# l<=x<=u
l = -1*ones(D)
u = 1*ones(D)

#Ax <= b  ==  Sum (-xi) <= 1
A = -1*ones((1,D))  #A = row vector of -1
b = 1

# Required functions
f = lambda x: sum(sin(x)) #objective function f
g = lambda x: cos(x)  #gradient of f
H = lambda x: diag(-sin(x)) # Hessian
bndH = lambda l,u: (diag(-ones(D)), diag(ones(D))) #elementwise lower and upper bound on the hessian
bndT = lambda l,u: (diagt(-ones(D)), diagt(ones(D))) #elementwise lower and upper bound on the derivative tensor

# Name objective function
f.__name__ = 'Sum of Sins'

# Run oBB
xs, fxs, tol, itr = obb(f, g, H, bndH, bndT, l, u, alg, mod, A=A, b=b, tol=tol, toltype=toltype, vis=vis, qpsolver=qpsolver)
