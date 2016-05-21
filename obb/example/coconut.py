# Example COCONUT RBF code for oBB
#from obb import obb_rbf_coconut

#import from the file project
from bb import obb_rbf_coconut

# Input Settings
# Algorithm (T1, T2_individual, T2_synchronised)
alg = 'T2_synchronised'

# Model type (q - norm quadratic, g/Hz/lbH/E0/Ediag - min eig. quadratic,
# c - norm cubic, gc - gershgorin cubic)
mod = 'c'

# Tolerance (can also be '12hr')
#tol = 1e-2
tol = '12hr'

# Tolerance type (r - relative, a - absolute)
toltype = 'r'

# Visualisation !!Requires matplotlib!! (0 - off, 1 - on)
vis = 1

# QP solver (cvxopt, quadprog)
qpsolver = 'cvxopt'

# Choose RBF approximation from COCONUT test
#f = 'hs045'
#f = 'biggsc4'
#f = 'biggs5'
#f = 'ex2_1_2'
#f = 'bt3'
#f = 'lsnnodoc'
#f = 'hs055'
#f = 'brownden'
#f = 'hatflda'

f = 'expfitb'
#f = 'hs086'

# Run oBB
xs, fxs, tol, itr = obb_rbf_coconut(f, alg, mod, tol=tol, toltype=toltype, vis=vis, qpsolver=qpsolver)
