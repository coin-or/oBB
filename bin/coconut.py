# Example COCONUT RBF code for oBB
from obb import obb_rbf_coconut

# Input Settings
# Algorithm (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
alg = 'T2_synchronised_rr'

# Model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag - min eig. quadratic,
# gc - gershgorin cubic)
mod = 'c'

# Tolerance (can also be '12hr')
tol = 1e-2

# Tolerance type (r - relative, a - absolute)
toltype = 'r'

# Visualisation !!Requires matplotlib!! (0 - off, 1 - on)
vis = 0

# QP solver (cvxopt, quadprog)
qpsolver = 'cvxopt'

# Choose RBF approximation from COCONUT test
f = 'hs041'

# Run oBB
xs, fxs, tol, itr = obb_rbf_coconut(f, alg, mod, tol=tol, toltype=toltype, vis=vis, qpsolver=qpsolver)
