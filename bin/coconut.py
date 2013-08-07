# Example COCONUT RBF code for oBB
from obb import obb_rbf_coconut

# Input Settings
# Algorithm (T1, T2_individual, T2_synchronised)
alg = 'T1'

# Model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag - min eig. quadratic, 
# gc - gershgorin cubic)
mod = 'c'

# Tolerance (positive number e.g. 1e-2 or '12hr')
tol = '12hr'

# Heuristic lattice (0 - off, 1 - on)
heur = 0

# Tolerance type (r - relative, a - absolute)
toltype = 'a'

# Visualisation !!Requires matplotlib!! (0 - off, 1 - on)
vis = 0

# Choose RBF approximation from COCONUT test
f = 'hs041'

# Run oBB
xs, fxs, tol, itr = obb_rbf_coconut(f, alg, mod, tol=tol, heur=heur, toltype=toltype, vis=vis)
