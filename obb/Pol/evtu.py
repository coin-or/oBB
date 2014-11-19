# Code for testing Evtushenko polynomials
# First argument is algorithm to run e.g. T2
# Second argument is function series e.g. 1
# Third argument is function realisation e.g. 2
# Fourth argument is model type:
# q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic,
# gc - gershgorin cubic
# Fifth (optional) argument is tolerance e.g. 1e-2
# Sixth (optional) argument is heuristic alg. i.e. (0 - off, 1 - on)
# Seventh (optional) argument is tolerance type i.e. (r - relative, a - absolute)

# Parse Inputs
alg = sys.argv[1]
ser = int(sys.argv[2])
rea = int(sys.argv[3])
regmod = sys.argv[4]
if(len(sys.argv) < 6):
    tol = 1e-2
    heur = 0
    toltype = 'r'
elif(len(sys.argv) < 7):
    tol = float(sys.argv[5])
    heur = 0
    toltype = 'r'
elif(len(sys.argv) < 8):
    tol = float(sys.argv[5])
    heur = int(sys.argv[6])
    toltype = 'r'
else:
    tol = float(sys.argv[5])
    heur = int(sys.argv[6])
    toltype = str(sys.argv[7])

# Load random seed for chosen function
from numpy import load
data = load('evtustates')
states = data['states']

# Get series and realisation
# Setup series
if(ser == 1):
    dim = 3
    deg = 4
elif(ser == 2):
    dim = 3
    deg = 6
elif(ser == 3):
    dim = 3
    deg = 8
elif(ser == 4):
    dim = 4
    deg = 4
elif(ser == 5):
    dim = 4
    deg = 6
else:
    raise RuntimeError('Series must be 1,2,3,4 or 5.')

# Set random number generator state
from random import setstate
setstate(states[ser-1,rea-1])

# Get bounds
from scipy.misc import comb
from numpy import ones
m = comb(deg-1+dim,dim,exact=1)
l = -m*ones(dim)
u = m*ones(dim)

# Generate Polynomial
from evtupoly import evtupoly
f,g,H,bndH,bndT = evtupoly(dim,deg)

# Name objective function
f.__name__ = 'Evtushenko Random Polynomial, Series '+str(ser)+', Realisation '+str(rea)

# Run oBB
xs, fxs, tol, itr = obb(f, g, H, bndH, bndT, l, u, alg, regmod, tol=tol, heur=heur, toltype=toltype)
