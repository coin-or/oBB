from __future__ import division
# RBF Code for Parallel TRS Overlap Paper
# First argument is algorithm to run e.g. run_T2
# Second argument is function to run e.g. ex2_1_1
# Third argument is model type: 
# q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic, 
# gc - gershgorin cubic
# Fourth (optional) argument is tolerance e.g. 1e-2
# Fifth (optional) argument is heuristic alg. i.e. (0 - off, 1 - on)
# Sixth (optional) argument is tolerance type i.e. (r - relative, a - absolute)
# !! Before running, check all algorithms for correct stopping rules !!

# Count func evals (1 - yes, 0 - no)
countf = 1

# Count subproblems solved (1 - yes, 0 - no)
countsp = 1

# MPI
from mpi4py import MPI

# MPI comm
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Default Tolerance
if(len(sys.argv) < 5):
	tol = 1e-2
	heur = 0
	toltype = 'r'
elif(len(sys.argv) < 6):
	tol = float(sys.argv[4])
	heur = 0
	toltype = 'r'
elif(len(sys.argv) < 7):
	tol = float(sys.argv[4])
	heur = int(sys.argv[5])
	toltype = 'r'		
else:
	tol = float(sys.argv[4])
	heur = int(sys.argv[5])
	toltype = str(sys.argv[6])

# Model type
regmod = sys.argv[3]
if(regmod in ['g','Hz','lbH','E0','Ediag','q']):
	istr = 'from lboundme_loc import bound'
elif(regmod in ['gc','c']):
	istr = 'from lboundgc_loc import bound'
else:
	raise RuntimeError('Model must be Norm Quadratic (q), Norm Cubic (c), Min Eigenvalue Type Quadratic (g, Hz, lbH, E0, Ediag) or Gershgorin Cubic (gc).')	
exec istr

# Bounding Method
from getkgc import getkgc as bndH
from getkhc import getkhc as bndT

# Algorithm
istr = 'from '+sys.argv[1]+' import runpar'
exec istr

# Set up function, gradient and Hessian
from rdl_funcs import frdl, grdl, hrdl

# Data loading
from numpy import load

# List of functions to test
fcn = sys.argv[2]
pth = '/exports/work/maths_oro/jfowkes/Coconut30/'+fcn

# Load data
data = load(pth)
N = data['N']
D = data['D']
a = data['a']
x = data['x']
l = data['l']
u = data['u']
tl = data['tl']
tl2 = data['tl2']
Ac = data['Ac']
lc = data['lc']
uc = data['uc']
if(Ac.shape[0] == 0):
	Ac = None
	lc = None
	uc = None

# Set up relevant global variables
import config 
config.D = D
config.N = N
config.a = a
config.tl = tl
config.tl2 = tl2
config.x = x

# Define decorator for counting function calls
def count_calls(fn):
	def wrapper(*args, **kwargs):
		wrapper.calls += 1
		return fn(*args, **kwargs)
	wrapper.calls = 0
	wrapper.__name__= fn.__name__
	return wrapper

# Optionally count fevals
if(countf == 1):

	@count_calls
	def frdlw(x):
		return frdl(x)
else:
	frdlw = frdl

# Optionally count subproblems solved
if(countsp == 1):
	
	@count_calls
	def boundw(*a,**kw):
		return bound(*a,**kw)
else:
	boundw = bound

# Set up circle class (needs to be in root file)
class circle:
	lbound = 0
	np = None
	name = None
	parents = None
	def __init__(self, xcc, rcc):
		self.xc = xcc
		self.r = rcc

# Master Process
if(rank == 0):

	# Imports
	from time import localtime, strftime
	
	# Get time
	t = localtime()

	# Get number of processors
	numprocs = comm.Get_size()	

	# Output basic info
	print('****************************************************')
	print('* Trust Region Optimization using Branch and Bound *')
	print('* Time: %02d:%02d:%02d  Date: %s   Ver: 0.2 *') % (t[3], t[4], t[5], strftime("%d %B %Y",t))
	print('****************************************************\n')

	# Output problem details
	print('Problem Details: \n----------------')
	print('Dimension: %d \nObjective Function: %s ') % (config.D, fcn)
	print('l: '),
	print l
	print('u: '),
	print u 
	print('Model Type: %s') % regmod
	print('Number of Processes: %i') % numprocs

	# Run TRS Overlap code
	print('\nStarting Optimization...')
	xs, fxs, tol, itr = runpar(frdlw, grdl, hrdl, (bndH,regmod), (bndT,regmod), l, u, boundw, circle, Ac, lc, uc, Tol=tol, Heur=heur, TolType=toltype)

# Worker processes
else:
            
	_, _, _, _ = runpar(frdlw, grdl, hrdl, (bndH,regmod), (bndT,regmod), l, u, boundw, circle, Ac, lc, uc, Tol=tol, Heur=heur, TolType=toltype)
      
if(countf == 1):
	print('Processor %i Number of function evaluations: %d ') % (rank,frdlw.calls)

if(countsp == 1):
	print('Processor %i Number of subproblems solved: %d ') % (rank,boundw.calls)

