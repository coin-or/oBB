# Main file for overlapping branch and bound
#
# Principal arguments are:
# f - function to optimize
# g - function gradient
# H - function Hessian
# bndH - Hessian bounding function
# bndT - derivative tensor bounding function
# l - function lower bound
# u - function upper bound
# alg - algorithm type (T1, T2_individual, T2_synchronised)
# mod - model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic, 
# 		   gc - gershgorin cubic)
# A - constraint matrix
# lc - constraint lower bound 
# uc - constraint upper bounds 
#
# Optional arguments are:
# tol - tolerance 
# heur - heuristic lattice (0 - off, 1 - on)
# toltype - tolerance type (r - relative, a - absolute)
# vis - visualisation (0 - off, 1 - on)
#
def obb(f, g, H, bndH, bndT, l, u, alg, mod, A=None, lc=None, uc=None, tol=1e-2, heur=0, toltype='r', vis=0, countf=1, countsp=1):

	# MPI
	from mpi4py import MPI

	# MPI comm
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	# Version number
	from version import __version__

	# Circle class 
	from circle import circle

	# Model type
	if(mod in ['g','Hz','lbH','E0','Ediag','q']):
		from lboundme_loc import bound
	elif(mod in ['gc','c']):
		from lboundgc_loc import bound
	else:
		raise RuntimeError('Model must be Norm Quadratic (q), Norm Cubic (c), Min Eigenvalue Type Quadratic (g, Hz, lbH, E0, Ediag) or Gershgorin Cubic (gc).')	

	# Algorithm
	if(alg == 'T1'):
		from T1 import runpar
	elif(alg == 'T2_individual'):
		from T2_individual import runpar
	elif(alg == 'T2_synchronised'):	
		from T2_synchronised import runpar
	else:
		raise RuntimeError('Algorithm must be T1, T2_individual or T2_synchronised.')	
				
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
		def fw(x):
			return f(x)
	else:
		fw = f

	# Optionally count subproblems solved
	if(countsp == 1):
	
		@count_calls
		def boundw(*a,**kw):
			return bound(*a,**kw)
	else:
		boundw = bound

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
		print('* Time: %02d:%02d:%02d  Date: %s   Ver: %s *') % (t[3], t[4], t[5], strftime("%d %B %Y",t),__version__)
		print('****************************************************\n')

		# Output problem details
		print('Problem Details: \n----------------')
		print('Dimension: %d \nObjective Function: %s ') % (len(l), f.__name__)
		print('l: '),
		print l
		print('u: '),
		print u 
		print('Model Type: %s') % mod
		print('Number of Processes: %i') % numprocs

		# Run TRS Overlap code
		print('\nStarting Optimization...')
		
	# All processes run selected algorithm
	xs, fxs, tol, itr = runpar(fw, g, H, (bndH,mod), (bndT,mod), l, u, boundw, circle, A=A, lc=lc, uc=uc, Tol=tol, Heur=heur, TolType=toltype, Vis=vis)
	  
	if(countf == 1):
		print('Processor %i Number of function evaluations: %d ') % (rank,fw.calls)

	if(countsp == 1):
		print('Processor %i Number of subproblems solved: %d ') % (rank,boundw.calls)
	
	return xs, fxs, tol, itr	
