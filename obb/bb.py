# Main file for overlapping branch and bound
#
# obb: Main Function
#
# Principal arguments are:
# f - function to optimize
# g - function gradient
# H - function Hessian
# bndH - Hessian bounding function
# bndT - derivative tensor bounding function
# l - function lower bound
# u - function upper bound
# alg - algorithm type (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
# mod - model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic,
#                  gc - gershgorin cubic)
# A - inequality constraint matrix
# b - inequality constraint upper bound
# E - equality constraint matrix
# d - equality constraint upper bound
#
# Optional arguments are:
# tol - tolerance
# heur - heuristic lattice (0 - off, 1 - on)
# toltype - tolerance type (r - relative, a - absolute)
# vis - visualisation (0 - off, 1 - on)
# qpsolver - QP solver (quadprog, cvxopt)
#
def obb(f, g, H, bndH, bndT, l, u, alg, mod, A=None, b=None, E=None, d=None, tol=1e-2, heur=0, toltype='r', vis=0, qpsolver='cvxopt', countf=1, countsp=1):

    # MPI
    from mpi4py import MPI

    # MPI comm
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Version number
    from version import __version__

    # Circle class
    from circle import circle

    # Algorithm
    if(alg == 'T1'):
        from T1 import runpar
    elif(alg == 'T2_individual'):
        from T2_individual import runpar
    elif(alg == 'T2_synchronised'):
        from T2_synchronised import runpar
    elif(alg == 'T2_synchronised_rr'):
        from T2_synchronised_rr import runpar
    else:
        raise RuntimeError('Algorithm must be T1, T2_individual, T2_synchronised or T2_synchronised_rr.')

    # Model type
    if(mod in ['g','Hz','lbH','E0','Ediag','q']):
        if(alg == 'T2_synchronised_rr'):
            from lboundme_loc_rr import bound
        else:
            from lboundme_loc import bound
    elif(mod in ['gc','c']):
        if(alg == 'T2_synchronised_rr'):
            from lboundgc_loc_rr import bound
        else:
            from lboundgc_loc import bound
    else:
        raise RuntimeError('Model must be Norm Quadratic (q), Norm Cubic (c), Min Eigenvalue Type Quadratic (g, Hz, lbH, E0, Ediag) or Gershgorin Cubic (gc).')

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
        print 'l:', l, '\nu:', u
        if(A is not None): print 'A:', A, '\nb:', b
        if(E is not None): print 'E:', E, '\nd:', d
        print('Model Type: %s') % mod
        print('Number of Processes: %i') % numprocs

        # Run TRS Overlap code
        print('\nStarting Optimization...')

    # All processes run selected algorithm
    xs, fxs, tol, itr = runpar(fw, g, H, (bndH,mod), (bndT,mod), l, u, boundw, circle, A=A, b=b, E=E, d=d, Tol=tol, Heur=heur, TolType=toltype, Vis=vis, qpsolver=qpsolver)

    if(countf == 1):
        print('Processor %i Number of function evaluations: %d ') % (rank,fw.calls)

    if(countsp == 1):
        print('Processor %i Number of subproblems solved: %d ') % (rank,boundw.calls)

    return xs, fxs, tol, itr

# obb_rbf: Main Function with RBF Layer
#
# Principal arguments are:
# f - function to approximate by RBF and optimize
# pts - points at which to sample function (to construct RBF)
# l - function lower bound
# u - function upper bound
# alg - algorithm type (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
# mod - model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic,
#                  gc - gershgorin cubic)
# A - constraint matrix
# lc - constraint lower bound
# uc - constraint upper bounds
#
# Optional arguments are:
# tol - tolerance
# heur - heuristic lattice (0 - off, 1 - on)
# toltype - tolerance type (r - relative, a - absolute)
# vis - visualisation (0 - off, 1 - on)
# qpsolver - QP solver (quadprog, cvxopt)
#
def obb_rbf(f, pts, l, u, alg, mod,  A=None, b=None, E=None, d=None, tol=1e-2, heur=0, toltype='r', vis=0, qpsolver='cvxopt', countf=1, countsp=1):

    # MPI
    from mpi4py import MPI

    # MPI comm
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Version number
    from version import __version__

    # Circle class
    from circle import circle

    # Algorithm
    if(alg == 'T1'):
        from T1 import runpar
    elif(alg == 'T2_individual'):
        from T2_individual import runpar
    elif(alg == 'T2_synchronised'):
        from T2_synchronised import runpar
    elif(alg == 'T2_synchronised_rr'):
        from T2_synchronised_rr import runpar
    else:
        raise RuntimeError('Algorithm must be T1, T2_individual, T2_synchronised or T2_synchronised_rr.')

    # Model type
    if(mod in ['g','Hz','lbH','E0','Ediag','q']):
        if(alg == 'T2_synchronised_rr'):
            from lboundme_loc_rr import bound
        else:
            from lboundme_loc import bound
    elif(mod in ['gc','c']):
        if(alg == 'T2_synchronised_rr'):
            from lboundgc_loc_rr import bound
        else:
            from lboundgc_loc import bound
    else:
        raise RuntimeError('Model must be Norm Quadratic (q), Norm Cubic (c), Min Eigenvalue Type Quadratic (g, Hz, lbH, E0, Ediag) or Gershgorin Cubic (gc).')

    # RBF Layer
    # Fit RBF surrogate
    from fit_rbf import fit_rbf
    fit_rbf(f,pts)

    # Set up RBF function, gradient and Hessian
    from rbf_funcs import frdl, grdl, hrdl

    # Import RBF Bounding Methods
    from rbf_bounds import getkgc as bndH
    from rbf_bounds import getkhc as bndT

    # Define decorator for counting function calls
    def count_calls(fn):
        def wrapper(*args, **kwargs):
            wrapper.calls += 1
            return fn(*args, **kwargs)
        wrapper.calls = 0
        wrapper.__name__= fn.__name__
        return wrapper

    # Optionally count RBF evals
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

    # Master Process
    if(rank == 0):

        # Imports
        from time import localtime, strftime
        from numpy import size

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
        print('Using RBF Layer, number of samples: %d') % size(pts,0)
        print('Dimension: %d \nObjective Function: %s ') % (len(l), f.__name__)
        print 'l:', l, '\nu:', u
        if(A is not None): print 'A:', A, '\nb:', b
        if(E is not None): print 'E:', E, '\nd:', d
        print('Model Type: %s') % mod
        print('Number of Processes: %i') % numprocs

        # Run TRS Overlap code
        print('\nStarting Optimization...')

    # All processes run selected algorithm
    xs, fxs, tol, itr = runpar(frdlw, grdl, hrdl, (bndH,mod), (bndT,mod), l, u, boundw, circle, A=A, b=b, E=E, d=d, Tol=tol, Heur=heur, TolType=toltype, Vis=vis, qpsolver=qpsolver)

    if(countf == 1):
        print('Processor %i Number of RBF surrogate evaluations: %d ') % (rank,frdlw.calls)

    if(countsp == 1):
        print('Processor %i Number of subproblems solved: %d ') % (rank,boundw.calls)

    return xs, fxs, tol, itr

# obb_rbf_coconut: Main Function with RBF Layer on COCONUT test set
#
# Principal arguments are:
# f - RBF approximation from COCONUT test set to optimize
# alg - algorithm type (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
# mod - model type (q - norm quadratic, c - norm cubic, g/Hz/lbH/E0/Ediag -  min eig. quadratic,
#                  gc - gershgorin cubic)
#
# Optional arguments are:
# tol - tolerance
# heur - heuristic lattice (0 - off, 1 - on)
# toltype - tolerance type (r - relative, a - absolute)
# vis - visualisation (0 - off, 1 - on)
# qpsolver - QP solver (quadprog, cvxopt)
#
def obb_rbf_coconut(f, alg, mod, tol=1e-2, heur=0, toltype='r', vis=0, qpsolver='cvxopt', countf=1, countsp=1):

    # MPI
    from mpi4py import MPI

    # MPI comm
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Version number
    from version import __version__

    # Circle class
    from circle import circle

    # Algorithm
    if(alg == 'T1'):
        from T1 import runpar
    elif(alg == 'T2_individual'):
        from T2_individual import runpar
    elif(alg == 'T2_synchronised'):
        from T2_synchronised import runpar
    elif(alg == 'T2_synchronised_rr'):
        from T2_synchronised_rr import runpar
    else:
        raise RuntimeError('Algorithm must be T1, T2_individual, T2_synchronised or T2_synchronised_rr.')

    # Model type
    if(mod in ['g','Hz','lbH','E0','Ediag','q']):
        if(alg == 'T2_synchronised_rr'):
            from lboundme_loc_rr import bound
        else:
            from lboundme_loc import bound
    elif(mod in ['gc','c']):
        if(alg == 'T2_synchronised_rr'):
            from lboundgc_loc_rr import bound
        else:
            from lboundgc_loc import bound
    else:
        raise RuntimeError('Model must be Norm Quadratic (q), Norm Cubic (c), Min Eigenvalue Type Quadratic (g, Hz, lbH, E0, Ediag) or Gershgorin Cubic (gc).')

    # RBF Layer
    # Load RBF surrogate
    from numpy import load, hstack, vstack, delete
    from pkg_resources import resource_stream
    pth = resource_stream('obb','coconut/'+str(f))

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
    A = data['A']
    b = data['b']
    E = data['E']
    d = data['d']

    # Set up relevant global variables
    import config
    config.D = D
    config.N = N
    config.a = a
    config.tl = tl
    config.tl2 = tl2
    config.x = x

    # Load tolerance for 12hrs on serial if requested
    if(tol == '12hr'):
        tol = load(resource_stream('obb','coconut_tol'))[str(f)]
        toltype = 'a'

    # Set up RBF function, gradient and Hessian
    from rbf_funcs import frdl, grdl, hrdl

    # Import RBF Bounding Methods
    from rbf_bounds import getkgc as bndH
    from rbf_bounds import getkhc as bndT

    # Define decorator for counting function calls
    def count_calls(fn):
        def wrapper(*args, **kwargs):
            wrapper.calls += 1
            return fn(*args, **kwargs)
        wrapper.calls = 0
        wrapper.__name__= fn.__name__
        return wrapper

    # Optionally count RBF evals
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

    # Master Process
    if(rank == 0):

        # Imports
        from time import localtime, strftime
        from numpy import size

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
        print('Using RBF Layer, number of samples: %d') % N
        print('Dimension: %d \nObjective Function: %s ') % (D, str(f))
        print 'l:', l, '\nu:', u
        if(A.shape[0] != 0): print 'A:', A, '\nb:', b
        if(E.shape[0] != 0): print 'E:', E, '\nd:', d
        print('Model Type: %s') % mod
        print('Number of Processes: %i') % numprocs

        # Run TRS Overlap code
        print('\nStarting Optimization...')

    # All processes run selected algorithm
    xs, fxs, tol, itr = runpar(frdlw, grdl, hrdl, (bndH,mod), (bndT,mod), l, u, boundw, circle, A=A, b=b, E=E, d=d, Tol=tol, Heur=heur, TolType=toltype, Vis=vis, qpsolver=qpsolver)

    if(countf == 1):
        print('Processor %i Number of RBF surrogate evaluations: %d ') % (rank,frdlw.calls)

    if(countsp == 1):
        print('Processor %i Number of subproblems solved: %d ') % (rank,boundw.calls)

    return xs, fxs, tol, itr
