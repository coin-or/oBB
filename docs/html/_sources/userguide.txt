User Guide
==========
This section describes the main interface to oBB and how to use it.

Global Optimization
-------------------
oBB is designed to solve the global optimization problem

  .. math::

    &\min_{x \in \mathcal{R}^n} f(x) \\
    \text{s.t. } \; &l \le x \le u \\
    \text{and } \; &lc \le Ax \le uc

where the domain :math:`\mathcal{D}` is convex, compact and the objective function :math:`f` has Lipschitz continuous gradient or Hessian. oBB does not need to know the Lipschitz constants explicitly, it merely requires the user to supply elementwise bounds on the Hessian or derivative tensor of the objective function (see `How to use oBB`_).

oBB uses local first or second order Taylor type approximations over balls within a parallel branch and bound framework. As with all branch and bound algorithms, the curse of dimensionality limits its use to low dimensional problems. The choice of whether to use first or second order approximations is down to the user  (see `How to use oBB`_). 

For an in-depth technical description of the algorithm see the tech-report [CFG2013]_ and the paper [FGF2012]_.

How to use oBB
--------------
oBB requires the user to write a python script file which defines the functions and parameters necessary to solve the problem. These are determined by the choice of a first or second order approximation model (see [CFG2013]_ for details):

* First order models: *q* - norm based,  *g*, *Hz*, *lbH*, *E0*, *Ediag* - minimum eigenvalue based
* Second order models: *c* - norm based, *gc* - minimum eigenvalue based

If using a first order model, the following functions are required:

* **f(x)** - returns objective function :math:`f` at point :math:`x` (scalar)
* **g(x)** - returns gradient :math:`\nabla_x f` of objective function at :math:`x` (numpy vector)

and the bounding function:

* **bndH(l,u)** - returns two numpy matrices of elementwise lower and upper bounds on the Hessian :math:`\nabla_{xx} f` of the objective function over :math:`[l,u]`

For a second order model we also require:

* **H(x)** - returns Hessian :math:`\nabla_{xx} f` of objective function at :math:`x` (numpy matrix)

and rather than bndH the bounding function:

* **bndT(l,u)** - returns two third order numpy tensors of elementwise lower and upper bounds on the derivative tensor :math:`\nabla_{xxx} f` of the objective function over :math:`[l,u]`

We now need only specify the type of parallel branch and bound algorithm to use (again, see [CFG2013]_ for details):

* *T1* - bounds in parallel
* *T2_individual*, *T2_synchronised* - tree in parallel

See `Example of Use`_ for an in-depth worked example in python.

Optional Parameters
-------------------
oBB allows the user to specify several optional parameters to control the behaviour of the algorithm:

* *tol* - objective function tolerance (e.g. 1e-2, the default)
* *toltype* - tolerance type (r - relative [default], a - absolute)
* *countf* - count objective function evaluations (0 - off, 1 - on [default])
* *countsp* - count subproblem evaluations (0 - off, 1 - on [default])

and if matplotlib (http://www.matplotlib.org/) is installed:

* vis - visualisation of the algorithm in 2D (0 - off [default], 1 - on)

Note that the linear constraint parameters :math:`A, lc` and  :math:`uc` are also optional as the optimization problem is only required to be bound constrained.

Example of Use
--------------
Suppose we wish to solve the following problem:

  .. math::

    &\min_{x \in \mathcal{R}^n} \sum_{i=1}^n \sin(x_i) \\
    \text{s.t. } \; &-1 \le x_i \le 1 \; \; \forall i=1,\dotsc,n \\
    \text{and }  \; &-1 \le \sum_{i=1}^nx_i \le 1

One can see that the gradient :math:`g`, Hessian :math:`H` and third order derivative tensor :math:`T` are given by

  .. math::
  
    &g(x) = ( \cos(x_1), \dotsc, \cos(x_n))^T \\
    &H(x) = diag( -\sin(x_1), \dotsc, -\sin(x_n)) \\
    &T(x) = diagt( -\cos(x_1), \dotsc, -\cos(x_n))

where :math:`diagt` is the tensor diagonal function (i.e. :math:`diagt(v)` places the vector :math:`v` on the diagonal of the tensor).

It is straightforward to obtain elementwise bounds on the Hessian matrix :math:`H` and third order derivative tensor :math:`T` as both :math:`sin` and :math:`cos` can be bounded below and above by -1 and 1 respectively. 

We can code this up in a python script file, let's call it sins.py as follows:  

  .. code-block:: python
  
    # Example code for oBB
    from obb import obb
    from numpy import sin,cos,diag,ones,zeros

    # Input Settings
    # Algorithm (T1, T2_individual, T2_synchronised)
    alg = 'T2_synchronised'

    # Model type (q - norm quadratic, g/Hz/lbH/E0/Ediag - min eig. quadratic, 
    # c - norm cubic, gc - gershgorin cubic)
    mod = 'c'

    # Tolerance
    tol = 1e-2

    # Tensor diagonal function
    def diagt(v):
	T = zeros((D,D,D))
	for i in range(0,D):
		T[i,i,i] = v[i]
	return T

    # Set up sum of sins test function
    # Dimension
    D = 2 
    # Constraints
    l = -1*ones(D)
    u = 1*ones(D)
    A = ones((1,D))
    lc = -1; uc = 1
    # Required functions
    f = lambda x: sum(sin(x))
    g = lambda x: cos(x)
    H = lambda x: diag(-sin(x))
    bndH = lambda l,u: (diag(-ones(D)), diag(ones(D)))
    bndT = lambda l,u: (diagt(-ones(D)), diagt(ones(D)))

    # Name objective function
    f.__name__ = 'Sum of Sins'

    # Run oBB
    xs, fxs, tol, itr = obb(f, g, H, bndH, bndT, l, u, alg, mod, A=A, lc=lc, uc=uc, tol=tol)

This file is included in oBB as sins.py.

Running the Algorithm
---------------------
To run the user-created python script file (e.g. sins.py) we need to execute it using MPI's mpiexec command, specifying the number of processor cores with the -n option. For example, to run oBB on four processor cores we simply execute the following shell command:

  .. code-block:: bash

     $ mpiexec -n 4 python sins.py

Note that if using the MPICH2 implementation of MPI we first need to start an mpd daemon in the background:

 .. code-block:: bash

    $ mpd &

but this is not necessary for other MPI implmentations, e.g. OpenMPI.

And that's all there is to it!

References
----------

.. [CFG2013]   
   C. Cartis, J. M. Fowkes and N. I. M. Gould. (2013) "Branching and Bounding Improvements for Global Optimization Algorithms with Lipschitz Continuity Properties", ERGO Technical Report, no. 13-010, pp. 1-33. http://www.maths.ed.ac.uk/ERGO/pubs/ERGO-13-010.html

.. [FGF2012]   
   J. M. Fowkes, N. I. M. Gould and C. L. Farmer. (2012) "A Branch and Bound Algorithm for the Global Optimization of Hessian Lipschitz Continuous Functions", Journal of Global Optimization, pp. 1-25. ISSN 0925-5001. http://dx.doi.org/10.1007/s10898-012-9937-9 
