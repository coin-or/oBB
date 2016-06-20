User Guide
==========
This section describes the main interface to oBB and how to use it.

Global Optimization
-------------------
oBB is designed to solve the global optimization problem

  .. math::

    &\min_{x \in \mathbb{R}^n} f(x) \\
    \text{s.t. } \; &l \le x \le u \\
    \text{and } \;  &Ax \le b \\
		    &Ex = d

where the objective function :math:`f` has Lipschitz continuous gradient or Hessian. oBB does not need to know the Lipschitz constants explicitly, it merely requires the user to supply elementwise bounds on the Hessian or derivative tensor of the objective function (see `How to use oBB`_). The linear inequality constraints :math:`Ax \le b` and equality constraints :math:`Ex = d` are optional but the bound constraints :math:`l \le x \le u` are required.

oBB uses local first or second order Taylor type approximations over balls within a parallel branch and bound framework. As with all branch and bound algorithms, the curse of dimensionality limits its use to low dimensional problems. The choice of whether to use first or second order approximations is down to the user  (see `How to use oBB`_). 

For an in-depth technical description of the algorithm see the tech-report [CFG2013]_ and the paper [FGF2013]_.

How to use oBB
--------------
oBB requires the user to write a python script file that defines the functions and parameters necessary to solve the global optimization problem and then passes them to the main **obb** function (see `Example of Use`_). The necessary functions are determined by the choice of a first or second order approximation model. The following approximation models can be passed to oBB's **obb** function using the **mod** argument (see [CFG2013]_ for details):

* First order models: **'q'** - norm based,  **'g'**, **'Hz'**, **'lbH'**, **'E0'**, **'Ediag'** - minimum eigenvalue based
* Second order models: **'c'** - norm based, **'gc'** - minimum eigenvalue based

If using a first order model, the user is required to write the following functions:

* **f(x)** - returns objective function :math:`f` at point :math:`x` (scalar)
* **g(x)** - returns gradient :math:`\nabla_x f` of objective function at :math:`x` (numpy 1D-array)

along with the bounding function:

* **bndH(l,u)** - returns two numpy 2D-arrays of elementwise lower and upper bounds on the Hessian :math:`\nabla_{xx} f` of the objective function over :math:`[l,u]`

For a second order model the user is additionally required to write the function:

* **H(x)** - returns Hessian :math:`\nabla_{xx} f` of objective function at :math:`x` (numpy 2D-array)

and rather than **bndH** the bounding function:

* **bndT(l,u)** - returns two numpy 3D-arrays of elementwise lower and upper bounds on the derivative tensor :math:`\nabla_{xxx} f` of the objective function over :math:`[l,u]`

The type of parallel branch and bound algorithm to use should be passed to oBB's **obb** function using the **alg** argument and can be one of the following (see [CFG2013]_ for details):

* **'T1'** - bounds in parallel
* **'T2_individual'**, **'T2_synchronised'** - tree in parallel
* **'T2_synchronised_rr'** - tree in parallel *with range reduction*

See `Example of Use`_ for an in-depth worked example in python.

Optional Arguments
------------------
oBB allows the user to specify several optional arguments that control the behaviour of the algorithm:

* **tol** - objective function tolerance (e.g. **1e-2**, the default)
* **toltype** - tolerance type (**'r'** - relative [default], **'a'** - absolute)
* **countf** - count objective function evaluations (**0** - off, **1** - on [default])
* **countsp** - count subproblem evaluations (**0** - off, **1** - on [default])

and if `matplotlib <http://www.matplotlib.org/>`_ is installed:

* **vis** - visualisation of the algorithm in 2D (**0** - off [default], **1** - on)

Note that the inequality constraint arguments :math:`A, b` and  equality constraint arguments :math:`E, d` are also optional as the optimization problem is only required to be bound constrained.

The user can also specify the QP solver that the algorithm calls to obtain feasible upper bounds (see [CFG2013]_ for details). At present the user can choose from `CVXOPT's qp solver <http://cvxopt.org/>`_ or the `QuadProg++ solver <http://github.com/mpy/PyQuadProg/>`_ using the optional argument:    

* **qpsolver** - QP solver to use (**'cvxopt'** - CVXOPT's qp [default], **'quadprog'** - QuadProg++)

Note that the QuadProg++ solver is faster as it is written in C++ but has very limited error handling and may not work in all cases. The CVXOPT solver is slower as it is written in Python but considerably more stable.

Example of Use
--------------
Suppose we wish to solve the following global optimization problem:

  .. math::

    &\min_{x \in \mathbb{R}^n} \sum_{i=1}^n \sin(x_i) \\
    \text{s.t. } \; &-1 \le x_i \le 1 \; \; \forall i=1,\dotsc,n \\
    \text{and }  \; &\sum_{i=1}^n -x_i \le 1

One can see that the gradient :math:`g`, Hessian :math:`H` and third order derivative tensor :math:`T` are given by

  .. math::
  
    &g(x) = ( \cos(x_1), \dotsc, \cos(x_n))^T \\
    &H(x) = diag( -\sin(x_1), \dotsc, -\sin(x_n)) \\
    &T(x) = diagt( -\cos(x_1), \dotsc, -\cos(x_n))

where :math:`diagt` is the tensor diagonal function (i.e. :math:`diagt(v)` places the vector :math:`v` on the diagonal of the tensor).

It is straightforward to obtain elementwise bounds on the Hessian matrix :math:`H` and third order derivative tensor :math:`T` as both :math:`sin` and :math:`cos` can be bounded below and above by -1 and 1 respectively. 

We can code this up in a python script file, let's call it sines.py as follows:  

  .. code-block:: python
  
      # Example code for oBB
      from obb import obb
      from numpy import sin, cos, diag, ones, zeros

      # Input Settings
      # Algorithm (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
      alg = 'T1'

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

      # Set up sum of sines test function
      # Dimension
      D = 2
      # Constraints
      l = -1*ones(D)
      u = 1*ones(D)
      A = -1*ones((1,D))
      b = 1
      # Required functions
      f = lambda x: sum(sin(x))
      g = lambda x: cos(x)
      H = lambda x: diag(-sin(x))
      bndH = lambda l,u: (diag(-ones(D)), diag(ones(D)))
      bndT = lambda l,u: (diagt(-ones(D)), diagt(ones(D)))

      # Name objective function
      f.__name__ = 'Sum of Sins'

      # Run oBB
      xs, fxs, tol, itr = obb(f, g, H, bndH, bndT, l, u, alg, mod, A=A, b=b, tol=tol)
      
This file is included in oBB as sines.py, to run it see `Running the Algorithm`_.

Running the Algorithm
---------------------
To run the user-created python script file (e.g. sines.py, see `Example of Use`_) we need to execute it using MPI's mpiexec command, specifying the number of processor cores with the -n option. For example, to run oBB on four processor cores we simply execute the following shell command:

  .. code-block:: bash

     $ mpiexec -n 4 python sines.py

Note that if using the MPICH implementation of MPI we first need to start an mpd daemon in the background:

 .. code-block:: bash

    $ mpd &

but this is not necessary for other MPI implementations, e.g. OpenMPI.

Using the RBF Layer
-------------------
oBB can optionally approximate the objective function :math:`f` by a radial basis function (RBF) surrogate and optimize the approximation instead (see [FGF2013]_ for details). The advantage of this approach is that the user merely needs to supply the objective function and a set of points at which it should be evaluated to construct the RBF approximation. The disadvantage is that the optimum found by the algorithm will only be close to the optimum of the objective function if it is sampled at sufficiently many points.

As before, the user is required to write a python script file that defines the functions and parameters necessary to solve the problem and then passes them to the **obb_rbf** function. In addition to the approximation model, algorithm type and objective function arguments described in `How to use oBB`_ only an :math:`n` by :math:`m` numpy array of :math:`m` points at which to sample the objective function needs to be passed to the **obb_rbf** function using the **pts** argument. 

For example, suppose we wish to solve an RBF approximation to the problem given in the `Example of Use`_ section:

  .. math::

    &\min_{x \in \mathbb{R}^n} \sum_{i=1}^n \sin(x_i) \\
    \text{s.t. } \; &-1 \le x_i \le 1 \; \; \forall i=1,\dotsc,n \\
    \text{and }  \; &\sum_{i=1}^n -x_i \le 1

We can code this up in a python script file, let's call it sines_rbf.py as follows:  

  .. code-block:: python
  
	# Example RBF Layer code for oBB
	from obb import obb_rbf
	from numpy import sin, ones
	from numpy.random import rand, seed

	# Input Settings
	# Algorithm (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
	alg = 'T1'

	# Model type (q - norm quadratic, g/Hz/lbH/E0/Ediag - min eig. quadratic, 
	# c - norm cubic, gc - gershgorin cubic)
	mod = 'c'

	# Tolerance
	tol = 1e-2

	# Set up sum of sines test function
	# Dimension
	D = 2
	# Constraints
	l = -1*ones(D)
	u = 1*ones(D)
	A = -1*ones((1,D))
	b = 1
	# Required functions
	f = lambda x: sum(sin(x))

	# Generate 10*D sample points for RBF approximation
	seed(5) # !!Sample points have to be the same on all processors!! 
	pts = rand(10*D, D)

	# Scale points so they lie in [l,u]
	for i in range(0,D):
	    pts[:,i] = l[i] + (u[i]-l[i])*pts[:,i]

	# Name objective function
	f.__name__ = 'RBF Sum of Sins'

	# Run oBB
	xs, fxs, tol, itr = obb_rbf(f, pts, l, u, alg, mod, A=A, b=b, tol=tol)
		
Note the use of **obb_rbf** instead of **obb** and the need for a random number seed so that the sample points are the same on all processors. This file is included in oBB as sines_rbf.py, to run it see `Running the Algorithm`_. 

RBF Layer for the COCONUT Test Set
----------------------------------
oBB comes with a set of pre-computed RBF approximations to selected functions from the `COCONUT test set <http://www.mat.univie.ac.at/~neum/glopt/coconut/Benchmark/Benchmark.html>`_ that were used to produce the numerical results in the paper [CFG2013]_. In order to optimize these approximations using oBB, the user is required to write a python script file that defines the desired function and tolerance and then passes them to the **obb_rbf_coconut** function (see [CFG2013]_ for a list of all 31 functions available). For example, to optimize an RBF approximation to the **'hs041'** function the user could write the following python script file, let's call it coconut.py: 

    .. code-block:: python
  
	# Example COCONUT RBF code for oBB
	from obb import obb_rbf_coconut

	# Input Settings
	# Algorithm (T1, T2_individual, T2_synchronised, T2_synchronised_rr)
	alg = 'T1'

	# Model type (q - norm quadratic, g/Hz/lbH/E0/Ediag - min eig. quadratic, 
	# c - norm cubic, gc - gershgorin cubic)	
	mod = 'c'

	# Tolerance (can also be '12hr')
	tol = 1e-2

	# Choose RBF approximation from COCONUT test
	f = 'hs041'

	# Run oBB
	xs, fxs, tol, itr = obb_rbf_coconut(f, alg, mod, tol=tol)
  
Note the use of **obb_rbf_coconut** as the calling function and the optional **'12hr'** tolerance setting which runs the algorithm to the absolute tolerance obtained by a serial code in twelve hours (see [CFG2013]_ for details). This file is included in oBB as coconut.py, to run it see `Running the Algorithm`_. 

Acknowledgements
----------------
This work was supported by EPSRC grants `EP/I028854/1 <http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/I028854/1>`_ (PI: `Dr Coralia Cartis <http://www.maths.ox.ac.uk/people/profiles/coralia.cartis>`_) and NAIS `EP/G036136/1 <http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/G036136/1>`_.
We are also grateful to `Prof Nick Gould <http://www.numerical.rl.ac.uk/people/nimg/>`_ for his help and advice during development, Mehdi Towhidi for providing us with his `PyQuadProg code <http://github.com/mpy/PyQuadProg>`_ and Alberto Guida for his range reduction strategy enhancement.

References
----------

.. [CFG2013]   
   Cartis, C., Fowkes, J. M. and Gould, N. I. M. (2013) 'Branching and Bounding Improvements for Global Optimization Algorithms with Lipschitz Continuity Properties', *ERGO Technical Report*, no. 13-010, pp. 1-33. http://www.maths.ed.ac.uk/ERGO/pubs/ERGO-13-010.html

.. [FGF2013]   
   Fowkes, J. M. , Gould,  N. I. M. and Farmer, C. L. (2013) 'A Branch and Bound Algorithm for the Global Optimization of Hessian Lipschitz Continuous Functions', *Journal of Global Optimization*, vol. 56, no. 4, pp. 1791-1815. http://dx.doi.org/10.1007/s10898-012-9937-9 
