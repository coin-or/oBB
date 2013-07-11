===
OBB
===
oBB is an algorithm for the parallel global optimization of functions with Lipchitz continuous gradient or Hessian.

Documentation
-------------

HTML documentation is available at http://packages.python.org/oBB

Requirements
------------
oBB requires the following software to be installed:

* Python 2.6 to 2.7 (http://www.python.org/)
* A working implementation of MPI-2 (e.g. OpenMPI, http://www.open-mpi.org/)

.. * NAG Fortran Library, tested with Mark 23 FLL6I23DCL  (http://www.nag.co.uk/)

Additionally, the following python packages should be installed (these will be installed automatically if using *pip*, see `Installation using pip (recommended)`_):

* NumPy 1.3.0 or higher (http://www.numpy.org/)
* MPI for Python 1.3 or higher (http://mpi4py.scipy.org/) 
* CVXOPT 1.1.3 or higher (http://cvxopt.org/)

Optionally, the following software may be manually installed for added functionality:

.. * SymPy 0.7.1 or higher (http://www.sympy.org/) - for automatically calculating derivatives

* matplotlib 1.1.0 or higher (http://www.matplotlib.org/) - for visualising the algorithm in 2D

Installation using pip (recommended)
--------------------------------------
For easy installation, use *pip* (http://www.pip-installer.org/)::

    $ [sudo] pip install obb

or alternatively *easy_install* (deprecated)::

    $ [sudo] easy_install obb
    
If you do not have root privileges or you want to install oBB for your private use, you can use::

    $ pip install --user obb
      
which will install oBB in your home directory.

Manual installation
-------------------
Alternatively, you can download the source code and unpack as follows::

    $ wget http://pypi.python.org/packages/source/o/oBB/oBB-X.X.tar.gz
    $ tar -xzvf oBB-X.X.tar.gz
    $ cd oBB-X.X

and then build and install manually using::

    $ python setup.py build
    $ python setup.py install

If you do not have root privileges or you want to install oBB for your private use, you can use::

    $ python setup.py install --user
    
instead.    

Uninstallation
--------------
If oBB was installed using *pip* you can uninstall as follows::

    $ [sudo] pip uninstall obb

If oBB was installed manually you have to remove the installed files by hand (located in your python site-packages directory).
