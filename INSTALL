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

Additionally, the following python packages should be installed (these will be installed automatically if using *pip*, see `Installation using pip`_):

* NumPy 1.3.0 or higher (http://www.numpy.org/)
* MPI for Python 1.3 or higher (http://mpi4py.scipy.org/) 
* CVXOPT 1.1.3 or higher (http://cvxopt.org/)

Optionally, matplotlib 1.1.0 or higher (http://www.matplotlib.org/) may be manually installed for visualising the algorithm in 2D.

Installation using pip
----------------------
For easy installation, use *pip* (http://www.pip-installer.org/) as root::

    $ [sudo] pip install obb

or alternatively *easy_install*::

    $ [sudo] easy_install obb
    
If you do not have root privileges or you want to install oBB for your private use, you can use::

    $ pip install --user obb
      
which will install oBB in your home directory.

Note that if an older install of oBB is present on your system you can use::

    $ [sudo] pip install --upgrade obb
      
to upgrade oBB to the latest version.

Manual installation
-------------------
Alternatively, you can download the source code and unpack as follows::

    $ wget http://pypi.python.org/packages/source/o/oBB/oBB-X.X.tar.gz
    $ tar -xzvf oBB-X.X.tar.gz
    $ cd oBB-X.X

and then build and install manually using::

    $ python setup.py build
    $ [sudo] python setup.py install

If you do not have root privileges or you want to install oBB for your private use, you can use::

    $ python setup.py install --user
    
instead.    

Testing
-------
oBB includes a command line test script to check that the installation was successfull. To run the test simply type the following into your shell::

    $ test_obb

This will run oBB using MPI on one processor core for a simple 2D sum of sines problem.

Note that if using the MPICH implementation of MPI you first need to start an mpd daemon in the background::

    $ mpd &

but this is not necessary for other MPI implementations, e.g. OpenMPI.

Uninstallation
--------------
If oBB was installed using *pip* you can uninstall as follows::

    $ [sudo] pip uninstall obb

If oBB was installed manually you have to remove the installed files by hand (located in your python site-packages directory).
