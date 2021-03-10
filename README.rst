=========================================================================
OBB: Overlapping Branch and Bound |License| |Build Status| |PyPI Version|
=========================================================================
oBB is an algorithm for the parallel global optimization of functions with Lipchitz continuous gradient or Hessian.

This is an implementation of the algorithm from our paper:
`Branching and Bounding Improvements for Global Optimization Algorithms with Lipschitz Continuity Properties <http://dx.doi.org/10.1007/s10898-014-0199-6>`_ 
C. Cartis, J. M. Fowkes and N. I. M. Gould. Journal of Global Optimization, vol. 61, no. 3, pp. 429–457, 2015.

The latest version contains an optional range reduction strategy that improves performance in many cases but may not always guarantee global optimality. For details please see the Master's thesis:  
`A Branch and Bound Algorithm for the Global Optimization and its Improvements <http://people.maths.ox.ac.uk/cartis/papers/Thesys_Alberto_Guida.pdf>`_ 
A. Guida. Master's Thesis, Faculty of Engineering, University of Florence, 2015.

Documentation
-------------
HTML documentation is available at http://packages.python.org/oBB

Requirements
------------
oBB requires the following software to be installed:

* Python 2.6/2.7 or Python 3 (http://www.python.org/)
* A working implementation of MPI-2 (e.g. OpenMPI, http://www.open-mpi.org/)

Additionally, the following python packages should be installed (these will be installed automatically if using *pip*, see `Installation using pip`_):

* NumPy 1.3.0 or higher (http://www.numpy.org/)
* MPI for Python 1.3 or higher (http://mpi4py.scipy.org/) 
* CVXOPT 1.1.3 or higher (http://cvxopt.org/)

Optionally, matplotlib 1.1.0 or higher (http://www.matplotlib.org/) may be manually installed for visualising the algorithm in 2D.

Installation using pip
----------------------
For easy installation, use *pip* (http://www.pip-installer.org/) as root::

    $ [sudo] pip install --pre obb

or alternatively *easy_install*::

    $ [sudo] easy_install obb
    
If you do not have root privileges or you want to install oBB for your private use, you can use::

    $ pip install --pre --user obb
      
which will install oBB in your home directory.

Note that if an older install of oBB is present on your system you can use::

    $ [sudo] pip install --pre --upgrade obb
      
to upgrade oBB to the latest version.

Manual installation
-------------------
Alternatively, you can download the source code and unpack as follows::

    $ wget https://pypi.io/packages/source/o/oBB/oBB-0.8b.zip
    $ unzip oBB-0.8b.zip
    $ cd oBB-0.8b

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

Bugs
----
Please report any bugs using GitHub's issue tracker.

License
-------
This algorithm is released under the GNU LGPLv3 license.

.. |License| image::  https://img.shields.io/badge/License-LGPL%20v3-blue.svg
             :target: https://www.gnu.org/licenses/lgpl-3.0
             :alt: GNU LGPL v3 License
.. |Build Status| image::  https://travis-ci.com/coin-or/oBB.svg?branch=master
                  :target: https://travis-ci.com/coin-or/oBB
.. |PyPI Version| image:: https://img.shields.io/pypi/v/oBB.svg
                  :target: https://pypi.python.org/pypi/oBB
