Installing oBB
==============

Requirements
------------
oBB requires the following software to be installed:

* `Python 2.6 to 2.7 <http://www.python.org/>`_
* A working implementation of MPI-2 (e.g. `OpenMPI <http://www.open-mpi.org/>`_ or `MPICH <http://www.mpich.org/>`_)

Additionally, the following python packages should be installed (these will be installed automatically if using `pip <http://www.pip-installer.org/>`_, see `Installation using pip`_):

* `NumPy 1.3.0 or higher <http://www.numpy.org/>`_ 
* `MPI for Python 1.3 or higher <http://mpi4py.scipy.org/>`_
* `CVXOPT 1.1.3 or higher <http://cvxopt.org/>`_ 

Optionally, `matplotlib 1.1.0 or higher <http://www.matplotlib.org/>`_ may be manually installed for visualising the algorithm in 2D.

Installation using pip
----------------------
For easy installation, use `pip <http://www.pip-installer.org/>`_:

 .. code-block:: bash

    $ [sudo] pip install obb

or alternatively *easy_install*:

 .. code-block:: bash

    $ [sudo] easy_install obb
    
If you do not have root privileges or you want to install oBB for your private use, you can use:

 .. code-block:: bash

    $ pip install --user obb
      
which will install oBB in your home directory.

Manual installation
-------------------
Alternatively, you can download the source code and unpack as follows:

 .. code-block:: bash

    $ wget http://pypi.python.org/packages/source/o/oBB/oBB-X.X.tar.gz
    $ tar -xzvf oBB-X.X.tar.gz
    $ cd oBB-X.X

and then build and install manually using:

 .. code-block:: bash

    $ python setup.py build
    $ [sudo] python setup.py install

If you do not have root privileges or you want to install oBB for your private use, you can use:

 .. code-block:: bash

    $ python setup.py install --user
    
instead.    

Uninstallation
--------------
If oBB was installed using `pip <http://www.pip-installer.org/>`_ you can uninstall as follows:

 .. code-block:: bash

    $ [sudo] pip uninstall obb

If oBB was installed manually you have to remove the installed files by hand (located in your python site-packages directory).

