from __future__ import print_function
#!/usr/bin/env python
# Setup script for oBB.
from setuptools import setup, Extension
from sys import version_info, exit
from numpy import get_include

# Make sure the correct version of Python is installed.
if (version_info[:2] < (2,6)):
    print("oBB requires Python 2.6/2.7 or Python 3. Python %d.%d detected" % version_info[:2])
    exit(-1)

# Get package version
exec(open('obb/version.py').read())

# Setup QuadProg++ and SLSQP Extension
ext_modules = [Extension('PyQuadProg',
                          sources=['pyquadprog/PyQuadProg.cpp','quadprog/QuadProg++.cc','quadprog/Array.cc'],
                          include_dirs=['pyquadprog','quadprog',get_include()],
                          language='c++'
                          ), 
               Extension('_nlopt',
                          sources=['slsqp/slsqp.c','nlopt/general.c','nlopt/optimize.c','nlopt/options.c','nlopt/stop.c','nlopt/timer.c','nlopt/nlopt_wrap.cxx'],
                          include_dirs=['slsqp','nlopt',get_include()],
                          language='c++'
                          )]

# Setup package
setup(
    name='oBB',
    version=__version__ ,
    description='Parallel global optimization of Hessian Lipschitz continuous functions.',
    author='J. Fowkes',
    author_email='jaroslav.fowkes@ed.ac.uk',
    packages=['obb'],
    scripts=['bin/sines.sh','bin/sines.py','bin/sines_rbf.sh','bin/sines_rbf.py','bin/coconut.sh','bin/coconut.py','bin/test_obb'],
    include_package_data=True,
    package_data={'obb': ['obb/coconut/*','obb/coconut_tol']},
    ext_modules=ext_modules,
    py_modules=['nlopt'],
    url='http://pypi.python.org/pypi/oBB/',
    license='LGPLv3',
    long_description=open('README.rst').read(),
    use_2to3 = True,
    install_requires=[
        "numpy >= 1.3.0",
        "mpi4py >= 1.3",
        "cvxopt >= 1.1.3",
        #"matplotlib >= 1.1.0",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    zip_safe=False)
