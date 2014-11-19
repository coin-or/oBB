#!/usr/bin/env python
# Setup script for oBB.
from setuptools import setup
#from numpy.distutils.core import setup, Extension
from sys import version_info, exit

# Make sure the correct version of Python is installed.
if (version_info[:2] < (2,6))or(version_info[0] == 3):
    print "oBB requires Python 2.6 to 2.7. Python %d.%d detected" % \
          version_info[:2]
    exit(-1)

# NAG Library QP Solver Fotran extension paths
#src = ['obb/qpsolver.pyf','obb/qpsolver.f90']
#extr_link_arg = []
#libs = ['$nag_dir/lib/libnag_mkl.so']
#lib_dirs = []
#inc_dirs=['$nag_dir/nag_interface_blocks']
#extr_obj=['$nag_dir/lib/libnag_mkl.so']

# Get package version
exec(open('obb/version.py').read())

# Setup package
setup(
    name='oBB',
    version=__version__ ,
    description='Parallel global optimization of Hessian Lipschitz continuous functions.',
    author='J. Fowkes',
    author_email='jaroslav.fowkes@ed.ac.uk',
    packages=['obb', 'obb.test'],
    scripts=['bin/sins.sh','bin/sins.py'],
    url='http://pypi.python.org/pypi/oBB/',
    license='LGPLv2',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.3.0",
        "mpi4py >= 1.3",
        "cvxopt >= 1.1.3",
        #"sympy >= 0.7.1",
        #"matplotlib >= 1.1.0",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 2.6",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    # NAG library QP solver Fotran extension
    #ext_modules=[Extension('qpsolver',
    #                           sources=src,
    #                           libraries=libs,
    #                           #library_dirs=[],
    #                           include_dirs=inc_dirs,
    #                           extra_objects=extr_obj,
    #                           ),
    #],
    zip_safe=False)
