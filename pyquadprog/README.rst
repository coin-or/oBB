PyQuadProg
=============
A Python interface to QuadProg++ (Luca Di Gaspero, http://www.diegm.uniud.it/digaspero/)
which solves dense convex QPs using Goldfarb and Idnani's method (1983).

Right now PyQuadProg only contains the basic functionalities.
Any feedback is highly appriciated (mehdi.towhidi@gerad.ca). Thank you!

Installation
============

STEP 1:
    Install QuadProg++. Run::

        $ ./configure
        $ make
        $ make install

STEP 2:
    Create an environment variable called QUADPROG_DIR pointing to your
    installation of Coin. For example::

        $ export QUADPROG_DIR=/Users/mehdi/quadprog

    You may also add this line to your ~/.bash_rc or ~/.profile to make
    it persistent.

STEP 3:
    Install PyQuadProg. Go to PyQuadProg's root directory and run::

        $ python setup.py install

While your compiling, if you get complains about some members being private,
please go to QuadProg++ source directory and open array.hh. Change the
private members of *Vector* and *Matrix* to public attributes
(e.g. by deleting the lines containing "private:").

Usage
=======

Please see test.py for a simple QP solution example.

Here is an almost exact copy of the comments in QuadProg++.hh (by Luca Di Gaspero, http://www.diegm.uniud.it/digaspero/)

The quadprog_solve() function implements the algorithm of Goldfarb and Idnani
for the solution of a (convex) Quadratic Programming problem
by means of an active-set dual method.

The problem is in the form:


.. :math:`min \frac{1}{2} xGx + g_0x`

min 0.5 * x G x + g0 x

subject to:

.. :math:`C_E^T x + c_e_0 = 0`

.. :math:`C_I^T x + c_i_0 \geq 0`

CE^T x + ce0 = 0

CI^T x + ci0 >= 0

The matrix and vectors dimensions are as follows:

.. :math:`G : n \times n`

.. :math:`g_0 : n`


.. :math:`C_E : n \times p`

.. :math:`c_e_0 : p`


.. :math:`C_I : n \times m`

.. :math:`c_i_0 : m`


.. :math:`x : n`


G: n * n

g0: n

CE: n * p

ce0: p

CI: n * m

ci0: m

x: n


The function will return the cost of the solution written in the x vector or
infinity if the problem is infeasible. In the latter case
the value of the x vector is not correct.


Reference
=============

D. Goldfarb, A. Idnani. A numerically stable dual method for solving
strictly convex quadratic programs. Mathematical Programming 27 (1983) pp. 1-33.
