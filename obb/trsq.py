from __future__ import division
def trsq(Lg, delta, xc, f, g):

    # import necessary functions
    from numpy import dot, inf
    from numpy.linalg import norm

    # our own functions
    from newton import newton

    # Define function to minimise
    fun = lambda x: f + dot(x,g) + (Lg / 2) * norm(x) ** 2
    # n.b. x -> x - xc for convenience

    # Case a) Trust-region inactive

    # Quadratic order approx. solution
    xq1 = -g / Lg

    # Check if q. ord. root is within trust region
    if (norm(xq1) < delta):
        bndq1 = fun(xq1)
        xbq1 = xq1
    else: # No solution
        bndq1 = inf
        xbq1 = inf

    # Case b) Trust-region active

    # Initial perturbation
    l = -Lg + 1e-5

    # Define nfq(l) to find quadratic approx. roots
    def nfq(l):

        # Find x(l)
        xl = -g / (l + Lg)

        # Calculate |xl|-delta (for newton stopping rule)
        xlmd = norm(xl) - delta

        # Calculate f(l) for p=-1
        fl = 1/norm(xl) - 1/delta

        # Find x'(l)
        xlp = g / ((l + Lg) ** 2)

        # Calculate f'(l) for p=-1
        flp = -dot(xl,xlp) * (dot(xl,xl) ** (-1.5))

        # Calculate increment
        dl = fl / flp

        # Set Delta
        Delta = delta

        return xlmd, dl, Delta

    # Run newton
    l = newton(nfq, l)

    # Given l, find xq2
    xq2 = -g / (l + Lg)
    bndq2 = fun(xq2)
    xbq2 = xq2

    # Return minimum of bndq1 and bndq2
    if (bndq1 < bndq2):
        bnd = bndq1
        xb = xbq1 + xc# since x -> x - xc
    else:
        bnd = bndq2
        xb = xbq2 + xc# since x -> x - xc

    return bnd, xb
