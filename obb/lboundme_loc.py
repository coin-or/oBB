from __future__ import division
def bound(c,(bndH,method),d1,f,g,d2,D):

    # Norm function
    from numpy import ones
    from numpy.linalg import norm

    # Our own functions
    from trsq import trsq
    from mest import mest

    # Calculate bounding box
    l = c.xc - c.r * ones(D)
    u = c.xc + c.r * ones(D)

    # Get bounds on RBF Hessian over [l,u]
    LH, UH = bndH(l, u)

    # Estimate smallest eigenvalue on [l,u]
    k = mest(LH,UH,method)

    # Evaluate function and gradient at xc
    fxc = f(c.xc)
    gxc = g(c.xc)

    # If gradient nonzero use trsq
    if(norm(gxc) > 1e-10):

        lbound,_ = trsq(k, c.r, c.xc, fxc, gxc)

    else: # Calculate lower bound when gradient zero

        if(k < 0):

            lbound = fxc + k*(c.r**2)

        else:
            lbound = fxc

    return lbound
