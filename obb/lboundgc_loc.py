from __future__ import division
def bound(c,d1,(bndT,method),f,g,H,D):

    # import necessary functions
    from numpy import ones

    # Our own functions
    from trsc import trsc
    from gcest import gcest

    # Get local Lipschitz constant
    l = c.xc - c.r * ones(D)
    u = c.xc + c.r * ones(D)

    # Get bounds on RBF third order derivative tensor over [l,u]
    LT, UT = bndT(l, u)

    # Estimate smallest eigenvalue on [l,u]
    k = gcest(LT,UT,method)

    # Evaluate function, gradient and Hessian at xc
    fxc = f(c.xc)
    gxc = g(c.xc)
    Hxc = H(c.xc)

    # Lower bound
    lbound,_ = trsc(k / 2, c.r, c.xc, fxc, gxc, Hxc, D)

    return lbound
