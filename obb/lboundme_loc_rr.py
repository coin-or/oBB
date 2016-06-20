from __future__ import division
def bound(c,(bndH,method),d1,f,g,d2,D,A=None,b=None,E=None,d=None,U=None,optimality=None,initialRadius=None):

    # Norm function
    from numpy import ones, sqrt
    from numpy.linalg import norm

    # Our own functions
    from trsq import trsq
    from mest import mest
    from domain_reduction import RR

    lr, ur = RR(c,bndH,f,g,D,A,b,E,d,U,optimality)

    if optimality is not None:

        # bounding box containing orginal circle
        l1 = c.xc - c.r * ones(D)
        u1 = c.xc + c.r * ones(D)

        # Radius and centre ball containing the reduced box
        r = norm(ur-lr)/2 # Radius
        xc = (ur+lr)/2 # Centre

        # bounds box containing ball around domain reduced
        l2 = xc - r * ones(D)
        u2 = xc + r * ones(D)

        check = 1
        for i in range(0,D):
            if l2[i]<l1[i] or u2[i]>u1[i]:
                check = 0
                break

        #if the box reduced is not completely contained in the orginal box
        if check == 0:

            # Revert to the orginal ball depending o the level of depth (= initialRadius/const)
            #if c.r >= initialRadius/8:
            #    l=l1
            #    u=u1
            #    xc = c.xc
            #    r = c.r

            #else:

            r = r/sqrt(D)
            l = xc - r * ones(D)
            u = xc + r * ones(D)

        else:
            l = l2
            u = u2
    else:
        l = lr
        u = ur
        xc = c.xc
        r = c.r

    # Get bounds on RBF Hessian over [l,u]
    LH, UH = bndH(l, u)

    # Estimate smallest eigenvalue on [l,u]
    k = mest(LH,UH,method)

    # Evaluate function and gradient at xc
    fxc = f(xc)
    gxc = g(xc)

    # If gradient nonzero use trsq
    if(norm(gxc) > 1e-10):

        lbound,_ = trsq(k, r, xc, fxc, gxc)

    else: # Calculate lower bound when gradient zero

        if(k < 0):

            lbound = fxc + k*(r**2)

        else:
            lbound = fxc

    return lbound,lr,ur
