from __future__ import division
def bound(c,(bndH,d1),(bndT,method),f,g,H,D,A=None,b=None,E=None,d=None,U=None,optimality=None,initialRadius=None):

    # Our own functions
    from trsc import trsc
    from gcest import gcest
    from domain_reduction import RR
    from numpy.linalg import norm
    from numpy import ones, sqrt

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


    # Get bounds on RBF third order derivative tensor over [l,u]
    LT, UT = bndT(l, u)

    # Estimate smallest eigenvalue on [l,u]
    k = gcest(LT,UT,method)

    lbound,_ = trsc(k / 2, r, xc, f(xc), g(xc), H(xc), D)

    return lbound,lr,ur
