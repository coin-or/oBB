from __future__ import division
def RR(c,bndH,f,g,D,A=None,b=None,E=None,d=None,U=None,optimality=None):

    from spline_alphaBB import convexUnderEstimator
    # import necessary functions
    from numpy import dot, array, fabs, zeros, ones, sqrt, delete, Infinity, insert
    from numpy.linalg import norm as linnorm

    import nlopt

    if c.lReduced is not None and c.uReduced is not None:
        l = c.lReduced
        u = c.uReduced
    else:
        # Calculate bounding box
        l = c.xc - c.r * ones(D)
        u = c.xc + c.r * ones(D)

    if optimality is None:
        RR = False
    else:
        RR = True
    if RR is True:
        #tollerance for stopping domain reduction strategy (if is optimality)
        stopTollerance = 0.1

        rel_tol = 1e-6
        if optimality is True:
            max_time = 5
        else:
            max_time = 1

        xc = c.xc
        r = c.r

    verbose = 0
    if verbose == 1:
        print ''
        print '#########start RR#########'
        print 'center initial ball ',xc
        print 'radius initial ball ',r
        print 'l initial ball ',c.xc - c.r * ones(D)
        print 'u initial ball ',c.xc + c.r * ones(D)
        print ''
        print 'l initial ball after intersection',l
        print 'u initial ball after intersection',u
        print ''

    if RR is True:

        while True:

            # Building the convex underestimator
            if optimality is True:
                convexUnder, gconvexUnder = convexUnderEstimator(f,g,bndH,D,l,u)

            # Start solving reduction's problem with NLopt solver

            lcopy = array(l)
            ucopy = array(u)

            for i in range(0,D):

                def myfunc(x,grad):

                    if grad.size >0:
                        c = zeros(D)
                        c[i] = 1
                        grad[:] = c

                    return x[i]

                def linconstraint(x,grad,M,v):
                    if grad.size>0:
                        grad[:] = M

                    return dot(M,x)-v

                def radcons(x,grad):
                    if grad.size>0:
                        for i in range(0,D):
                            grad[i] = 2*(x[i]-xc[i])

                    return linnorm(x-xc)**2 - r**2

                def convexUndercons(x,grad):
                    if grad.size>0:
                        grad[:]=gconvexUnder(x)

                    return convexUnder(x)-U

                opt = nlopt.opt(nlopt.LD_SLSQP,D)
                opt.set_lower_bounds(l)
                opt.set_upper_bounds(u)
                opt.set_min_objective(myfunc)

                opt2 = nlopt.opt(nlopt.LD_SLSQP,D)
                opt2.set_upper_bounds(u)
                opt2.set_max_objective(myfunc)

                if E.any():
                    for j in range(0,len(E)):
                        opt.add_equality_constraint(lambda x, grad:linconstraint(x,grad,E[j],d[j]),1e-8)
                        opt2.add_equality_constraint(lambda x, grad:linconstraint(x,grad,E[j],d[j]),1e-8)
                if A.any():
                    for j in range(0,len(A)):
                        opt.add_inequality_constraint(lambda x,grad:linconstraint(x,grad,A[j],b[j]), 1e-8)
                        opt2.add_inequality_constraint(lambda x,grad:linconstraint(x,grad,A[j],b[j]), 1e-8)

                opt.add_inequality_constraint(lambda x,grad: radcons(x,grad), 1e-8)
                opt2.add_inequality_constraint(lambda x,grad: radcons(x,grad), 1e-8)

                if optimality is True:
                    opt.add_inequality_constraint(lambda x,grad: convexUndercons(x,grad),1e-8)
                    opt2.add_inequality_constraint(lambda x,grad: convexUndercons(x,grad),1e-8)

                opt.set_xtol_rel(rel_tol)
                opt.set_maxtime(max_time)
                try:
                    sol = opt.optimize(l)

                    li = opt.last_optimum_value()
                    if li>=l[i] and li<=u[i] and opt.last_optimize_result()>0:
                        l[i]=li

                    if verbose==1:
                        print 'component l',i
                        print 'sol',sol
                        print 'value of the optimum',li
                        if optimality is True:
                            print 'value convexUnder',convexUnder(sol)
                            print 'value f', f(sol)
                        print 'distance',linnorm(xc-sol)
                        print 'result code',opt.last_optimize_result()
                        print ''
                except:
                    pass

                opt2.set_lower_bounds(l)
                opt2.set_xtol_rel(rel_tol)
                opt2.set_maxtime(max_time)
                try:
                    sol = opt2.optimize(u)

                    ui = opt2.last_optimum_value()
                    if ui <=u[i] and ui>=l[i] and opt2.last_optimize_result()>0:
                        u[i]=ui

                    if verbose==1:
                        print 'component u',i
                        print 'value of the optimum',ui
                        print 'sol',sol
                        if optimality is True:
                            print 'value convexUnder',convexUnder(sol)
                            print 'value f', f(sol)
                        print 'distance',linnorm(xc-sol)
                        print 'result code',opt2.last_optimize_result()
                        print ''
                except:
                    pass


            if optimality is False:
                break

            if fabs(linnorm(l-lcopy))<=stopTollerance and fabs(linnorm(u-ucopy))<=stopTollerance :
                break

            if verbose == 1:
                print 'intermediate l reduced ',l
                print 'intermediate u reduced ',u
                print ''

    if verbose == 1:
        print 'final l reduced ',l
        print 'final u reduced ',u
        print '#########endRR#########'

    return l,u
