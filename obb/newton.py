# -*- coding: utf-8 -*-
# Newton's Method for finding the root of a scalar function
def newton(nf, x):

    # Inputs
    # nf - calculates f(x)/df(x) for given x
    # x - initial guess

    # Outputs
    # x - solution (hopefully)

    # Parameters
    #itr = 1 # Iteration count
    #mitr = 1000 # Max. iterations
    #Ke = 0.29 # Tolerance
    Ke = 0.01
    #eps = 1e-15 # Tolerance
    fx, dx, Delta = nf(x)# Increment

    #while((abs(dx) > eps)&&(iter < miter))
    while ((abs(fx) > Ke * Delta)):
        #and(x < 1e10))

        # Update iteration count
        #itr = itr + 1

        # Increment x
        x = x - dx

        # Calculate new increment
        fx, dx, Delta = nf(x)

    return x

    # % Check if newton runs for ages
    # if(itr >= mitr):
    #     print "Fix me!"
    #     print x
    #     print Delta
    #     pause;
