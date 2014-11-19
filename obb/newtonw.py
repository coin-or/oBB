# -*- coding: utf-8 -*-
# Newton's Method for finding the root of a scalar function
def newtonw(nf, x, xl):

    # Inputs
    # nf - calculates f(x)/df(x) for given x
    # x - initial guess
    # xl - lower bound on x

    # Outputs
    # x - solution (hopefully)

    # Parameters
    #Ke = 0.29 # Tolerance
    Ke = 0.01
    fx, dx, Delta = nf(x)# Increment

    while((abs(fx) > Ke * Delta) and (x > xl)):

        # Increment x
        x = x - dx

        # Calculate new increment
        fx, dx, Delta = nf(x)

    return x
