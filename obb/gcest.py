from __future__ import division
# Use the third order derivative tensor Gershgorin estimation method
def gcest(LT,UT,method):

    # Get dimension
    D = LT.shape[0]

    #    # Positivity/ Negativity tensor checks
    #    from numpy import all
    #    print('LT non-positive: %s') % all(LT <= 0)
    #    print('UT non-negative: %s') % all(UT >= 0)
    #
    #    # Ediag check      ef
    #    print('MA non-positive: %s') % all((UT+LT)/2 <= 0)
    #    print('RA non-negative: %s') % all((UT-LT)/2 >= 0)
    #               Gershgorin
    #    mRA = (UT-LT)/2. # radius tensor.
    #    mMA = (UT+LT)/2. # midpoint tensor
    #    for i in range(0,D):
    #       for j in range(0,D):
    #               for k in range(0,D):
    #                       if((i==j)and(j==k)):
    #                               mRA[i,j,k] = 0
    #                               mMA[i,j,k] = LT[i,j,k]
    #    print('mMA non-positive: %s') % all(mMA <= 0)
    #    print('mRA non-negative: %s') % all(mRA >= 0)
    #
    #    # lbH check (equivalent to Gersh like quad?)
    #    NRA = (LT-UT)/2.
    #    rs = (NRA.sum(axis=1)).sum(axis=1)
    #    A = (LT+UT)/2.
    #    for i in range(0,D):
    #       for j in range(0,D):
    #               for k in range(0,D):
    #                       if((i==j)and(j==k)):
    #                               A[i,j,k] = LT[i,j,k] + (rs[i] - NRA[i,j,k])
    #    print('lbH non-positive: %s') % all(A <= 0)

    # Select estimation method (gc, c)
    # Gershgorin for Tensors
    if(method == 'gc'):

        # Imports
        from numpy import maximum, zeros

        # Calculate max absolute value of bounds
        VT = maximum(abs(LT),abs(UT))

        # Get row plane sums
        rs = (VT.sum(axis=1)).sum(axis=1)

        # Tensor diagonal function
        def diagt(T):
            v = zeros(D)
            for i in range(0,D):
                v[i] = T[i,i,i]
            return v

        # Calculate lower bounds on Gershgorin disks
        G = diagt(LT) - (rs-diagt(VT))

        # Calculate Gershgorin lower bound
        k = min(G)

        # If k negative ok, if k positive need other bound
        if(k < 0):
            pass
            #print('k ok, negative')
        else:
            #print('k positive, using other bound.')
            k = (D**(-0.5))*k
        return k

    # Lh = norm_F(VT) so return -Lh
    elif(method == 'c'):

        # Imports
        from numpy import maximum, sqrt, sum

        # Calculate max absolute value of bounds
        VT = maximum(abs(LT),abs(UT))

        # Calculate frobenius norm of VT
        return -sqrt(sum(sum(sum(VT ** 2))))

    else:

        raise RuntimeError('Method must be one of gc, c.')
