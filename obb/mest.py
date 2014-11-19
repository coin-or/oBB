from __future__ import division
# Use one of the Hessian eigenvalue estimation methods from Floudas
def mest(LM,UM,method):

    # Get dimension
    D = LM.shape[0]

    # Select estimation method (g, Hz, lbH, E0, Ediag, q)
    # Gershgorin
    if(method == 'g'):

        # Imports
        from numpy import maximum, diag

        # Calculate max absolute value of bounds
        VM = maximum(abs(LM),abs(UM))

        # Get row sums
        rs = VM.sum(axis=1)

        # Calculate lower bounds on Gershgorin disks
        G = diag(LM) - (rs-diag(VM))

        # Calculate Gershgorin lower bound
        k = min(G)

        return k

    # Hertz
    elif(method == 'Hz'):

        # Imports
        from numpy import array, zeros, inf
        from numpy.linalg import eig
        from itertools import product

        # Calculate orthants
        X = array([p for p in product([-1,1], repeat=int(D))])

        # Calculate min eig. of vertex matrices
        lmin = inf
        for k in range(0,2**(D-1)):

            A = zeros((D,D))
            for i in range(0,D):
                for j in range(0,D):

                    if(i == j):
                        A[i,j] = LM[i,j]
                    else:
                        if(X[k,i]*X[k,j] >= 0):
                            A[i,j] = LM[i,j]
                        else:
                            A[i,j] = UM[i,j]

            LV,_ = eig(A)
            l1 = min(LV)
            if(l1 < lmin):
                lmin = l1

        return lmin

    # Lower bounding Hessian
    elif(method == 'lbH'):

        # Imports
        from numpy import diag
        from numpy.linalg import eig

        # Calculate negative radius matrix
        NRA = (LM-UM)/2.

        # Get row sums of negative radius matrix
        rs = NRA.sum(axis=1)

        # Let A be midpoint matrix on offdiagonals
        A = (LM+UM)/2.

        # Set diagonals
        A.flat[::D+1] = diag(LM) + (rs - diag(NRA))

        # Calculate min eig.
        LV,_ = eig(A)
        l1 = min(LV)

        return l1

    # E-matrix, E = 0
    elif(method == 'E0'):

        # Imports
        from numpy import diag
        from numpy.linalg import eig

        # Calculate modified radius and midpoint matrices
        mRA = (UM-LM)/2. # radius matrix.
        mRA.flat[::D+1] = 0 # set diagonal to 0
        mMA = (UM+LM)/2. # midpoint matrix
        mMA.flat[::D+1] = diag(LM) # set diagonal to diag(LM)

        # Calculate min eig. mMA
        LV,_ = eig(mMA)
        l1 = min(LV)

        # Calculate spectral radius of mRA
        LV,_ = eig(mRA)
        p = max(abs(LV))

        return l1 - p

    # E-matrix, E = diag(RA) (at least as good as Mori and Kokame)
    elif(method == 'Ediag'):

        # Imports
        from numpy.linalg import eig

        # Calculate radius and midpoint matrices
        RA = (UM-LM)/2. # radius matrix.
        MA = (UM+LM)/2. # midpoint matrix

        # Calculate min eig. MA
        LV,_ = eig(MA)
        l1 = min(LV)

        # Calculate spectral radius of RA
        LV,_ = eig(RA)
        p = max(abs(LV))

        return l1 - p

    # Lg = norm_F(VM) so return -Lg
    elif(method == 'q'):

        # Imports
        from numpy import maximum, sqrt, sum

        # Calculate max absolute value of bounds
        VM = maximum(abs(LM),abs(UM))

        # Calculate frobenius norm of VM
        return -sqrt(sum(sum(VM ** 2)))

    else:

        raise RuntimeError('Method must be one of g, Hz, lbH, E0, Ediag, q.')
