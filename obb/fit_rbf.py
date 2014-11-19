# Fit Radial Basis Functions Surrogate
def fit_rbf(fun, x):

    # Load relevant functions
    from numpy import size, zeros, ones, diag, dot, hstack, vstack, array
    from numpy.linalg import norm, solve

    # Get sizes
    N = size(x,0) # Number of samples
    D = size(x,1) # Problem Dimension

    # Generate initial sample
    z = zeros(N)
    for i in range(0,N):
        z[i] = fun(x[i,:])

    # Fit initial linear + rdb surrogate

    # Calculate Polynomials
    P = zeros((N, D + 1))
    for i in range(0,N):
        P[i,:] = hstack(( array([1]), x[i,:] ))

    # Scaling parameters
    tl = ones(D)

    # Basis function
    phi = lambda x: norm(dot(diag(tl),x)) ** 3

    # Calculate rdb's
    Q = zeros((N, N))
    for i in range(0,N):
        for j in range(0,N):
            Q[i,j] = phi(x[i,:]-x[j,:])

    # Add additional constraints
    b = hstack(( z, zeros(D+1) ))

    # Assemble A
    A = vstack(( hstack(( Q, P )), hstack(( P.T, zeros((D+1,D+1)) )) ))
    a = solve(A,b)

    # Write relevant global variables
    import config
    config.N = N
    config.D = D
    config.a = a
    config.x = x
    config.tl = tl
    config.tl2 = config.tl**2
