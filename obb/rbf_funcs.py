from __future__ import division
# RBF function, gradient and Hessian
def frdl(s):

    # necessary functions
    from numpy import diag, dot, hstack, array
    from numpy.linalg import norm

    # load data
    from config import N, D, a, x, tl

    phi = lambda x: norm(dot(diag(tl),x)) ** 3

    z = dot(hstack((array([1]), s[0:D])),a[N:N+D+1])

    for i in range(0,N):
        z = z + a[i] * phi(s - x[i,:])

    return z

def grdl(s):

    # necessary functions
    from numpy import diag, dot
    from numpy.linalg import norm

    # load data
    from config import N, D, a, x, tl

    gphi = lambda s: 3 * dot(dot(diag(tl), diag(tl)),s) * norm(dot(diag(tl),s))

    z = a[N+1:N+D+1]

    for i in range(0,N):
        z = z + a[i] * gphi(s - x[i,:])

    return z

def hrdl(s):

    # necessary functions
    from numpy import diag, dot, zeros
    from numpy.linalg import norm

    # load data
    from config import N, D, a, x, tl, tl2

    def hphi(s):

        nrm_st = norm(dot(diag(tl),s))
        h = zeros((D,D))
        for i in range(0,D):
            for j in range(0,D):
                h[i,j] = 3 * tl2[i] * tl2[j] * s[i] * s[j] / nrm_st
                if (i == j):
                    h[i,j] = h[i,j] + 3 * tl2[i] * nrm_st

        return h

    z = zeros((D,D))
    for p in range(0,N):
        z = z + a[p] * hphi(s - x[p,:])

    return z
