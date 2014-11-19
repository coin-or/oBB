from __future__ import division
def trsc(sigma, delta, xc, f, g, H, D):

    # import necessary functions
    from numpy import diag, min, argmin, identity, dot, inf, sqrt
    from numpy.linalg import eig, solve, norm

    # our own functions
    from newton import newton
    from newtonw import newtonw

    # Define function to minimise
    fun = lambda x: f + dot(x,g) + 0.5 * dot(x,dot(H,x)) + (sigma / 3) * norm(x) ** 3
    # n.b. x -> x - xc for convenience

    # Return roots of quadratic equation
    def quadeq(a, b, c):

        x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a)
        x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a)

        return x1, x2

    # Case a) Trust-region inactive

    # Given w, get eig. dec.
    LV, UV = eig(H)

    # Get minimum eigenvalue
    l1 = min(LV)
    il = argmin(LV)
    LV = diag(LV) # diagonalise LV

    # Define nf1(w) to find roots
    def nf1(w):

        # Find x(w)
        xw = solve((LV + w * identity(D)),dot(-UV.T,g))

        # Calculate |xw|-w/sigma (for newton stopping rule)
        xwmws = norm(xw) - w / sigma

        # Calculate f(w)for p=-1
        fw = 1 / norm(xw) - sigma / w

        # Find x'(w)
        xwp = solve((LV + w * identity(D)),(-xw))

        # Calculate f'(w) for p=-1
        fwp = -dot(xw,xwp) * (dot(xw,xw) ** (-1.5)) + sigma / (w * w)

        # Calculate increment
        dw = fw / fwp

        # Set w/sigma
        Wsigma = w / sigma

        return xwmws, dw, Wsigma

    # No interior solution for sigma < 0 and H < 0
    if (not ((sigma < 0) and (l1 < 0))):

        # If sigma is negative
        if (sigma < 0):

            # Initial perturbation
            w = -1e-5

            # Lower bound on w
            wl = -l1

        else:

            # Find singularity
            sp = max(0, -l1)

            # Initial perturbation
            w = sp + 1e-5

            # Lower bound on w
            wl = 0

        # Find x(w)
        xw = solve((LV + w * identity(D)),dot(-UV.T,g))

        # Check for Hard Case
        if (norm(xw) <= w / sigma):
            #disp('Hard Case! TR Inactive')

            u = UV[:,il]

            # Find quadratic roots
            [alpha1, alpha2] = quadeq(dot(u,u), 2 * dot(xw,u), dot(xw,xw) - (w / sigma) ** 2)
            xbh1 = xw + alpha1 * u
            bndh1 = fun(xbh1)
            xbh2 = xw + alpha2 * u
            bndh2 = fun(xbh2)

            # Return minimum of bndh1 and bndh2
            if (bndh1 < bndh2):
                bnd1 = bndh1
                xb1 = xbh1
            else:
                bnd1 = bndh2
                xb1 = xbh2

        else:

            # Run newton
            w = newtonw(nf1, w, wl)

            # Given w, find x
            x = solve((H + w * identity(D)),(-g))


            # Check if root is within trust region
            if (norm(x) < delta):
                bnd1 = fun(x)
                xb1 = x
            else:        # No solution
                bnd1 = inf
                xb1 = inf

    else:# No solution
        bnd1 = inf
        xb1 = inf

    # Case b) Trust-region active

    # Initial perturbation
    l = -l1 - sigma * delta + 1e-5

    # Define nf2(l) to find roots
    def nf2(l):

        # Find x(l)
        xl = solve((LV + (l + sigma * delta) * identity(D)),dot(-UV.T,g))


        # Calculate |xl|-delta (for newton stopping rule)
        xlmd = norm(xl) - delta

        # Calculate f(l) for p=-1
        fl = 1 / norm(xl) - 1 / delta

        # Find x'(l)
        xlp = solve((LV + (l + sigma * delta) * identity(D)),(-xl))


        # Calculate f'(l) for p=-1
        flp = -dot(xl,xlp) * (dot(xl,xl) ** (-1.5))

        # Calculate increment
        dl = fl / flp

        # Set Delta
        Delta = delta

        return xlmd, dl, Delta

    # Find x(l)
    xl = solve((LV + (l + sigma * delta) * identity(D)),dot(-UV.T,g))


    # Check for hard case
    if (norm(xl) <= delta):
        #disp('Hard Case! TR Active')

        u = UV[:,il]

        # Find quadratic roots
        [alpha1, alpha2] = quadeq(dot(u, u), 2 * dot(xl,u), dot(xl,xl) - (delta * delta))
        xbh1 = xl + alpha1 * u
        bndh1 = fun(xbh1)
        xbh2 = xl + alpha2 * u
        bndh2 = fun(xbh2)

        # Return minimum of bndh1 and bndh2
        if (bndh1 < bndh2):
            bnd2 = bndh1
            xb2 = xbh1
        else:
            bnd2 = bndh2
            xb2 = xbh2

    else:

        # Run newton
        l = newton(nf2, l)

        # Given l, find x
        x = solve((H + (l + sigma * delta) * identity(D)),(-g))

        bnd2 = fun(x)
        xb2 = x

    # Return minimum of bnd1 and bnd2
    if (bnd1 < bnd2):
        bnd = bnd1
        xb = xb1 + xc# since x -> x - xc
    else:
        bnd = bnd2
        xb = xb2 + xc# since x -> x - xc

    return bnd, xb
