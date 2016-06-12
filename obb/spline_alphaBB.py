from __future__ import division
def convexUnderEstimator(f,g,bndH,D,l,u):

    from numpy import zeros, fabs
    from numpy.linalg import lstsq
    from itertools import repeat

    ndiv = D
    inc = (u-l)/ndiv

    lp = zeros((ndiv+1,D))
    lp[0] = l
    for z in range(1,ndiv+1):
        lp[z] = lp[z-1]+inc

    alpha = zeros((ndiv,D))
    beta = zeros((ndiv,D))
    gamma = zeros((ndiv,D))

    # for each sub-box computing alpha
    for z in range(0,ndiv):

        # Get bounds on RBF Hessian over [l,u]
        LH, UH = bndH(lp[z],lp[z+1])

        for i in range(0,D):

            #di = lp[z+1][i]-lp[z][i]
            term=0
            for j in range(0,i):
                #dj = lp[z+1][j]-lp[z][j]
                term += max(fabs(LH[i,j]),fabs(UH[i,j]))#*(dj/di)
            for j in range(i+1,D):
                #dj = lp[z+1][j]-lp[z][j]
                term += max(fabs(LH[i,j]),fabs(UH[i,j]))#*(dj/di)

            alpha[z][i]= 0.5*max(0,-LH[i,i]+term)

    # Building the linear system
    for i in range(0,D):

        # Building the part of the linear system relative to the interval [li,ui]
        a = []
        b = []


        elements = []
        elements.append(lp[0][i])#beta1
        elements.extend(repeat(0,ndiv-1))#beta2...betaN
        elements.append(1)#gamma1
        elements.extend(repeat(0,ndiv-1))#gamma2...gammaN

        a.append(elements)
        b.append(0)

        elements = []
        elements.extend(repeat(0,ndiv-1))#beta1...betaN-1
        elements.append(lp[ndiv][i])#betaN
        elements.extend(repeat(0,ndiv-1))#gamma1...gammaN-1
        elements.append(0)#gammaN

        a.append(elements)
        b.append(0)

        for z in range(1,ndiv):

            elements = []
            elements.extend(repeat(0,z-1))#beta1...beta_i-1
            elements.append(lp[z][i])#betai
            elements.append(-lp[z][i])#beta_i+1
            elements.extend(repeat(0,ndiv-(z+1)))#beta_i+2...betaN
            elements.extend(repeat(0,z-1))#gamma1...gamma_i-1
            elements.append(1)#gammai
            elements.append(-1)#gamma_i+1
            elements.extend(repeat(0,ndiv-(z+1)))#gamma_i+2...gammaN

            a.append(elements)
            b.append(0)

            elements = []
            elements.extend(repeat(0,z-1))#beta1...beta_i-1
            elements.append(1)#betai
            elements.append(-1)#beta_i+1
            elements.extend(repeat(0,ndiv-(z+1)))#beta_i+2...betaN
            elements.extend(repeat(0,ndiv))#gamma1...gammaN


            a.append(elements)
            b.append( (lp[z][i]-lp[z-1][i])*alpha[z-1][i]+(lp[z+1][i]-lp[z][i])*alpha[z][i])

        sol= lstsq(a,b)
        #print 'sol',sol

        for z in range(0,ndiv):
            beta[z][i]= sol[0][z]
            gamma[z][i]= sol[0][z+ndiv]

    def q(x):
        value = 0
        for j in range(0,D):
            for z in range(0,ndiv):
                if x[j]<=lp[z][j] and x[j]<=lp[z+1][j]:
                    value += alpha[z][j]*(x[j]-lp[z][j])*(lp[z+1][j]-x[j]) + beta[z][j]*x[j] + gamma[z][j]

        return value

    def gq(x):
        gq = zeros(D)
        for j in range(1,D):
            for z in range(0,ndiv):
                if x[j]<=lp[z][j] and x[j]<=lp[z+1][j]:
                    gq[j] = alpha[z][j]*(lp[z][j]+lp[z+1][j]-2*x[j]) + beta[z][j]
        return gq

    convexUnder = lambda x:f(x)-q(x)
    gConvexUnder = lambda x: g(x)-gq(x)

    return convexUnder, gConvexUnder
