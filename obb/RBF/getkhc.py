from __future__ import division
# Find lower and upper bounds on cubic RBF third order derivative tensor over box
def getkhc(l, u):

    # necessary functions
    from numpy import diag, dot, zeros, array, sqrt, sum, isnan, vstack, amax, amin, maximum
    from numpy.linalg import norm

    # load data
    from config import N, D, a, x, tl, tl2

    # Simple function
    spls = lambda s,j: tl2[j]*s[j]/norm(dot(diag(tl),s))

    # Simple function for getlhc
    def spl(s): 
	
	nrm_st = norm(dot(diag(tl),s))
	y = zeros(D)
	for i in range(0,D):
		y[i] = tl2[i]*s[i]/nrm_st
	y[isnan(y)] = 0
	return y

    # Recursively find maximum of spl on box
    def recurse(cd, d, v, i):

	# Set up necessary variables
	vl = zeros(cd)
	vu = zeros(cd)

	# Retain previous values
	vl[0:cd - 1] = v[0:cd - 1]
	vu[0:cd - 1] = v[0:cd - 1]
	vl[cd-1] = l[cd-1]
	vu[cd-1] = u[cd-1]

	# If recursion has come to an end 
	if (cd == d):
	    ml = spl(vl - x[i,:])
	    mu = spl(vu - x[i,:])
	# Otherwise recurse     
	else:
	    mll, mul = recurse(cd + 1, d, vl, i)
	    mlu, muu = recurse(cd + 1, d, vu, i)
	    ml = amin(vstack([mll,mul,mlu,muu]),axis=0)
	    mu = amax(vstack([mll,mul,mlu,muu]),axis=0)
	return ml, mu

    # Get bounds on norm over box

    # Storage
    cml = zeros((N,D))
    cmu = zeros((N,D))

    # Bound for each centre
    for i in range(0,N):

	# Coupled simple bounds at vertices
	cvl, cvu = recurse(1, D, array([1]), i)

        # Include centres
	ccl = zeros(D)
	ccu = zeros(D)

	# Optimal norm depends on direction
	for j in range(0,D):

		# Modified xm    
		xmj = x[i,:].copy()   

		# Lower bound
		xmj[j] = l[j] 
		ccl[j] = spls(xmj-x[i,:],j)
		
		# Upper bound
		xmj[j] = u[j] 
		ccu[j] = spls(xmj-x[i,:],j)

	# Replace NaNs with zero
	ccl[isnan(ccl)] = 0	    
	ccu[isnan(ccu)] = 0

	# Componentwise bounds
        cml[i,:] = amin(vstack([cvl,ccl,cvu,ccu]),axis=0)
        cmu[i,:] = amax(vstack([cvl,ccl,cvu,ccu]),axis=0)


    # Calculate lower and upper bounds on third order RBF tensor
    LT = zeros((D, D, D))
    UT = zeros((D, D, D))
    for i in range(0,D):
        for j in range(0,D):
            for k in range(0,D):
                for p in range(0,N):
                    # Positive coefficient
                    if (a[p] > 0):
                        UT[i, j, k] = UT[i, j, k] - a[p] * 3 * cml[p, i] * cml[p, j] * cml[p, k]
                        if ((i == j) and (j == k)):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 9 * tl2[i] * cmu[p, i] 
                        elif (i == j):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[i] * cmu[p, k]
                        elif (j == k):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[j] * cmu[p, i] 
                        elif (i == k):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[k] * cmu[p, j] 
			    
			LT[i, j, k] = LT[i, j, k] - a[p] * 3 * cmu[p, i] * cmu[p, j] * cmu[p, k] 
                        if ((i == j) and (j == k)):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 9 * tl2[i] * cml[p, i]
                        elif (i == j):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[i] * cml[p, k]
                        elif (j == k):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[j] * cml[p, i]
                        elif (i == k):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[k] * cml[p, j]    
                        
                    # Negative coefficient
                    else:
                        UT[i, j, k] = UT[i, j, k] - a[p] * 3 * cmu[p, i] * cmu[p, j] * cmu[p, k] 
                        if ((i == j) and (j == k)):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 9 * tl2[i] * cml[p, i]
                        elif (i == j):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[i] * cml[p, k]
                        elif (j == k):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[j] * cml[p, i]
                        elif (i == k):
                            UT[i, j, k] = UT[i, j, k] + a[p] * 3 * tl2[k] * cml[p, j]

			LT[i, j, k] = LT[i, j, k] - a[p] * 3 * cml[p, i] * cml[p, j] * cml[p, k]
                        if ((i == j) and (j == k)):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 9 * tl2[i] * cmu[p, i] 
                        elif (i == j):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[i] * cmu[p, k]
                        elif (j == k):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[j] * cmu[p, i] 
                        elif (i == k):
                            LT[i, j, k] = LT[i, j, k] + a[p] * 3 * tl2[k] * cmu[p, j]

    # Return bounds
    return LT,UT
    