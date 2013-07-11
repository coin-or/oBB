# Generate all possible exponents
def generate_exponents(n_variables, order):
    """
    Find the exponents of a multivariate polynomial expression of order
    `order` and `n_variable` number of variables. 
    """
    pattern = [0] * n_variables
    for current_sum in range(0, order+1):
        pattern[0] = current_sum
        yield tuple(pattern)
        while pattern[-1] < current_sum:
            for i in range(2, n_variables + 1):
                if 0 < pattern[n_variables - i]:
                    pattern[n_variables - i] -= 1
                    if 2 < i:
                        pattern[n_variables - i + 1] = 1 + pattern[-1]
                        pattern[-1] = 0
                    else:
                        pattern[-1] += 1
                    break
            yield tuple(pattern)
        pattern[-1] = 0			
	
# Generate random polynomial	 
def randpoly(n,d):
	"""
	A random polynomial in n variables
	of degree at most d is
	returned as a tuple of two lists:
	coefficients and exponents.
	"""
	from random import uniform
	E = [g for g in generate_exponents(n,d)]
	C = [uniform(0,10) for i in range(len(E))]
	return (C,E)
	 
		 
def strmon(S,c,e):
	"""
	Returns a string representation of a
	monomial using the list of symbols in S,
	the coefficient c, and the exponents e.
	"""
	r = (' + ' + str(c))
	for i in range(len(e)):
		if e[i] > 0: r += ('*' + S[i])
		if e[i] > 1: r += ('**' + str(e[i]))
	return r
			 
def strpoly(S,C,E):
	"""
	Returns a string representation of
	a polynomial in the variables in S with
	coefficients in C and exponents in E.
	"""
	p = ""
	for i in range(len(C)):
		p += strmon(S,C[i],E[i])
	return p	 
				 
# Interval arithmetic functions:
# Powers:
# For n zero:
# [a,b]^0 = [a^0,b^0]
# For n odd:
# [a,b]^n = [a^n,b^n]
# For n even:
# [a,b]^n = [a^n,b^n] if a >= 0
# [a,b]^n = [b^n,a^n] if b < 0 
# [a,b]^n = [0,max{a^n,b^n}] o/w
def intpow(a,b,n):
	
	# n zero
	if(n == 0):
		return a**n,b**n
	# n odd
	elif(n%2):
		return a**n,b**n
	# n even
	else:
		if(a >= 0):
			return a**n,b**n
		elif(b < 0):
			return b**n,a**n	
		else:
			return 0,max(a**n,b**n)		

# Multiplication
# [a,b]x[c,d] = [min{ac,ad,bc,bd},max{ac,ad,bc,bd}]
def intmult(a,b,c,d):
	
	return min(a*c,a*d,b*c,b*d),max(a*c,a*d,b*c,b*d)				 
			
def evtupoly(n,d):
	"""
	Generates a random polynomial in n 
	variables with largest degree d.
	"""
	from scipy.misc import comb
	from sympy import var, Poly, Matrix
	from numpy import zeros, empty
	
	m = comb(d-1+n,n,exact=1)
	print 'number of terms: ', m, '+', n, 'leading terms'
	print ''
	(C,E) = randpoly(n,d-1)
	for i in range(0,n):
		C.append(10)
		t = [0 for j in range(0,n)]
		t[i] = d
		E.append(tuple(t))
	print 'coefficients:'
	print C 
	print ''
	print 'exponents:'
	print E
	print ''
	
	# Set up variables
	S = ['x' + str(i) for i in range(1,n+1)]
	
	# Construct symbolic polynomial
	p = Poly(strpoly(S,C,E),var(S))
   	print 'sympy expression:'
	print p
	print ''
	
	# Construct polynomial function
	def pf(x):
		return float(p.eval(zip(S,x)))
	
	# Calculate symbolic gradient
	grad = Matrix(n,1, lambda i,j: p.diff(S[i]))
	print 'sympy gradient:'
	print grad
	print ''
	
	# Construct gradient function
	def pg(x):
		gval = zeros(n)
		for i in range(0,n):
			gval[i] = float(grad[i].eval(zip(S,x)))
		return gval	
	
	# Calculate symbolic Hessian
	H = Matrix(n,n, lambda i,j: (p.diff(S[i])).diff(S[j]))
	print 'sympy Hessian:'
	print H
	print ''
	
	# Construct Hessian function
	def pH(x):
		Hval = zeros((n,n))
		for i in range(0,n):
			for j in range(0,n):
				Hval[i,j] = float(H[i,j].eval(zip(S,x)))
		return Hval
	
	# Calculate symbolic derivative Tensor
	T = empty((n,n,n),dtype='object') 
	for i in range(0,n):
			for j in range(0,n):
				for k in range(0,n):
					T[i,j,k] = ((p.diff(S[i])).diff(S[j])).diff(S[k])
	
	print 'sympy Tensor:'
	print T
	print ''
	
	# Construct Hessian bounding function
	def bndH(l,u):
		"""
		Bound Hessian over l_i <= x_i <= u_i
		"""
		LH = zeros((n,n))
		UH = zeros((n,n))
		
		# Lower bound
		# See which entries have even exponent
		for i in range(0,n):
			for j in range(0,n):
				
				# Get exponents and coefficients
				E = H[i,j].monoms()
				C = H[i,j].coeffs()
				
				# For each term
				for t in range(0,len(E)):
					
					# Containers
					L = 1
					U = 1
					
					# For each variable
					for p in range(0,n):
						
						# Get interval approx.
						Lp,Up = intpow(l[p],u[p],E[t][p])
					
						# Multiply onto previous				
						L,U = intmult(L,U,Lp,Up)
					
					LH[i,j] = LH[i,j] + C[t]*L
					UH[i,j] = UH[i,j] + C[t]*U
		
		return LH, UH		

	# Construct derivative tensor bounding function
	def bndT(l,u):
		"""
		Bound derivative Tensor over l_i <= x_i <= u_i
		"""
		LT = zeros((n,n,n))
		UT = zeros((n,n,n))
		
		# Lower bound
		# See which entries have even exponent
		for i in range(0,n):
			for j in range(0,n):
				for k in range(0,n):
				
					# Get exponents and coefficients
					E = T[i,j,k].monoms()
					C = T[i,j,k].coeffs()
					
					# For each term
					for t in range(0,len(E)):
						
						# Containers
						L = 1
						U = 1
						
						# For each variable
						for p in range(0,n):
							
							# Get interval approx.
							Lp,Up = intpow(l[p],u[p],E[t][p])
						
							# Multiply onto previous				
							L,U = intmult(L,U,Lp,Up)
						
						LT[i,j,k] = LT[i,j,k] + C[t]*L
						UT[i,j,k] = UT[i,j,k] + C[t]*U
		
		return LT, UT

	# Return f(x), g(x), H(x) and H, T bounding functions
	return pf, pg, pH, bndH, bndT
				 