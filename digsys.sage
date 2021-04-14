import numpy as np

testBound = 1000

class RotationDS:
	
	'''Constructs a digit system for a rotation matrix'''
	
	def __init__(self, base_matrix, rem='normal', req0=False):
	
		'''a constructor for the digit system class
		   arguments:
				base_matrix - rational square matrix M that serves as a basis for the digit system
				rem - type of remaindes for the residue set; default is 'normal', also posible 'reduced'
				req0 - is 0-vector required to be included into the digit set; default is False
		   returns:
				initialized instance of a digit system (M, D), in which every integer vector has a finite radix expansion
		'''
		self.rems = rem
		self.need0 = req0
		
		'''Base matrix and it's dimensions'''
		self.base = (Matrix(base_matrix)).change_ring(QQ)
		#to do: add check .base_ring == QQ
	
		self.dimr = self.base.nrows()
		self.dimc = self.base.ncols()
		self.dim  = min(self.dimr, self.dimc)

		#to do: add check dimr == dimc

		self.zero = vector(ZZ, [0]*self.dim)
		
		'''Prepare transform matrices T and Q for rotation form'''
		self._init_transforms_()
		
		'''LCM of a base matrix denomiators'''
		self._initcden_()
	
		'''Integral part of a base matrix, its Smith Normal Form and invariants'''
		self.A = (self.cden*self.base).change_ring(ZZ)
		self.detA = self.A.det()
		self.adjA = self.A.adjugate()
		
		self.S, self.U, self.V = self.A.smith_form()
		
		self.inv_U = det(self.U) * (self.U).adjugate()
		self.inv_V = det(self.V) * (self.V).adjugate(

		)
		
		'''Smith invariants and basis matrix L for auxiliary lattice'''
		self.divs = [s//gcd(self.cden, s) for s in self.S.diagonal()]
		self.L = self.inv_U * diagonal_matrix(ZZ, self.divs) * (self.inv_V)
		self.invL = self.L.inverse()
		
		'''Hash vectors and residue vectors mod L'''
		self._inithash_()
		self.hvecs = [vector(ZZ, tuple) for tuple in np.ndindex(*self.divs)]
		self.res = {self.vhash(vec): self.inv_U*self.redvec(vec) for vec in self.hvecs}
		
		
		#standard hulls for residue vectors
		self.hull = {h: self.enclose(self.L, self.res[h]) for h in self.res.keys()}
		
		#make a digit set
		self._build_digits_()
	
	def _initcden_(self):
	
		''' Finds LCM denominator of a base matrix entries '''
		
		self.cden  = 1
		for rw in self.base.rows():
			for el in rw:
				self.cden = lcm(self.cden, el.denominator())		
		
	def _inithash_(self):
		
		''' Initialized a sequence of factors for the hash function of vectors '''	
		
		self.hfacts = [1];
		for d in reversed(self.divs[1:]):
			self.hfacts.insert(0, d*self.hfacts[0])
			
	def realpart(self, vec):
		
		'''returns a real part of a complex vector'''
		
		return vector([real(t) for t in vec.list()])
	
	def imagpart(self, vec):
		
		'''returns a real part of a complex vector'''
		
		return vector([imag(t) for t in vec.list()])
		
	def unique(self, lst):
		
		'''Returns list of unique elements of a list'''
		
		unq = []
		for itm in lst:
			if not itm in unq:
				unq += [itm]
		return unq
		
	def _init_transforms_(self):
		
		'''Prepares transformation matrices that bring M into rotation form'''
		
		D, P = self.base.eigenmatrix_right()
		cols = P.columns()
		rcls = [self.realpart(c) for c in cols]
		ucols = self.unique(rcls)
		ncols = []
		for u in ucols:
			ncols += [u]
			v = self.imagpart(cols[rcls.index(u)])
			if v + v != v:
				ncols += [v]
		self.T = matrix(ncols).transpose()
		self.T.change_ring(AA)
		self.Tinv = self.T.inverse()
		
		self.Q = self.Tinv.transpose()*self.Tinv
		
	def invnorm(self, vec):
		
		'''Invariant under the multiplication by M norm of R^d'''
		
		return(norm(self.Tinv*vec))
		
	def redrem(a, b):
		r = a % b
		if (2*r > b):
			r -= b
		return r
	
	def redvec(self, vec):
		'''Reduces entries vec[i] of vec to min{r[i], r[i]-d[i]}, where x[i] = r[i] mod Smith invariant divs[i]'''
		rvec = vector(vec)
		if self.rems == 'reduced':
			for i in range(self.dim):
				rvec[i] = self.redrem(vec[i], self.divs[i])
		return rvec
			
		
	def hvec(self, vec):
		
		''' Returns the hash vector of a residue class mod L'''
		
		tvec = self.U * vec
		for i in range(self.dim):
			tvec[i] %= self.divs[i]
			
		return tvec

	def vhash(self, vec):
		
		'''Calculates hash code of a residue class mod Smith invariants'''
		
		return sum(vec[i]*self.hfacts[i] for i in range(self.dim))
	
	def hcode(self, vec):
		
		'''Hash code of a residue class mod auxiliary lattice L'''
		
		hv = self.hvec(vec)
		return self.vhash(hv)
	
	def resvec(self, vec):
	
		'''returns the residue vector mod lattice L'''
		return self.res[self.hcode(vec)]
		
	def enclose(self, L, x):
		
		'''Builds enclosing convex hull of x from lattice L points'''
		
		y = L.inverse()*x
		d = len(y)
		nvert = vector(ZZ, [t.round('toward') for t in y]) #nearest to y vertex of grid Z^d
		z = y - nvert
		signs = [1]*d		#quadrant location
		tip = [0]*d			#point of origin
		cone = identity_matrix(ZZ, d).columns()	#simplex cone
		for i in range(d):
			if z[i] < 0:
					signs[i] = -1
			else:
				if z[i] == 0:			   #if z is on the boundary, push tip of the cone backwards in that coordinate
					tip[i] = -1
		cone.insert(0, vector(ZZ,tip))     #insert the (modified) tip into cone
		sm = diagonal_matrix(ZZ, signs)
		
		return [L*(nvert+sm*vec) for vec in cone]
			
	def get_matrixform(self, V):
		
		'''prepares M, b for the linear programming problem M*x<=b in R^d'''
		
		M = matrix([self.Q*v  for v in V])
		b = vector([v*(self.Q*v)/2 for v in V])
		return M, b

	def get_polyhedron(self, M, b):
	
		'''solves linear programming problem M*x<=b in R^d'''
		
		lp = MixedIntegerLinearProgram(solver="PPL")
		x = lp.new_variable(integer=True); x
		lp.add_constraint(M*x <= b)
		lp.set_objective(None)
		return lp.polyhedron()

	def get_solutions(self, V, L, r):
		
		'''Calculates the feasible region of M*x<=b in R^d and all lattice points x = r (mod L)  in that region'''
		
		M, b = self.get_matrixform(V)
		feas = self.get_polyhedron(M, b) 										  #feasible region F
		sols = [L*x+r for x in self.get_polyhedron(M*L, b-M*r).integral_points()]  #lattice points r (mod L) inside F
		return 	(feas, sols)
	
	def nearest(self, x, V):
		'''Finds and element of V that is closest to x in M-invariant norm'''
		min_dist = oo
		for v in V:
			dist = self.invnorm(x-v)
			if dist < min_dist:
				min_dist = dist
				minv = v
		return minv
		
	def ndig(self,x):
		'''Returns the digit that is nearest to x in M-invariant norm and has same residue mod L'''
		return self.nearest(x, self.dig[self.hcode(x)])

	def F(self, x):
		'''Returns a pair (F(x), digit(x)), where digit(x)= x mod L and F(x)=M^{-1}(x-digit(x))'''
		cdig = self.ndig(x)                                        #current digit
		qvec = self.cden*(self.adjA*(x-cdig))                       #we want integer vector remain integer
		return (vector(ZZ, [v // self.detA for v in qvec]), cdig)

	def digfun(self, x):
		''' Modified digit function: if x in not in attractor, finds nearest digit with same residue mod L, otherwise returns x'''
		if x in self.alldigs:
			return x
		else:
			return self.ndig(x)
		    
	def Phi(self, x):
		'''If x is not in attractor, calculates F(x)=M^{-1}(x-digit(x)). Otherwise returns x.'''
		if x in self.alldigs:
			return (self.zero, x)
		else:
			return self.F(x)
		
	def orbit(self, x, bound):
		'''Builds an orbit of x under the mapping x->F(x)'''
		orb = [x]
		y = (self.F(x))[0]
		while (y not in orb) and (len(orb) < bound):
			orb.append(y)
			y = (self.F(y))[0]
		return (orb, y)

	def expand(self, x):
		'''Returns list of digits in a digital expansion of a vector in the lattice L'''
		qvec = x
		cdig = None
		diglist = []
		qvclist = []
		cont = True
		while cont:
			qvec, cdig = self.Phi(qvec)
			cont = qvec not in qvclist
			qvclist.append(qvec)
			diglist.append(cdig)
			cont = cont and (qvec != self.zero)
		return qvclist, diglist

	def digexp(self, x):
		'''Returns list of digit numbers in a digital expansion of a vector in the lattice L'''
		diglist = self.expand(x)[1]
		return [self.alldigs.index(cdig) for cdig in diglist]

	def assemble(self, diglist):
		s = vector(QQ, self.zero)
		#print('s=' +str(s))
		for dg in reversed(diglist):
			s = self.base*s + dg
			#print('digit=' + str(dg) + ', s=' +str(s))
		return s

	def digassm(self, idxlist):
		s = vector(QQ, self.zero)
		for idx in reversed(idxlist):
			s = self.base*s + self.alldigs[idx]
		return s
	      

	def _init_attractor_(self):
		'''Builds attractor set of the DS'''
		self.attractor = []
		self.attr_orbit = {}
		checked = []
		for x in self.allreps:
			if x not in checked:
				xorbit, last = self.orbit(x, testBound)
				self.attractor += [last]
				self.attr_orbit[tuple(last)] = xorbit[xorbit.index(last):]
				checked += xorbit
		self.attractor = self.unique(self.attractor)

		#removes attractor points whose orbit caintains the digit
		self.tidy_attractor = []

		for x in self.attractor:
			xlst= [y for y in self.attr_orbit[tuple(x)] if y in self.dig[self.hcode(y)]]
			if len(xlst) == 0:
				self.tidy_attractor.append(x)
		
		#identify all possible distinct digits
		self.alldigs = []
		for h in self.res.keys():
			self.alldigs += self.dig[h]
		self.alldigs += self.tidy_attractor
		if self.need0:
		    self.alldigs.append(self.zero)
		self.alldigs = self.unique(self.alldigs)

	def _build_digits_(self):
		'''Builds digit set'''
		self.dig = {}
		self.reps = {} #repeller points
		self.regs = {} #repelling regions
		self.allreps = []
		for h in self.res.keys():
			self.dig[h] = [self.res[h] - vec for vec in self.hull[h]]
			self.regs[h], self.reps[h] = self.get_solutions(self.dig[h], self.L, self.res[h])  
			self.allreps += self.reps[h]
		
		self._init_attractor_()

class TwistedDS:

	def __init__(self, A, B, C):

		self.baseA = matrix(A)
		self.baseA.inverse()
		self.dimA = self.baseA.ncols()
		
		self.baseB = matrix(B)
		self.baseBinv = self.baseB.inverse()
		self.dimB = self.baseB.ncols()
		
		self.baseC = matrix(ZZ, C)

		self.baseO = zero_matrix(ZZ, self.baseA.nrows(), self.baseB.ncols())

		#self.baseN = zero_matrix(ZZ, self.baseB.nrows, self.baseA.ncols)
		#self.baseN[0][-1]=1

		self.baseM = block_matrix([[self.baseA, self.baseO], [self.baseC, self.baseB]])
		self.baseMadj = self.baseM.adjugate()
		self.baseMinv = self.baseM.inverse() 

		self.rdsA = RotationDS(self.baseA, rem='normal', req0=True)
		self.rdsB = RotationDS(self.baseB, rem='normal', req0=True)
		
		self.zero = self.blocksum(self.rdsA.zero, self.rdsB.zero)
		self._build_digits_()


	def project(self, vec):
		return (vec[:self.dimA], vec[self.dimA:])

	def blocksum(self, v, w):
		return vector(v.list()+w.list())

	def digB(self, z):
		x, y = self.project(z)
		return self.rdsB.digfun(y-self.baseC*(self.rdsA.Phi(x)[0]))

	def PhiB(self, z):
		x, y = self.project(z)
		return self.rdsB.Phi(y-self.baseC*(self.rdsA.Phi(x)[0]))		
		

	def digM(self, z):
		return self.blocksum(self.rdsA.digfun(self.project(z)[0]), self.digB(z))

	def PhiM(self, z):
		return (self.blocksum(self.rdsA.Phi(self.project(z)[0])[0], self.PhiB(z)[0]), self.digM(z))

	def _build_digits_(self):
	    self.alldigs = [self.blocksum(dig1, dig2) for dig1 in self.rdsA.alldigs for dig2 in self.rdsB.alldigs]

#------ Provides orbit finding and digit expansion capabilities for Twisted Digit Functions class ------

	def orbit(self, x, bound):
		'''Builds an orbit of x under the mapping x->F(x)'''
		orb = [x]
		y = (self.PhiM(x))[0]
		while (y not in orb) and (len(orb) < bound):
			orb.append(y)
			y = (self.PhiM(y))[0]
		return (orb, y)

	def expand(self, x):
		'''Returns list of digits in a digital expansion of a vector in the lattice L'''
		qvec = x
		cdig = None
		diglist = []
		qvclist = []
		cont = True
		while cont:
			qvec, cdig = self.PhiM(qvec)
			cont = qvec not in qvclist
			qvclist.append(qvec)
			diglist.append(cdig)
			cont = cont and (qvec != self.zero)
		return qvclist, diglist

	def digexp(self, x):
		'''Returns list of digit numbers in a digital expansion of a vector in the lattice L'''
		diglist = self.expand(x)[1]
		return [self.alldigs.index(cdig) for cdig in diglist]

	def assemble(self, diglist):
		s = vector(QQ, self.zero)
		#print('s=' +str(s))
		for dg in reversed(diglist):
			s = self.baseM*s + dg
			#print('digit=' + str(dg) + ', s=' +str(s))
		return s

	def digassm(self, idxlist):
		s = vector(QQ, self.zero)
		for idx in reversed(idxlist):
			s = self.baseM*s + self.alldigs[idx]
		return s

