""" Kummer arithmetic without objects. Stolen from Michael Stoll's magma code.
   TODO: 
   (1) Lifting from Kummer to Jacobian. Need to be careful about precision
   (2) General surjection J-->P. Currently only for pairs of points
   (3) Take care of "is zero vs. is weakly zero" on local fields.
"""

#other comments on sage translation: 
#no "universe" in sage for a list, so atm just take parent of first list item
#and assert all other elements are in there..
#precision of the field I have replaced with precision of the element.


load("../QCExample/kummer_data.sage")


def auxf(f):
	#Not really sure yet what this function is doing but
	#Note that I have shortened the length of the output 
	#list since I am assuming curve of the form y^2 = quintic.
	f0 = f[0]
	f1 = f[1]
	f2 = f[2]
	f3 = f[3]
	f4 = f[4]
	f5 = f[5]
	P4.<k1,k2,k3,k4> = PolynomialRing(f.parent().base_ring())
	P2.<s1,s2> = PolynomialRing(f.parent().base_ring())
	return [f5*s1*s2^2 + 2*f4*s2^2 + f3*s1*s2 + 2*f2*s2 + f1*s1 + 2*f0,
       [f0,f1,f2,f3,f4,f5],
       k1^3*k4 + f2*k1^4 + f3*k1^3*k2 + f4*k1^2*k2^2 + f5*k1*k2*(k2^2-k1*k3),
       -k1^2*k2*k4 + f5*k1*k3*(k1*k3-2*k2^2)
        - 2*f4*k1^2*k2*k3 - f3*k1^3*k3 + f1*k1^4,
       k1^2*k3*k4 + f0*k1^4 + f4*k1^2*k3^2 + f5*k1*k2*k3^2] #5 outputs 


def Points(f,s0):
	#Return the indexed set of points on J mapping to P on
	#Kummer surface of J;
	#but instead of inputting J I am going to input f
	#and P already as a list.
	aux = auxf(f)
	R = f.parent().base_ring()
	PR.<z> = PolynomialRing(R)
	if s0[0]!=0:
		#point on J must have the form <x^2 - tx + n, a*x + b, 2>
		#with s1 = [1, t, n, c]
		s1 = [c/s0[0] for c in s0]
		pola = PR([s1[2], -s1[1], s1[0]])
		G0 = aux[2](s1)
		#FB:NB: s1 is normalised, 
		#under any reasonable assumption on the fi
		#then if working over a Laurent series ring
		#G0 has non-negative valuation, so can just work
		#over power series rings?
		#same for p2
		p1 = PR([-G0,0,1]) #a is a root of p1
		#Going to assume char(R) != 2.
		disc1 = p1.discriminant()
		if disc1 == 0: #unique value for a
			a = 0
			G1 = aux[4](s1)
			p2 = PR([-G1,0,1]) #b is a root of p2
			disc2 = p2.discriminant()
			if disc2 == 0: #unique value for b
				b = 0
				return [[pola, PR([b,a])]]
			else: #disc2 != 0
				test = disc2.is_square()
				if test: #two values for b
					re = 0
					im = disc2.sqrt()/2
					return [[pola, PR([im,a]), 2], [pola, PR([-im,a]),2]]
				else:
					return []
		else: #disc1!=0
			#FB: I used the following before realising I could
			#just work with power series directly.
			#if isinstance(disc1.parent(), LaurentSeriesRing):
				#print("valuation is", disc1.valuation())
				#disc1_val = disc1.valuation()
				#if disc1_val % 2 != 0:
					#return []
				#elif disc1_val >= 0:
					#disc1 = disc1.power_series()
					#test = disc1.is_square()
					#if test: #then two values for a
						#im = disc1.sqrt()/2
						#a1 = im
						#a2 = -im
						#G2 = aux[3](s1)
						#b1,b2 = tuple([G2/(2*a) for a in [a1,a2]]) #have 2*a*b - G2 = 0.
						#the denominator can't vanish because disc1 doesn't.
						#return [[pola, PR([b1,a1]),2], [pola,PR([b2,a2]),2]]
					#else:
						#return []
				#else:
					#t = disc1.parent().gens()[0]
					#disc1_norm = t^(-disc1_val)*disc1
					#test = disc1_norm.power_series().is_square()
					#if test: #then two values for a
						#im = disc1.sqrt()/2*t^(ZZ(disc_1_val/2))
						#a1 = im
						#a2 = -im
						#G2 = aux[3](s1)
						#b1,b2 = tuple([G2/(2*a) for a in [a1,a2]]) #have 2*a*b - G2 = 0.
						##the denominator can't vanish because disc1 doesn't.
						#return [[pola, PR([b1,a1]),2], [pola,PR([b2,a2]),2]]
					#else:
						#return []
			#else:
			test = disc1.is_square()
			if test: #then two values for a
				im = disc1.sqrt()/2
				a1 = im
				a2 = -im
				G2 = aux[3](s1)
				b1,b2 = tuple([G2/(2*a) for a in [a1,a2]]) #have 2*a*b - G2 = 0.
				#the denominator can't vanish because disc1 doesn't.
				return [[pola, PR([b1,a1]),2], [pola,PR([b2,a2]),2]]
			else:
				return []
	
	elif s0[1] !=0:
		#point has form <x-t, a*x^3+b, 1> with s1 = [0, 1, t, c]
		t = s0[2]/s0[1]
		pola = PR(-t,1)
		p1 = PR([0,0,1])
		#in our case the discriminant of p1 will always be 0
		p2 = PR([-PR(aux[1])(t), 0, 1]) #y^2 - f(t)
		disc2 = p2.discriminant()
		if disc2 == 0:#unique value for b #could probs simplify: just ask for constant term to be 0
			return [[pola, 0]]
		else: 
			test = disc2.is_square()
			if test: #two values for b
				im = disc.sqrt()/2
				return [[pola, im], [pola, -im]]
			else:
				return []	

	else:
		return [0]

def normalize(xs):
	# Choose right normalization method
	# based on ground field.
	F = xs[0].parent()
	#FB: added next line - as if for instance you have a power series
	#want to normalise as a Laurent series.
	F = FractionField(F)
	assert [xs[i] in F for i in range(len(xs))], "look back at this"
	if isinstance(F,LaurentSeriesRing) or isinstance(F,sage.rings.padics.local_generic.LocalGeneric):
		#normalise such that coordinates are integral, with one being a unit,
		#and such that the first unit is 1.
		#First step: Make minimal valuation zero.
		v = min([min(x.valuation(), x.precision_absolute()) for x in xs if x!=0])
		a = F.uniformizer()^(-v)
		xs1 = [a*x for x in xs]
		# Second step: normalise first unit.
		i = 0
		while xs1[i].valuation() > 0:
			i += 1
		return [x/xs1[i] for x in xs1]
	i = 0
	while xs[i] == 0:
		i += 1
	return [x/xs[i] for x in xs]		
		

def double(xs, deltas):
	# Returns the double of the point P
	s = [d(xs) for d in deltas]
	return normalize(s)


def select_coordinate(xs):
	F = xs[0].parent()
	assert [xs[i] in F for i in range(len(xs))], "look back at this"
	if isinstance(F,LaurentSeriesRing) or isinstance(F,sage.rings.padics.local_generic.LocalGeneric):
		i = 0
		while xs[i].valuation() > 0:
			i += 1
	else:
		i = 0
		while xs[i] == 0:
			i += 1
	return i	


def pseudo_add(xP, xQ, xPminusQ, bbmatrix):
	#Given the images on the Kummer surface of points P, Q, P-Q on the
	#Jacobian, returns the image of P+Q
	xPQ = xP + xQ
	i = select_coordinate(xPminusQ)
	c1 = xPminusQ[i]
	c2 = bbmatrix[i][i](xPQ)
	L = [c1*bbmatrix[i][j](xPQ) - xPminusQ[j]*c2 for j in range(0,4)]
	L[i] = c1*c2
	return normalize(L)


def scalar_multiple(xs, n, deltas, bbmatrix):
	#The n-th multiple of P on the Kummer surface K
	P = xs
	F = xs[0].parent()
	assert [xs[i] in F for i in range(len(xs))], "look back at this"
	Zero = [F(0),F(0),F(0),F(1)]
	if n == 0:
		return Zero
	m = abs(n)
	Px = Zero
	Py = xs
	Pz = xs
	# invariants: Px = ax*P, Py = ay*P, Pz = az*P with
	# ax + ay = az and ax + m*az = n .
	while True:
		if m % 2 == 1:
			Px = pseudo_add(Px, Pz, Py, bbmatrix)
		else:
			Py = pseudo_add(Py, Pz, Px, bbmatrix)
		m = m // 2
		if m == 0:
			return Px
		Pz = double(Pz, deltas)


def kummer_point(P, Q, f):
	# Image of {P,Q} on Kummer. Here P and Q are points on y^2=f(x), given as pairs (or
	# triples(s), if infinity). assumes odd degree, but this isn't crucial.
	# TODO: Use Mumford rep as input.
	assert f.degree() == 5
	f0,f1,f2,f3,f4,f5 = tuple(f.coefficients(sparse = False))
	F = (P[0]).parent() #FB: might need to change
	if len(P) == len(Q) and P[0] == Q[0] and P[1] == -Q[1]:
		return [F(0),F(0),F(0),F(1)]
	#TODO: Watch out for p-adic precision issues in checking equality!
	if P == [0,1,0]:#swap P and Q #FB in sage point at infinity is 0,1,0
		R = P
		P = Q
		Q = R
	#Now P is affine
	x1,y1 = tuple(P)
	if Q == [0,1,0]:
		return normalize([0,1,x1,f5*x1^2])
	#Q is affine
	x2,y2 = tuple(Q)
	k1 = 1
	k2 = x1+x2
	k3 = x1*x2
	if P!=Q:
		k4 = (2*f0+f1*(x1+x2)+2*f2*x1*x2+f3*(x1+x2)*x1*x2+2*f4*x1^2*x2^2+f5*(x1+x2)*x1^2*x2^2-2*y1*y2)/(x1-x2)^2
	else:
		PF.<z> = PolynomialRing(F)
		eqn = kummer_equation(f, [k1,k2,k3,z])
		#FB: root finding does not work (in Sage) for polynomials
		#over Laurent series ring.
		#However, in this case, the polynomial is linear
		#so replace commented lines with root finding by -const term/coeff of z.
		#rts = eqn.roots()
		#assert len(rts) == 1
		#k4 = rts[0][0]
		assert eqn.degree() == 1
		k4 = -eqn[0]/eqn[1]
	return normalize([k1,k2,k3,k4])		


def kummer_data(f):
	F = f.base_ring()
	A4.<x1,x2,x3,x4> = PolynomialRing(F)
	xs = [x1,x2,x3,x4]
	eqn = kummer_equation(f,xs)
	deltas = kummer_deltas(f,xs)
	A8.<x1,x2,x3,x4,y1,y2,y3,y4> = PolynomialRing(F)
	xs = [x1,x2,x3,x4]
	ys = [y1,y2,y3,y4]
	bbmatrix = kummer_bbmatrix(f, xs, ys)
	return (eqn, deltas, bbmatrix)


"""
PQ.<x> = PolynomialRing(Rationals())
f = x^5-x^4+x^3+x^2- 2*x + 1
C = HyperellipticCurve(f)
P = C(0,1)
wP = C(0,-1)
O = C(0,1,0)
Q = C(-1,-1)
wQ = C(-1,1)
R = C(2,5)
wR = C(2,-5)
S = C(1,-1)
wS = C(1,1)


eqn, deltas, bbmatrix =kummer_data(f);
p = [P[0],P[1]]
q = [Q[0],Q[1]]
r = [R[0],R[1]]
s = [S[0],S[1]]
A = kummer_point(p,q,f)
B = kummer_point(r,s,f)
a2 = double(A,deltas)
a3 = scalar_multiple(A,3,deltas,bbmatrix);

"""

