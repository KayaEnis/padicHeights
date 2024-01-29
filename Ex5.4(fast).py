##########################################################################################

# Same as Ex5.4.py, except that we use the code Gajović and the third-named author instead 
# of Balakrishnan's code for computing Coleman integrals of differentials of the third kind.

# SETUP

_.<x> = PolynomialRing(QQ);
bx = x^5 - 4*x^4 - 48*x^3 + 64*x^2 + 512*x + 256 # monic version of [LMFDB,125237.a.125237.1]
b5, b4, b3, b2, b1, b0 = bx.list() # bx = x^5 + b1*x^4 + b2*x^3 + b3*x^2 + b4*x + b5
C = HyperellipticCurve(bx)
cup_prod = C.cup_product_matrix()

pt1, pt2   = C(-4,16), C(-4,-16) # know points (LMFDB)
pt3, pt4   = C(-3,5), C(-3,-5)
pt5, pt6   = C(0,16), C(0,-16)
pt7, pt8   = C(4,16), C(4,-16)
pt9, pt10  = C(8,16), C(8,-16)
pt11, pt12 = C(16,784), C(16,-784)
pt13, pt14 = C(36,7184), C(36,-7184)
badpoints = [pt11,pt12] # points in Weierstrass discs
nicepoints = [pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10,pt13,pt14]  # points in non-Weierstrass discs
points = badpoints + nicepoints # all in the component v0 = π⁻¹(U0)

N, rang = 16, 4 # truncation level and (2*genus)
p, prec = 7, 8 # reduction at p is an elliptic curve with an ordinary double point
k = GF(p)
K =  Qp(p,prec)
def Log(z): return K(z).log(p_branch=0, change_frac=True) # Iwasawa branch

# factorization of bx over K
bxK = bx.change_ring(K)
f1 = 1
for factor in list(bxK.factor()):
    factor_wo_mult = factor[0]
    if factor_wo_mult.change_ring(k) == (x^2+4*x+4).change_ring(k):
        f0 = factor_wo_mult
    else:
        f1 *= factor_wo_mult
# bxK = f0*f1 and f0 = (x-5)^2 , f1 = (x-2)(x^2+x+3) (modp)

def coeff(f):
    B = f(0)
    A = f(1) - B - 1
    C = B - A^2/4
    return K(A), K(B), K(C)
A0, B0, C0 = coeff(f0) # f0 = x^2 + A0*x + B0, C0  = B0 - A0^2/4

##########################################################################################

# A) The component v0 = π⁻¹(U0)

Ctilde = HyperellipticCurve(f1) # elliptic curve over K with good ordinary reduction
pole, wpole = Ctilde.lift_x(-A0/2, all = True) # points on Ctilde with x-coord's -A0/2
f1der = f1.derivative()

def NC0(z): # points in the new coordinates
    x, y = K(z[0]), K(z[1])
    Z = x + A0/2
    l = Z*K(1+C0/Z^2).sqrt()
    return Ctilde(x,y/l)

_.<t> = PowerSeriesRing(K, 't', default_prec = 2*N)
l = t*(1+C0*t^2)^(-1/2)

# 1) omega_i, i=0,1,2,3

# pole reduction: x-coordinate is -A0/2
x,y = Ctilde.local_coord(pole,2*N)
dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
w = [l.truncate(N)(1/(x+A0/2))*x^i*dx_2y for i in range(rang)] # original forms around the pole
d0 = [w[i].residue()/(1/(x+A0/2)*dx_y).residue() for i in range(rang)] # coefficients of 1/(x+A0/2)*dx/y
w = [w[i] - d0[i]/(x+A0/2)*dx_y for i in range(rang)]
ED = [0,0] + [((x+A0/2)*f1der(x)-2*i*f1(x))/(x+A0/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
_.<u,v> = QQ[]
F0 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F0[i] += FF[i][j]*u^(N-2-j)*v

# omega_3 requires a suitable multiple of dy
x, y = Ctilde.local_coordinates_at_infinity(2*N)
dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
F03new = F0[3] + (l.truncate(N)(1/(x+A0/2))*x^3*dx_2y - F0[3](1/(x+A0/2),y).derivative() - d0[3]/(x+A0/2)*dx_y).list()[0]/(y.derivative()).list()[0]*v
F0.append(F03new)
F0.remove(F0[3]) # final meromorphic functions

c01 = [(l.truncate(N)(1/(x+A0/2))*x^i*dx_2y - F0[i](1/(x+A0/2),y).derivative() - d0[i]/(x+A0/2)*dx_y).list()[0]/(x*dx_2y).list()[0] for i in range(rang)] # coefficients of x*dx/2y
c00 = [(l.truncate(N)(1/(x+A0/2))*x^i*dx_2y - F0[i](1/(x+A0/2),y).derivative() - d0[i]/(x+A0/2)*dx_y - c01[i]*x*dx_2y).list()[0]/(dx_2y).list()[0] for i in range(rang)] # coefficients of dx/2y

# omega_i = dF0[i] + d0[i]/(x+A0/2)*dx/y + c00[i]*dx/2y + c01[i]*x*dx/2y

def Int01(i,z1,z2): # integral of omega_i from z1 to z2
    S, R = NC0(z1), NC0(z2)
    x1, y1 = S[0], S[1]
    x2, y2 = R[0], R[1]
    exact_part = F0[i](1/(x2+A0/2),y2) - F0[i](1/(x1+A0/2),y1)
    basis_values = Ctilde.coleman_integrals_on_basis(S,R)
    second_kind_part = c00[i]*basis_values[0] + c01[i]*basis_values[1]
    third_kind_part = d0[i]/pole[1]*Coleman_Integral_third_kind_antisymmetric(pole,R,S)
    return exact_part + second_kind_part + third_kind_part

# 2) omega = b/(x-a)*dx/y where pt = (a,b) is in v0

def base00(pt):
    btilde0 = NC0(pt)[1] # coefficient of 1/(x-a)*dx/y
    return btilde0

def base01(pt):
    a, b = pt[0], pt[1]
    btilde0 = base00(pt)
    x, y = Ctilde.local_coord(pole,2*N)
    dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
    w = l.truncate(N)(1/(x+A0/2))*b/(x-a)*dx_y - btilde0/(x-a)*dx_y
    dpt0 = w.residue()/(1/(x+A0/2)*dx_y).residue() # coefficient of 1/(x+A0/2)*dx/y
    w = w - dpt0/(x+A0/2)*dx_y
    ED = [0,0] + [((x+A0/2)*f1der(x)-2*i*f1(x))/(x+A0/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
    FF = []
    while w.valuation() < 0:
        val = w.valuation()
        FF.append(w.list()[0]/ED[-val].list()[0])
        w = w - (w.list()[0]/ED[-val].list()[0])*ED[-val]
    _.<u,v> = QQ[]
    Fpt0 = 0  # meromorphic function
    for i in range(N-2):
        Fpt0 += FF[i]*u^(N-2-i)*v
    return [dpt0,Fpt0]

def base02(pt):
    a, b = pt[0], pt[1]
    btilde0 = base00(pt)
    dpt0,Fpt0 = base01(pt)
    x, y = Ctilde.local_coordinates_at_infinity(2*N)
    dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
    w = l.truncate(N)(1/(x+A0/2))*b/(x-a)*dx_y - btilde0/(x-a)*dx_y - dpt0/(x+A0/2)*dx_y - Fpt0(1/(x+A0/2),y).derivative()
    cpt01 = w.list()[0]/(x*dx_2y).list()[0] # coefficient of x*dx/2y
    cpt00 = (w - cpt01*x*dx_2y).list()[0]/(dx_2y).list()[0] # coefficient of dx/2y
    return [cpt00,cpt01]

def base0(pt):
    btilde0 = base00(pt)
    dpt0,Fpt0 = base01(pt)
    cpt00,cpt01 = base02(pt)
    return [btilde0,dpt0,Fpt0,cpt00,cpt01]

# omega = btilde0/(x-a)*dx/y + dpt0/(x+A0/2)*dx/y + dFpt0 + cpt00*dx/2y + cpt01*x*dx/2y

def Int02(pt,z1,z2): # integral of b/(x-a)*dx/y from z1 to z2. z1 and z2 must be DIFFERENT from pt and w(pt)
    a = pt[0]
    S, R = NC0(z1), NC0(z2)
    x1, y1 = S[0], S[1]
    x2, y2 = R[0], R[1]
    btilde0,dpt0,Fpt0,cpt00,cpt01 = base0(pt)
    exact_part = Fpt0(1/(x2+A0/2),y2) - Fpt0(1/(x1+A0/2),y1)
    basis_values = Ctilde.coleman_integrals_on_basis(S,R)
    second_kind_part = cpt00*basis_values[0] + cpt01*basis_values[1]
    third_kind_part_1 = dpt0/pole[1]*Coleman_Integral_third_kind_antisymmetric(pole,R,S)
    third_kind_part_2 = Coleman_Integral_third_kind_antisymmetric(Ctilde(a,btilde0),R,S)
    return exact_part + second_kind_part + third_kind_part_1 + third_kind_part_2

def Int03(pt,z1,z2): # integral of omega = (y+b)/(x-a)*dx/2y from z1 to z2. z1 and z2 must be DIFFERENT from pt and w(pt) (I can easily remove this)
    a = pt[0]
    x1 = NC0(z1)[0]
    x2 = NC0(z2)[0]
    first_term = Log(x2-a) - Log(x1-a)
    second_term = Int02(pt,z1,z2)
    return 1/2*(first_term + second_term)

##########################################################################################

# B) LOCAL HEIGHT PAIRING AT P

# Blakestad's constants
# computed using canonical_constants in canonical_subspace.py
alpha_Bl = 3*7 + 2*7^2 + 2*7^4 + 3*7^5 + O(7^6)
beta_Bl  = 5 + 5*7^2 + 2*7^3 + 3*7^4 + 5*7^5 + O(7^6)
gamma_Bl = 7 + 3*7^2 + 2*7^3 + 5*7^4 + 5*7^5 + O(7^6)
delta_Bl = 7 + 7^2 + 3*7^3 + 2*7^5 + O(7^6)

# W = Blakestad's subspace is generated by eta0 and eta1, where
# eta0 = k0*omega0 + k1*omega1 + k2*omega2
# eta1 = l0*omega0 + l1*omega1 + l2*omega2 + l3*omega3
k0 = K(-(b2 - b1*beta_Bl + 3*gamma_Bl))
k1 = K(-beta_Bl)
k2 = K(1)
l0 = K(-(3*b1*b2 - b1*alpha_Bl + 3*b3 + 3*delta_Bl))
l1 = K(-alpha_Bl)
l2 = K(3*b1)
l3 = K(3)
# W is isotropic: if N = cup_prod, then eta0 ∪ eta1 =
# k0*(l0*N[0][0] + l1*N[1][0] + l2*N[2][0] + l3*N[3][0]) +
# k1*(l0*N[0][1] + l1*N[1][1] + l2*N[2][1] + l3*N[3][1]) +
# k2*(l0*N[0][2] + l1*N[1][2] + l2*N[2][2] + l3*N[3][2]) = O(7^6)

# z1 = (a,b), z2 = (c,d), D = (z1)-(z2) ==> omega = ((y+b)/(x-a)-(y+d)/(x-c))*dx/2y
# is a form of the third kind such that Res(omega) = D.
# Ψ(omega) = d0*omega0 + d1*omega1 + e0*eta0 + e1*eta1 for some d0,d1,e0,d1.
# omegaD := omega - d0*omega0 - d1*omega1 ==> omegaD is the unique form
# of the third kind such that Res(omegaD) = D and Ψ(omegaD) is in W.

def omegaD_coeff(z1,z2):
    global_symbols = Matrix([[Int01(0,z1,z2)],[Int01(1,z1,z2)],[Int01(2,z1,z2)],[Int01(3,z1,z2)]])
    c0 = (-cup_prod.inverse()*global_symbols)[0][0]
    c1 = (-cup_prod.inverse()*global_symbols)[1][0]
    c2 = (-cup_prod.inverse()*global_symbols)[2][0]
    c3 = (-cup_prod.inverse()*global_symbols)[3][0]
    e0 = c3/l3
    e1 = (c2-e0*l2)/k2
    d0 = c0 - e0*l0 - e1*k0
    d1 = c1 - e0*l1 - e1*k1
    return [d0,d1]

def local_height_at_p(z1,z2,z3,z4): # h_p((z1)-(z2),(z3)-(z4))
    d0,d1 = omegaD_coeff(z1,z2)
    return Int03(z1,z4,z3) - Int03(z2,z4,z3) - d0*Int01(0,z4,z3) - d1*Int01(1,z4,z3)

##########################################################################################

# C) EXAMPLE

P, wP = pt1, pt2
Q     = pt4
R, wR = pt5, pt6
S     = pt7
T, wT = pt9, pt10
U, wU = pt13, pt14

print('------------------------------------------------------------------------------------')
print('local height of (P)-(wP) and (R)-(wR) at p is ' + str(local_height_at_p(P,wP,R,wR)))
print('------------------------------------------------------------------------------------')
print('local height of (Q)-(S) and (T)-(U) at p is '   + str(local_height_at_p(Q,S,T,U)))
print('------------------------------------------------------------------------------------')
print('local height of (Q)-(S) and (wT)-(wU) at p is ' + str(local_height_at_p(Q,S,wT,wU)))
print('------------------------------------------------------------------------------------')

##########################################################################################
