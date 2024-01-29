################################################################################

# SETUP

_.<x> = PolynomialRing(QQ)
bx = x^5 + 5*x^4 - 168*x^3 + 1584*x^2 - 10368*x + 20736 # from [BMS16]
b5, b4, b3, b2, b1, b0 = bx.list() # bx = x^5 + b1*x^4 + b2*x^3 + b3*x^2 + b4*x + b5
C = HyperellipticCurve(bx)
cup_prod = C.cup_product_matrix()

pt1, pt2   = C(-12,720), C(-12,-720)
pt3, pt4   = C(-8,528), C(-8,-528)
pt5, pt6   = C(0,144), C(0,-144)
pt7, pt8   = C(8,80), C(8,-80)  
pt9, pt10  = C(12,432), C(12,-432)
pt11, pt12 = C(36,7920), C(36,-7920)
points0 = [pt3,pt4,pt5,pt6,pt9,pt10,pt11,pt12] # all in the component v0 = π⁻¹(U0)
points1 = [pt1,pt2,pt7,pt8] # all in the component v1 = π⁻¹(U1)
points = points0 + points1

N, rang = 16, 4 # truncation level and (2*genus)
p, prec = 5, 9 # reduction at p is a projective line with two ordinary double points
k     = GF(p)
K     = Qp(p,prec)
L.<a> = K.ext(x^2-2) # unramified extension
M.<b> = L.ext(x^2-p) # ramified extension defined by x^2 - p

def Log(z): # Iwasawa branch
    if z == 0:
        return "Log is not defined at 0"
    elif z in L:
        Lz = L(z)
        val = Lz.valuation()
        return log(Lz/p^val)
    else:
        Mz = M(z)
        val = Mz.valuation()
        return log(Mz/b^val)

# factorization of bx over K
f2 = x^2 + 12*x - 144
h = x^3 - 7*x^2 + 60*x - 144 # b/f2
roots_of_h = h.change_ring(K).roots(multiplicities=False)
for root in roots_of_h:
    if k(root) == 1:
        alpha = root
f0 = x - alpha
roots_of_h.remove(alpha)
f1 = prod(x-root for root in roots_of_h)
# bx = f0*f1*f2, and f0 = x-1, f1 = (x-3)^2, f2 = (x-4)^2 mod p

def coeff(f):
    B = f(0)
    A = f(1) - B - 1
    C = B - A^2/4
    return K(A), K(B), K(C)
A1, B1, C1 = coeff(f1) # f1 = x^2 + A1*x + B1, C1  = B1 - A1^2/4
A2, B2, C2 = coeff(f2) # f2 = x^2 + A2*x + B2, C2  = B2 - A2^2/4

##########################################################################################

# A) THE COMPONENT v0 = π⁻¹(U0)

f0der = f0.derivative()
def NC0(z): # points in the new coordinates
    x, y = z[0], z[1]
    Z1 = x + A1/2
    l1 = Z1*(1+C1/Z1^2).sqrt()
    Z2 = x + A2/2
    l2 = Z2*(1+C2/Z2^2).sqrt()
    l = l1*l2
    return (x,y/l)

_.<t> = PowerSeriesRing(K, 't', default_prec = N)
l01 = t*(1+C1*t^2)^(-1/2)
l02 = t*(1+C2*t^2)^(-1/2)
def ell0(x): # series in the pullback: omega -- > ell0(x)*omegatilde
    return l01.truncate()(1/(x+A1/2))*l02.truncate()(1/(x+A2/2))

# 1) omega_i, i=0,1,2,3

# first pole reduction: x-coordinate is -A1/2
_.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
x = t - A1/2
y = f0(x).sqrt()
dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
w = [ell0(x)*x^i*dx_2y for i in range(rang)] # original forms around the first pole
d01 = [w[i].residue()/(1/(x+A1/2)*dx_y).residue() for i in range(rang)] # coefficients of 1/(x+A1/2)*dx/y
w = [w[i] - d01[i]/(x+A1/2)*dx_y for i in range(rang)]
ED = [0,0] + [((x+A1/2)*f0der(x)-2*i*f0(x))/(x+A1/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
_.<u,v> = QQ[]
F01 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F01[i] += FF[i][j]*u^(N-2-j)*v

# second pole reduction: x-coordinate is -A2/2
_.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
x = t - A2/2
y = f0(x).sqrt()
dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
w = [ell0(x)*x^i*dx_2y - F01[i](1/(x+A1/2),y).derivative() - d01[i]/(x+A1/2)*dx_y for i in range(rang)]
d02 = [w[i].residue()/(1/(x+A2/2)*dx_y).residue() for i in range(rang)] # coefficients of 1/(x+A2/2)*dx/y
w = [w[i] - d02[i]/(x+A2/2)*dx_y for i in range(rang)]
ED = [0,0] + [((x+A2/2)*f0der(x)-2*i*f0(x))/(x+A2/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
FF = [[] for i in range(rang)]
for i in range(rang):
    while w[i].valuation() < 0:
        val = w[i].valuation()
        FF[i].append(w[i].list()[0]/ED[-val].list()[0])
        w[i] = w[i] - (w[i].list()[0]/ED[-val].list()[0])*ED[-val]
_.<u,v> = QQ[]
F02 = [0 for i in range(rang)] # meromorphic functions
for i in range(rang):
    for j in range(N-2):
        F02[i] += FF[i][j]*u^(N-2-j)*v

dinf020 = ((ell0(x)*x^2*dx_2y - F01[2](1/(x+A1/2),y).derivative() - d01[2]/(x+A1/2)*dx_y - F02[2](1/(x+A2/2),y).derivative() - d02[2]/(x+A2/2)*dx_y)/(dx_2y)).list()[0]
dinf031 = ((ell0(x)*x^3*dx_2y - F01[3](1/(x+A1/2),y).derivative() - d01[3]/(x+A1/2)*dx_y - F02[3](1/(x+A2/2),y).derivative() - d02[3]/(x+A2/2)*dx_y)/(dx_2y)).list()[1]
dinf030 = ((ell0(x)*x^3*dx_2y - F01[3](1/(x+A1/2),y).derivative() - d01[3]/(x+A1/2)*dx_y - F02[3](1/(x+A2/2),y).derivative() - d02[3]/(x+A2/2)*dx_y)/(dx_2y) - dinf031*x).list()[0]

# omega_i = dF01[i] + d01[i]/(x+A1/2)*dx/y + dF02[i] + d02[i]/(x+A2/2)*dx/y, i = 0,1
# omega_2 = dF01[2] + d01[2]/(x+A1/2)*dx/y + dF02[2] + d02[2]/(x+A2/2)*dx/y + dinf020*dx/2y
# omega_3 = dF01[3] + d01[3]/(x+A1/2)*dx/y + dF02[3] + d02[3]/(x+A2/2)*dx/y + dinf030*dx/2y + dinf031*x*dx/2y

def Int01(i,z1,z2): # integral of omega_i from z1 to z2
    x1, y1 = NC0(z1)
    x2, y2 = NC0(z2)
    exact_part_1 = F01[i](1/(x2+A1/2),y2) - F01[i](1/(x1+A1/2),y1)
    exact_part_2 = F02[i](1/(x2+A2/2),y2) - F02[i](1/(x1+A2/2),y1)
    exact_part = exact_part_1 + exact_part_2
    D1 = L(-alpha-A1/2).sqrt()
    D2 = L(-alpha-A2/2).sqrt()
    third_kind_part_1 = d01[i]/D1*(Log((y2-D1)/(y2+D1)) - Log((y1-D1)/(y1+D1)))
    third_kind_part_2 = d02[i]/D2*(Log((y2-D2)/(y2+D2)) - Log((y1-D2)/(y1+D2)))
    third_kind_part = third_kind_part_1 + third_kind_part_2
    if i == 2:
        inf_part_0 = dinf020*(y2 - y1)
        inf_part_1 = 0
    elif i == 3:
        inf_part_0 = dinf030*(y2 - y1)
        inf_part_1 = dinf031*((y2^3/3 + alpha*y2) - (y1^3/3 + alpha*y1))
    else:
        inf_part_0 = 0
        inf_part_1 = 0
    inf_part = inf_part_0 + inf_part_1
    return exact_part + third_kind_part + inf_part

# 2) omega = b/(x-a)*dx/y where pt = (a,b) is in v0

def base00(pt):
    btilde0 = NC0(pt)[1] # coefficient of 1/(x-a)*dx/y
    return btilde0

def base01(pt):
    a, b = pt[0], pt[1]
    btilde0 = base00(pt)
    _.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
    x = t - A1/2
    y = f0(x).sqrt()
    dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
    w = ell0(x)*b/(x-a)*dx_y - btilde0/(x-a)*dx_y
    dpt01 = w.residue()/(1/(x+A1/2)*dx_y).residue() # coefficient of 1/(x+A1/2)*dx/y
    w = w - dpt01/(x+A1/2)*dx_y
    ED = [0,0] + [((x+A1/2)*f0der(x)-2*i*f0(x))/(x+A1/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
    FF = []
    while w.valuation() < 0:
        val = w.valuation()
        FF.append(w.list()[0]/ED[-val].list()[0])
        w = w - (w.list()[0]/ED[-val].list()[0])*ED[-val]
    _.<u,v> = QQ[]
    Fpt01 = 0  # meromorphic function
    for i in range(N-2):
        Fpt01 += FF[i]*u^(N-2-i)*v
    return [dpt01,Fpt01]

def base02(pt):
    a, b = pt[0], pt[1]
    btilde0 = base00(pt)
    dpt01,Fpt01 = base01(pt)
    _.<t> = PowerSeriesRing(L, 't', default_prec = 2*N)
    x = t - A2/2
    y = f0(x).sqrt()
    dx_y, dx_2y = x.derivative()/y, x.derivative()/(2*y)
    w = ell0(x)*b/(x-a)*dx_y - btilde0/(x-a)*dx_y - dpt01/(x+A1/2)*dx_y - Fpt01(x+A1/2,y).derivative()
    dpt02 = w.residue()/(1/(x+A2/2)*dx_y).residue() # coefficient of 1/(x+A2/2)*dx/y
    w = w - dpt02/(x+A2/2)*dx_y
    ED = [0,0] + [((x+A2/2)*f0der(x)-2*i*f0(x))/(x+A2/2)^(i+1)*dx_2y for i in range(1,N-1)] # exact differentials
    FF = []
    while w.valuation() < 0:
        val = w.valuation()
        FF.append(w.list()[0]/ED[-val].list()[0])
        w = w - (w.list()[0]/ED[-val].list()[0])*ED[-val]
    Fpt02 = 0  # meromorphic function
    for i in range(N-2):
        Fpt02 += FF[i]*u^(N-2-i)*v
    return [dpt02,Fpt02]

def base0(pt):
    btilde0 = base00(pt)
    dpt01,Fpt01 = base01(pt)
    dpt02,Fpt02 = base02(pt)
    return [btilde0,dpt01,Fpt01,dpt02,Fpt02]

# omega = btilde0/(x-a)*dx/y + dpt01/(x+A1/2)*dx/2y + dFpt01 + dpt02/(x+A2/2)*dx/2y + dFpt02

def Int02(pt,z1,z2): # integral of omega = (y+b)/(x-a)*dx/2y from z1 to z2
    a, b = pt[0], pt[1]
    btilde0,dpt01,Fpt01,dpt02,Fpt02 = base0(pt)
    x1, y1 = NC0(z1)
    x2, y2 = NC0(z2)
    first_part = 2*(Log(y2 - btilde0) - Log(y1 - btilde0))
    exact_part_1 = Fpt01(1/(x2+A1/2),y2) - Fpt01(1/(x1+A1/2),y1)
    exact_part_2 = Fpt02(1/(x2+A2/2),y2) - Fpt02(1/(x1+A2/2),y1)
    exact_part = exact_part_1 + exact_part_2
    D1 = L(-alpha-A1/2).sqrt()
    D2 = L(-alpha-A2/2).sqrt()
    third_kind_part_1 = dpt01/D1*(Log((y2-D1)/(y2+D1)) - Log((y1-D1)/(y1+D1)))
    third_kind_part_2 = dpt02/D2*(Log((y2-D2)/(y2+D2)) - Log((y1-D2)/(y1+D2)))
    third_kind_part = third_kind_part_1 + third_kind_part_2
    second_part = exact_part + third_kind_part
    return 1/2*(first_part + second_part)

##########################################################################################

# B) THE COMPONENT v1 = π⁻¹(U1) (the following is based on the parametrization of v1 as an annulus)

def NC1(z): # points in the new coordinates
    x, y = z[0], z[1]
    l0 = (-alpha).sqrt()*M(1-x/alpha).sqrt()
    Z2 = x + A2/2
    l2 = Z2*(1+C2/Z2^2).sqrt()
    l = l0*l2
    return (x,y/l)

# 1) omega_i, i=0,1,2,3

r1,r2 = f1.roots(multiplicities=False) # f1 = (x-r1)*(x-r2)
D = (r1 - r2)/2
_.<t> = LaurentSeriesRing(L, 't')
l10 = 0
for i in range(N):
    l10 += binomial(-1/2,i)*((t-D)^2/(2*(r1-alpha)*t))^i
l10 = (L(r1-alpha).sqrt())^(-1)*l10
const = r1^2 + 12*r1 - 144 # (r1-theta1)*(r1-theta2) where f2 = (x-theta1)(x-theta2)
l12 = 0
for i in range(N/2):
    l12 += binomial(-1/2,i)*((2*r1+12)/const*(t-D)^2/(2*t) + (t-D)^4/(4*const*t^2))^i
l12 = (const.sqrt())^(-1)*l12
omega = [l10*l12*((t-D)^2/(2*t)+r1)^i*(1/(2*t)) for i in range(rang)]
res = [omega[i].residue() for i in range(rang)]
omega_wo_res = [omega[i] - res[i]*t^-1 for i in range(rang)]

def Int11(i,z1,z2): # integral of omega_i from z1 to z2
    x1, y1 = NC1(z1)
    x2, y2 = NC1(z2)
    t1 = x1 + y1 - (r1+r2)/2
    t2 = x2 + y2 - (r1+r2)/2
    int_wo_res = omega_wo_res[i].integral()
    return int_wo_res(t2) - int_wo_res(t1) + res[i]*(Log(t2) - Log(t1))

# 2) omega = b/(x-a)*dx/y where pt = (a,b) is NOT in v1

def wab(pt):
    a, b = pt[0], pt[1]
    wab = 0
    for i in range(N):
        wab += (-1)^i*((t-D)^2/(2*(r1-a)*t))^i
    wab = b/(r1-a)*wab # b/(x-a) in terms of t
    return wab

def Int12(pt,z1,z2): # integral of omega = b/(x-a)*dx/y from z1 to z2
    x1, y1 = NC1(z1)
    x2, y2 = NC1(z2)
    t1 = x1 + y1 - (r1+r2)/2
    t2 = x2 + y2 - (r1+r2)/2
    omega = l10*l12*wab(pt)*(1/t)
    res = omega.residue()
    omega_wo_res = omega - res*t^-1
    int_wo_res = omega_wo_res.integral()
    return int_wo_res(t2) - int_wo_res(t1) + res*(Log(t2) - Log(t1))

def Int13(pt,z1,z2): # integral of omega = (y+b)/(x-a)*dx/2y from z1 to z2
    a = pt[0]
    x1 = z1[0]
    x2 = z2[0]
    first_term = Log(x2-a) - Log(x1-a)
    second_term = Int12(pt,z1,z2)
    return 1/2*(first_term + second_term)

##########################################################################################

# C) GENERAL INTEGRATION (endpoints might lie in different components)

CM = C.change_ring(M)
Pref1 = CM(b+3, bx(b+3).sqrt()) # reference point in e_1
Pref2 = CM(b+3,-bx(b+3).sqrt()) # reference point in e_2

def Int1(i,z1,z2): # integral of omega_i from z1 to z2
    if z1 in points0 and z2 in points0:
        return Int01(i,z1,z2)
    elif z1 in points1 and z2 in points1:
        return Int11(i,z1,z2)
    elif z1 in points0 and z2 in points1:
        path_integral   = Int01(i,z1,Pref2) + Int11(i,Pref2,z2)
        period_integral = Int01(i,Pref1,Pref2) + Int11(i,Pref2,Pref1)
        return path_integral - 1/2*period_integral
    elif z1 in points1 and z2 in points0:
        path_integral   = Int11(i,z1,Pref1) + Int01(i,Pref1,z2)
        period_integral = Int01(i,Pref1,Pref2) + Int11(i,Pref2,Pref1)
        return path_integral - 1/2*period_integral
    else:
        return 'NotImplemented'

def Int2(pt,z1,z2): # integral of omega = (y+b)/(x-a)*dx/2y from z1 to z2, where pt=(a,b) is IN v0
    if z1 in points0 and z2 in points0:
        return Int02(pt,z1,z2)
    elif z1 in points1 and z2 in points1:
        return Int13(pt,z1,z2)
    elif z1 in points0 and z2 in points1:
        path_integral   = Int02(pt,z1,Pref2) + Int13(pt,Pref2,z2)
        period_integral = Int02(pt,Pref1,Pref2) + Int13(pt,Pref2,Pref1)
        return path_integral - 1/2*period_integral
    elif z1 in points1 and z2 in points0:
        path_integral   = Int13(pt,z1,Pref1) + Int02(pt,Pref1,z2)
        period_integral = Int02(pt,Pref1,Pref2) + Int13(pt,Pref2,Pref1)
        return path_integral - 1/2*period_integral
    else:
        return 'NotImplemented'

#########################################################################################

# D) LOCAL HEIGHT PAIRING AT P

# Blakestad's constants
alpha_Bl = 1 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 5^5 + 4*5^6 + O(5^7)
beta_Bl  = 2 + 4*5^2 + 3*5^3 + 5^4 + 2*5^5 + O(5^7)
gamma_Bl = 2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 3*5^6 + O(5^7)
delta_Bl = 2 + 4*5^4 + 4*5^5 + 4*5^6 + O(5^7)

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
# k2*(l0*N[0][2] + l1*N[1][2] + l2*N[2][2] + l3*N[3][2]) = O(5^7)

# z1 = (a,b), z2 = (c,d), D = (z1)-(z2) ==> omega = ((y+b)/(x-a)-(y+d)/(x-c))*dx/2y
# is a form of the third kind such that Res(omega) = D.
# Ψ(omega) = d0*omega0 + d1*omega1 + e0*eta0 + e1*eta1 for some d0,d1,e0,d1.
# omegaD := omega - d0*omega0 - d1*omega1 ==> omegaD is the unique form
# of the third kind such that Res(omegaD) = D and Ψ(omegaD) is in W.

def omegaD_coeff(z1,z2):
    global_symbols = Matrix([[Int1(0,z1,z2)],[Int1(1,z1,z2)],[Int1(2,z1,z2)],[Int1(3,z1,z2)]])
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
    return Int2(z1,z4,z3) - Int2(z2,z4,z3) - d0*Int1(0,z4,z3) - d1*Int1(1,z4,z3)

##########################################################################################

# E) EXAMPLE

P, wP = pt1, pt2
Q, wQ = pt3, pt4
R     = pt6
S     = pt9
T, wT = pt11, pt12

# checking quadraticity of the GLOBAL height:
# D1 = (Q)-(wQ), D2 = (R) - (P), D3 = (S) - (wP), D4 = (T) - (wT)
# D1 = 4D2 & D4 = 6D3 ===> 6h(D1,D3) = h(D1,D4) = h(D4,D1) = 4h(D4,D2)
# 1) D1 = (Q)-(wQ), D3 = (S)-(wP)
local_ht_at_p = local_height_at_p(Q,wQ,S,wP)
local_hts_away_from_p = 0 # computed using https://github.com/emresertoz/neron-tate
global_height_1 = local_ht_at_p + local_hts_away_from_p
# 2) D4 = (T)-(wT), D2 = (R)-(P)
local_ht_at_p = local_height_at_p(T,wT,R,P)
local_hts_away_from_p = -2*Log(2) + Log(3) # computed using https://github.com/emresertoz/neron-tate
global_height_2 = local_ht_at_p + local_hts_away_from_p
# check if zero:
6*global_height_1 - 4*global_height_2

##########################################################################################

