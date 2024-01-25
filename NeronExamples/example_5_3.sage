load("../padic_heights_g2/heights_and_abelian_integrals.sage")
R.<x> = PolynomialRing(Rationals())
f = x^5+5*x^4-168*x^3+1584*x^2-10368*x+20736
[b5, b4, b3, b2, b1, b0] = f.coefficients(sparse = False)
p = 5
n = 7
# computed using canonical_subspace
alpha =  1 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 5^5 + 4*5^6 + O(5^7)
beta = 2 + 4*5^2 + 3*5^3 + 5^4 + 2*5^5 + O(5^7)
gamma = 2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 3*5^6 + O(5^7)
delta = 2 + 4*5^4 + 4*5^5 + 4*5^6 + O(5^7)

c11 = 2*b1*b2 - b1*alpha + b1^2*beta + 3*delta - 3*b1*gamma + 3*b3 
c12 = b2 - b1*beta + 3*gamma
c21 = b2 + alpha - b1*beta 
assert c12 == c21, "Not isotropic? There must be an error."
c22 = beta

sigma_can = sigma_function(f, n + 3, cij = [c11,c12,c22]) #+3 here because leading term of sigma is T1, then terms of degree 3, so the absolute precision can go up to 3*min(ord(T1), ord(T2)) + n

# computing m such that mD=0 for all D in J(Fp)
fac = f.change_ring(GF(p)).factor()
len(fac)
l1 = fac[1][0]
m1 = fac[2][0]
h2 = fac[0][0]
assert not l1.resultant(h2).is_square()
assert not m1.resultant(h2).is_square()
#So according to Bruin-Stoll, the order is (p+1)^2 = 36

H = HyperellipticCurve(f)
J = Jacobian(H)
P = H(-12, 720)
wP = H(-12, -720)
Q = H(-8,528)
wQ = H(-8,-528)
R = H(0,-144)
wR = H(0,144)
S = H(12,432)
T = H(36, 7920)
wT = H(36, -7920)
U = H(8,80)
wU = H(8, -80)

# We compute h5(D1,D3) and h5(D4,D2) using Corollary 5.32 of [Bia23]


# D1 = (Q) − (w(Q)), D2 = (R) − (P ), D3 = (S) − (w(P )), D4 = (T ) − (w(T )).
u1 = J(Q)+J(R)
u2 = J(wQ)+J(R)

htQS = height_at_p(J(S)-u1, sigma_can, m = 36)
htQwP = height_at_p(J(wP)-u1, sigma_can, m = 36)
htwQS = height_at_p(J(S)-u2, sigma_can, m = 36)
htwQwP = height_at_p(J(wP)-u2, sigma_can, m = 36)
htD1D3 = -1/2*(htQS-htQwP-htwQS+htwQwP)

print("h5(D1,D3)", htD1D3)

#
#htQS = height_at_p(J(Q)-J(S), sigma_can, m = 36)
#htQwP = height_at_p(J(Q)-J(wP), sigma_can, m = 36)
#htwQS = height_at_p(J(wQ)-J(S), sigma_can, m = 36)
#htwQwP = height_at_p(J(wQ)-J(wP), sigma_can, m = 36)
#htD1D3 = -1/2*(htQS-htQwP-htwQS+htwQwP)
#
#print("h5(D1,D3)", htD1D3)

u1 = J(T)+J(Q)
u2 = J(wT)+J(Q)

htTR = height_at_p(J(R)-u1, sigma_can, m = 36)
htTP = height_at_p(J(R)-u2, sigma_can, m = 36)
htwTR = height_at_p(J(P)-u1, sigma_can, m = 36)
htwTP = height_at_p(J(P)-u2, sigma_can, m = 36)
htD4D2 = -1/2*(htTR-htTP-htwTR+htwTP)

print("h5(D4,D2)", htD4D2)



