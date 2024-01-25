load("../padic_heights_g2/heights_and_abelian_integrals.sage")
R.<x> = PolynomialRing(Rationals())
f = x^5 - 4*x^4 - 48*x^3 + 64*x^2 + 512*x + 256
[b5, b4, b3, b2, b1, b0] = f.coefficients(sparse = False)
p = 7
n = 6
# computed using canonical_subspace
alpha =  3*7 + 2*7^2 + 2*7^4 + 3*7^5 + O(7^6)
beta = 5 + 5*7^2 + 2*7^3 + 3*7^4 + 5*7^5 + O(7^6)
gamma = 7 + 3*7^2 + 2*7^3 + 5*7^4 + 5*7^5 + O(7^6)
delta = 7 + 7^2 + 3*7^3 + 2*7^5 + O(7^6)

c11 = 2*b1*b2 - b1*alpha + b1^2*beta + 3*delta - 3*b1*gamma + 3*b3 
c12 = b2 - b1*beta + 3*gamma
c21 = b2 + alpha - b1*beta #there was a typo in Cliff's thesis for c21
assert c12 == c21, "Not isotropic? There must be an error."
c22 = beta

print("start sigma computation")
sigma_can = sigma_function(f, n + 3, cij = [c11,c12,c22]) #+3 here because leading term of sigma is T1, then terms of degree 3, so the absolute precision can go up to 3*min(ord(T1), ord(T2)) + n
print("finished sigma computation")

# computing m such that mD=0 for all D in J(Fp)
fac = f.change_ring(GF(p)).factor()
l1 = fac[1][0]
m1 = fac[2][0]
h2 = fac[0][0]
print(fac)
#So in terms of Figure 1 of Bruin--Stoll we should be in case 2, second last column.
#so need to multiply by (p-1)#E(Fp). What is E? E is y^2 = m1h2.
pol = m1*h2
E = EllipticCurve([0,pol[2],0,pol[1],pol[0]])
m = (p-1)*len(E.points())
print(m)

H = HyperellipticCurve(f)
J = Jacobian(H)
P = H(-4, 16)
wP = H(-4, -16)
Q = H(-3, -5)
wQ = H(-3, 5)
R = H(0, 16)
wR = H(0, -16)
S = H(4, 16)
T = H(8, 16)
wT = H(8, -16)
U = H(36, 7184)
wU = H(36, -7184)

# D1 = (P ) − (w(P )),D2 = (R) − (w(R)), D3 = (Q) − (S)
# D4 = (T ) − (U ),D5 = (w(T )) − (w(U )).
# We compute h7(D1,D2), h7(D3,D4) and h7(D3,D5)  using Corollary 5.32 of [Bia23]


u1 = J(P)+J(Q)
u2 = J(wP)+J(Q)

htPR = height_at_p(J(R)-u1, sigma_can, m = m)
htPwR = height_at_p(J(wR)-u1, sigma_can, m = m)
htwPR = height_at_p(J(R)-u2, sigma_can, m = m)
htwPwR = height_at_p(J(wR)-u2, sigma_can, m = m)
htD1D2 = -1/2*(htPR-htPwR-htwPR+htwPwR)

print("h7(D1,D2)", htD1D2)


u1 = J(Q)+J(R)
u2 = J(S)+J(R)

htQT = height_at_p(J(T)-u1, sigma_can, m = m)
htST = height_at_p(J(T)-u2, sigma_can, m = m)
htQU = height_at_p(J(U)-u1, sigma_can, m = m)
htSU = height_at_p(J(U)-u2, sigma_can, m = m)
htD3D4 = -1/2*(htQT-htST-htQU+htSU)

print("h7(D3,D4)", htD3D4)

htQwT = height_at_p(J(wT)-u1, sigma_can, m = m)
htSwT = height_at_p(J(wT)-u2, sigma_can, m = m)
htQwU = height_at_p(J(wU)-u1, sigma_can, m = m)
htSwU = height_at_p(J(wU)-u2, sigma_can, m = m)
htD3D5 = -1/2*(htQwT-htSwT-htQwU+htSwU)

print("h7(D3,D5)", htD3D5)
