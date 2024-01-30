##########################################################################################

import time
_.<x> = PolynomialRing(QQ)

def my_local_coords(f, p, n, prec):
    """
    Compute a local coordinate at infinity on the odd degree hyperelliptic 
    curve y^2 = f(x) modulo p^n to prec digits of precision   
    """
    pol = f.change_ring(Integers(p^n))
    g = pol.degree()//2
    
    K = LaurentSeriesRing(Integers(p^n), name='t', default_prec=prec+2)
    t = K.gen()
    L = PolynomialRing(K,'x')
    x = L.gen()
    i = 0 
    w = (x**g/t)**2-pol
    wprime = w.derivative(x)
    x = t**-2
    for i in range((RR(log(prec+2)/log(2))).ceil()):
        x = x - w(x)/wprime(x)
    y = x**g/t
    return x+O(t**(prec+2)) , y+O(t**(prec+2))

def blakestad_zeta(b, p, n):
    r"""
    Compute Blakestad's zeta functions zeta_1 and zeta_2 on the genus 2 curve
    C: y^2 = b(x) modulo p^n. Here b is a monic quintic over the integers and 
    p is a prime of ordinary reduction 
    """
    bx,by = my_local_coords(b, p, n, 3*p^n)
    if by[-5] == 1:
        by = -by
    assert bx.parent().gens()[0] == -bx^2/by
    by = -by
    x_powers = [bx^i for i in range(5)]
    i = ZZ(mod(3*p^n, 5)/2)
    j = (3*p^n-2*i) // 5 # 3p^n = 2i+5j
    yj = by^j
    ysquareinverse = (by^2)^(-1)
    rho3pn = x_powers[i]*yj
    rhopn = by.parent().zero()
    k = p^n
    while j > 0: 
        while i > 0:
            i -= 1
            m = 2*i+5*j
            if mod(m,500) == 0:
                print("m", m)
            if m == k:
                rhopn = x_powers[i]*yj
            if rho3pn[-m] != 0 or rhopn[-m] != 0:
                x_powers[i] = x_powers[i].truncate_laurentseries(4+5*j)
                # Truncate to get intermediate result up to O(t^4)
                prod = x_powers[i]*yj
                rho3pn -= rho3pn[-m]*prod
                if m < k and rhopn[-m] != 0:
                    rhopn -= rhopn[-m]*prod
        i = 5
        # Truncate to get intermediate result up to O(t^4)
        yj *= ysquareinverse.truncate_laurentseries(-yj.valuation()+12) 
        j -= 2  #No even terms

    An = rho3pn[-3]
    Bn = rho3pn[-1]
    In = rho3pn[1]
    Rn = rho3pn[3]
    Cn = rhopn[-3]
    Dn = rhopn[-1]
    Jn = rhopn[1]
    Sn = rhopn[3]
    M = Matrix(2,2,[An,Bn,Cn,Dn])
    N = Matrix(2,4,[[An,Bn,In,Rn],[Cn,Dn,Jn,Sn]])
    assert M.determinant().valuation(p) == 0, "Hasse--Witt matrix not invertible mod p"
    return M^(-1)*vector([rho3pn,rhopn]), rhopn, rho3pn 

def canonical_constants(b, p, n):
    r"""
    Compute the canonical complementary subspace of H^1_dR for
    C: y^2 = b(x) modulo p^n. Here b is a monic quintic over the integers and 
    p is a prime of ordinary reduction.
    """
    H = HyperellipticCurve(b)
    zeta, _, _ = blakestad_zeta(b, p, n) 
    alphan = Qp(p,n)(zeta[0][1])
    deltan = Qp(p,n)(zeta[0][3])
    betan  = Qp(p,n)(zeta[1][1])
    gamman = Qp(p,n)(zeta[1][3])
    assert alphan == 3*gamman
    b0,b1,b2,b3,b4,b5 = b.coefficients()
    # eta1 = [3*x^3 + 3*b1*x^2 - alphan*x - (3*b1*b2-b1*alphan+3*b3+3*deltan)] times dx/2y
    # eta2 = [x^2 - betan*x - (b2-b1*betan+3*gamman)] times dx/2y
    f11 = 2*b1*b2 - b1*alphan + b1^2*betan + 3*deltan - 3*b1*gamman + 3*b3
    f12 = b2 - b1*betan + 3*gamman
    f22 = betan
    return [alphan, betan, gamman, deltan], [f11, f12, f22]

##########################################################################################

t = time.time()
# bx = x^5 + 5*x^4 - 168*x^3 + 1584*x^2 - 10368*x + 20736
# p = 5
# n = 7
bx = x^5 - 4*x^4 - 48*x^3 + 64*x^2 + 512*x + 256 
p = 7
n = 4  # Used higher precision for the example, but that takes very long...
L, _ = canonical_constants(bx, p, n)
t = time.time() - t

print('------------------------------------------------------------------------------------')
print('Blakestad`s constants up to prec ' + str(n))
print('alpha = ' + str(L[0]))
print('beta  = ' + str(L[1]))
print('gamma = ' + str(L[2]))
print('delta = ' + str(L[3]))
print('------------------------------------------------------------------------------------')
print('Time: ' + str(t))
print('------------------------------------------------------------------------------------')

##########################################################################################
