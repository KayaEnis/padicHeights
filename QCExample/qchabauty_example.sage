# This is code for Example 5.5. We use quadratic Chabauty as in §4.1.1 to
# compute the integral points on y^2 = x^5+x^3-2x+1
# Tested on Sage 9.3 and Magma V2.27-5
# To run the computation, load this file into SageMath.

from sage.rings.padics.precision_error import PrecisionError
load("../QCExample/kummer_formulas.sage")
load("../padic_heights_g2/heights_and_abelian_integrals.sage")
import logging
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    tqdm = lambda x,**kwargs : x

t0 = cputime()

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)
logger = logging.getLogger(__name__)



R.<x> = PolynomialRing(Rationals())
f = x^5 + x^3 - 2*x + 1
[b5, b4, b3, b2, b1, b0] = f.coefficients(sparse = False)
prec = 9 
disc_prec = 15 
p = 5
logger.info('Computing sigma')

if not "naive_sigma" in locals():
    naive_sigma = sigma_function(f, prec + 5) #TO DO: precision
logger.info('Computed sigma')
t1 = cputime()
print("time to compute sigma", t1-t0)
# Generators of J(Q)
H = HyperellipticCurve(f)
J = Jacobian(H)
gen_1 = J(H(0,1)) + J(H(1,-1))
gen_2 = J(H(1,-1))
# J(Q) = <gen_1,gen_2>, See https://www.lmfdb.org/Genus2Curve/Q/85280/d/682240/1

# From gen_1, gen_2, construct generators
# gen_1_prime, gen_2_prime of
# a finite index subgroup of J(Q)
# so that gen_1_prime, gen_2_prime, gen_1_prime + gen_2_prime are not on Theta:
splitting_indices = [[1,1], [1,-2]]
gen_1_prime = splitting_indices[0][0]*gen_1 + splitting_indices[0][1]*gen_2
gen_2_prime = splitting_indices[1][0]*gen_1 + splitting_indices[1][1]*gen_2
index_matrix = Matrix(splitting_indices).transpose()
gens_formal_gp_dict = {"gen_1_prime" : {"coef": 36}, "gen_2_prime" : {"coef": 18}, "gen_1_prime+gen_2_prime": {"coef": 36}}




def compute_neron_at_p_gens():
    HQp = H.change_ring(Qp(p,prec+1))
    JQp = Jacobian(HQp)
    lambda_p_gens_dict = gens_formal_gp_dict.copy()

    for pt, info in lambda_p_gens_dict.items():
        info["value"] = height_at_p(eval(pt), naive_sigma.change_ring(Qp(p,prec+1)), m=info["coef"])
    return lambda_p_gens_dict


def compute_neron_at_2_gens():
    q = 2
    HQq = H.change_ring(Qp(q,10))
    JQq = Jacobian(HQq)
    _.<x> = PolynomialRing(Qp(q,10))
    fq = f.change_ring(GF(q))
    lambda_2_gens_dict = {pt: {"coef": 4} for pt in ["gen_1_prime", "gen_2_prime", "gen_1_prime+gen_2_prime"]}
    for pt, info in lambda_2_gens_dict.items():
        for _ in range(1, 5):
            try:
                info["value"] = height_away_p_rat(eval(pt), q, m = _*info["coef"])
                info["coef"] = _ * info['coef']
                break
            except AssertionError:
                continue
    return lambda_2_gens_dict


def compute_global_heights_gens():
    lambda_p_gens_dict = compute_neron_at_p_gens()

    # Néron functions away from p:
    # Images on the Kummer surface:
    # K(gen_1_prime) = (1 : -1 : 0 : 2), K(gen_2_prime) = (1 : 1 : 0 : -2),
    # K(gen_1_prime + gen_2_prime) = (1 : 0 : -1 : 1)
    # So possible non trivial contributions only at the bad primes (Bia23, Rmk 4.6(ii))
    # Bad primes (besides 5) are: 2, 13, 41
    # Néron fcts at 13 and 41 are trivial.

    lambda_2_gens_dict = compute_neron_at_2_gens()

    global_heights_gens = {}
    for pt in lambda_p_gens_dict.keys():
        global_heights_gens[pt] = lambda_p_gens_dict[pt]["value"] + log(Qp(p,prec)(2)) * lambda_2_gens_dict[pt]["value"]
    return global_heights_gens


def compute_abelian_logs_gens():
    abelian_logs_gens = gens_formal_gp_dict.copy()
    des_prec = prec + max([info["coef"].valuation(p) for info in abelian_logs_gens.values()])
    adj_prec = adjusted_prec_Log(des_prec, p)
    abelian_logs_formal_series,_,_ = strict_log(f, adj_prec)
    for pt, info in abelian_logs_gens.items():
        info["value"] = abelian_logs(eval(pt), p, prec, m = info["coef"], LOG = abelian_logs_formal_series)

    assert [abelian_logs_gens["gen_1_prime"]["value"][i] +  abelian_logs_gens["gen_2_prime"]["value"][i] ==  abelian_logs_gens["gen_1_prime+gen_2_prime"]["value"][i] for i in [0,1]], "log not linear!!!!! ???"
    return abelian_logs_gens, abelian_logs_formal_series


def compute_qc_constants(global_heights_gens, abelian_logs_gens):
    _squared_logs = []
    _global_heights = []
    for pt in global_heights_gens.keys():
        value = abelian_logs_gens[pt]["value"]
        _squared_logs.append([value[i]*value[j] for i in [0,1] for j in [0,1] if  i<=j])
        _global_heights.append(global_heights_gens[pt])

    return Matrix(_squared_logs)**(-1) * vector(_global_heights)


def compute_T1_T2(P, fs, p, n):
    g = 2
    b5, b4, b3, b2, b1, b0 = fs
    pgis = [-P[0][g - i - 1] for i in range(g)]
    pggis = [2*P[1][g - i - 1 ] for i in range(g)]
    #now add p11, p111, p112
    pgis.append( -pgis[0]^3 - pgis[0]*pgis[1] - fs[4]*pgis[0]^2 - fs[3]*pgis[0] - fs[2] + pggis[0]^2 / 4) #p11, using f_6 in Grant
    pggis.append(pgis[1]*pggis[0] - pgis[0]*pggis[1]) #p112, using f_3 in Grant
    pggis.append(2*(pgis[0] + fs[4])*pggis[2] - (pgis[1] + fs[3])*pggis[1] - pgis[2]*pggis[0]) #p111, using f_4 in Grant
    #The orders are [p22, p12, p11] and [p222, p122, p112, p111].
    T1 = - 2*(pgis[2].change_ring(Qp(p, n))/pggis[3].change_ring(Qp(p, n)))
    pofP = pgis[2]*pgis[0] - pgis[1]^2
    XofP = 1/2*(pofP + b2*pgis[1] - b4)
    T2 = - XofP.change_ring(Qp(p, n))/(1/2*pggis[3].change_ring(Qp(p, n)))
    return T1, T2


def compute_formal_parameters_disc_multiple(parametrised_x, parametrised_y, coef, kummer_deltas, kummer_bbmatrix):
    parametrised_disc = [PowerSeriesRing(Rationals(), "t")(_) for _ in [parametrised_x, parametrised_y]]
    kummer_twice_disc = kummer_point(parametrised_disc, parametrised_disc, f)
    kummer_multiple = scalar_multiple(kummer_twice_disc, coef, kummer_deltas, kummer_bbmatrix)
    kummer_multiple = [_.power_series() for _ in kummer_multiple]
    jacobian_multiple = Points(f.change_ring(parametrised_x.parent()), kummer_multiple)[0]
    return compute_T1_T2(jacobian_multiple, [b5, b4, b3, b2, b1, b0] , p, disc_prec)


def compute_disc_expansions(disc_centres, coefs=[], alphas=[], abelian_logs_formal_series =[]):
    L.<t> = LaurentSeriesRing (Qp(p , prec+5))
    HL = H.change_ring(L)
    JL = Jacobian(HL)
    _, kummer_deltas, kummer_bbmatrix = kummer_data(f)
    disc_power_series = []
    for i, disc_centre in tqdm(enumerate(disc_centres)):
        coef = coefs[i]
        parametrised_x, parametrised_y = H.local_coord(disc_centre, prec = disc_prec)

        # Division polynomial
        jacobian_twice_disc = JL(HL(parametrised_x, parametrised_y)) - JL(HL(parametrised_x, -parametrised_y))
        division_value_coef = division_value_new(jacobian_twice_disc, [b5, b4, b3, b2, b1, b0] , coef)

        # Sigma function
        T1_disc_multiple, T2_disc_multiple = compute_formal_parameters_disc_multiple(parametrised_x, parametrised_y, coef, kummer_deltas, kummer_bbmatrix)
        sigma_value = naive_sigma(T1_disc_multiple, T2_disc_multiple)

        sigma_over_division = (sigma_value/division_value_coef).power_series()
        correction_term = (4*parametrised_y^2).change_ring(Qp(p,prec))
        nu_term = -2/coef^2*(log((sigma_over_division/sigma_over_division[0])) + log(sigma_over_division[0])) + log(correction_term/correction_term[0]) + log(correction_term[0])

        # Abelian integrals
        abelian_logs_disc = [_(T1_disc_multiple, T2_disc_multiple)/coef for _ in abelian_logs_formal_series]
        squared_logs_disc = [abelian_logs_disc[i]*abelian_logs_disc[j] for i in [0,1] for j in [0,1] if  i<=j]

        # QC function
        qc_function = nu_term - sum([alphas[i]*squared_logs_disc[i] for i in range(len(alphas))])
        disc_power_series.append(qc_function(p*t))
    return disc_power_series

def compute_multiple(P, m):
    # Compute m*D and (-m)*D for m an integer and D = [P-infty] 
    # Find multiple on the Kummer, then lift
    _, kummer_deltas, kummer_bbmatrix = kummer_data(f)
    kpt = kummer_point([P[0], P[1]], [0,1,0], f)
    kummer_multiple = scalar_multiple(kpt, m, kummer_deltas, kummer_bbmatrix)
    return Points(f.change_ring(Qp(p, prec)), kummer_multiple)


def compute_abelian_logs(D, p, prec, m, LOG):
    fs = f.coefficients(sparse = False)
    adj_prec = adjusted_prec_Log(prec, p)
    T1_D, T2_D = T1_T2(D, fs, p, prec)
    #prec + n- ord(n) >= prec iff n - ord(n) >= 0. Always.
    assert T1_D.valuation() > 0, "Not in the formal group"
    assert T2_D.valuation() > 0, "Not in the formal group"
    assert adj_prec <= min(LOG[0].precision_absolute(), LOG[1].precision_absolute()) , "Input power series precision is too low."

    LOG1D = (LOG[0].truncate(adj_prec)(T1_D, T2_D).constant_coefficient())/m + O(p^prec)
    LOG2D = (LOG[1].truncate(adj_prec)(T1_D, T2_D).constant_coefficient())/m + O(p^prec)
    return LOG1D, LOG2D



def compute_roots(HK, disc_power_series, disc_centres, Gamma): 

    # roots of disc_power_series[i] - Gamma[j] for i,j in [0,1,2]
    # adapted from https://github.com/bianchifrancesca/QC_bielliptic
    points = [[] for i in range(len(Gamma))]
    K = Qp(p, prec)
    for i in range(len(disc_power_series)):
        Q = disc_centres[i]
        rho = disc_power_series[i]
        xx, yy = HK.local_coord(Q, prec = disc_prec)
        t = rho.parent().gens()[0]
        M = rho.prec()
        for l in range(len(Gamma)):
            omega = Gamma[l]
            rhoomega = rho - omega
            k = min(rhoomega[i].valuation(p) for i in range(M))
            N = prec + k
            rhoomega = (p^(-k)*rhoomega).truncate(M) 
            rhoomega_new = 0
            deg = rhoomega.degree()
            for j in range(M):
                if rhoomega[j].valuation() >= N-k:
                    continue
                rhoomega_new += (rhoomega[j] + O(p^(N-k)))*t^j
            rhoomega = rhoomega_new
            NNk = min([rhoomega[i].precision_absolute() for i in range(rhoomega.degree())])
            if NNk < N-k:
                N = NNk + k

            rhoomega_val = rhoomega.valuation()
            assert rhoomega_val <= 1, "multiple root -- increase precision!"
            if rhoomega_val > 0:
                  val_rho = N - k
                  val_der = rhoomega.derivative()(O(p^(N - k))).valuation()
                  assert val_der < (val_rho/2).ceil(), "Hensel not applicable."
                  points[l].append(HK(xx(O(p^(val_rho - val_der + 1))), yy(O(p^(val_rho - val_der + 1)))))


            roots = list(gp.polrootspadic(rhoomega/(rhoomega.parent().0^rhoomega_val), p, 1))
            roots = [p*r + O(p^(N-k+1)) for r in roots if r.valuation(p) >= 0]

            for r in roots:
                val_rho = min([rhoomega(K(sage_eval('%s'%r))/p).valuation(), N-k])
                val_der = (rhoomega.derivative()(K(sage_eval('%s'%r))/p)).valuation()
                assert val_der < (val_rho/2).ceil(), "Hensel not applicable."
                roots[roots.index(r)] = r + O(p^(val_rho - val_der + 1))

            new_points = [HK(xx(K(sage_eval('%s'%t0))), yy(K(sage_eval('%s'%t0)))) for t0 in roots]
            points[l].extend(new_points)

    rational_points = []
    rational_non_integral_points = []
    integral_points = []
    extra_points = []
    for l in range(len(Gamma)):
        for P in points[l]:
            if P == H(0, 1, 0) or P == HK(0, 1, 0):
                continue
            try:
                RP = H.lift_x(QQ(P[0]))
                if RP[1] - P[1] == 0:
                    rational_points.append(RP)
                elif RP[1] + P[1] == 0:
                    RP = H(RP[0], -RP[1])
                    rational_points.append(RP)
                else:
                    extra_points.append(P)
            except ValueError:
                pol = algdep(P[0], 1)
                pol = PolynomialRing(QQ, "x")(pol)
                try:
                    RP = H.lift_x(pol.roots()[0][0])
                    if RP[1] - P[1] == 0:
                        rational_points.append(RP)
                    elif RP[1] + P[1] == 0:
                        RP = H(RP[0], -RP[1])
                        rational_points.append(RP)
                    else:
                        extra_points.append(P)
                except ValueError:
                        extra_points.append(P)
    
    for P in rational_points:
        if P[0] in ZZ:
            integral_points.append(P)
        else:
            rational_non_integral_points.append(P)
    return integral_points, rational_non_integral_points, extra_points

        #assert(len(list(Set(integral_points[l])))) == len(integral_points[l])
        #assert(len(list(Set(extra_points[l])))) == len(extra_points[l])


def main():
    global_heights_gens = compute_global_heights_gens()
    abelian_logs_gens, abelian_logs_formal_series = compute_abelian_logs_gens()
    alphas = compute_qc_constants(global_heights_gens, abelian_logs_gens)
    # Discs to consider:
    # - the singular point (3,0) doesn't lift to any Qp-point
    # - left with (0, 1), (1, 1), (4, 1) up to hyperelliptic involution
    # - they each contain at least one integral point, which we'll take as disc centre
    disc_centres = [H(0,1), H(1,1), H(-1,-1)]
    logger.info('Now power series computations.')
    disc_power_series = compute_disc_expansions(disc_centres, coefs=[36, 36, 36],alphas=alphas,abelian_logs_formal_series = abelian_logs_formal_series)
    K = Qp(p, prec)
    HK = H.change_ring(K)
    # The following points cover all components of the minimal regular
    # model at 2, except the component of infty, where nu_2 is 0
    # See qc_example.m
    gamma1 = height_away_p_rat(2*J(H(1,1)), 2, m=36)-2
    gamma2 = height_away_p_rat(2*J(H(-1,1)), 2, m=36)-2
    gamma3 = height_away_p_rat(2*J(H(4/9, 339/9^3)), 2, m=36)-2
    assert gamma1 == gamma2
    Gamma = [0, -gamma1*log(K(2)), -gamma3*log(K(2))] 
    assert Gamma == [0, 8/3*log(K(2)), 2*log(K(2))]
    
    logger.info('Now compute roots.')
    integral_points, rational_non_integral_points, extra_points = compute_roots(HK, disc_power_series, disc_centres, Gamma)

    logger.info('Found roots. Compute Mordell-Weil sieve coefficients')
    m = 36
#    adj_prec = adjusted_prec_Log(N + k + m.valuation(p),p)
#    abelian_logs_formal_series,_,_ = strict_log(f, adj_prec) 
    logs_gens = [abelian_logs_gens[a]['value'] for a in ['gen_1_prime', 'gen_2_prime']]
    adj = max([a.valuation(p) for a in logs_gens[0] + logs_gens[1] if a!=0])
    M = (Matrix(Qp(p,prec - adj), logs_gens)^(-1)).transpose()

    extra_coeffs = []
    integral_coeffs = []
    assert len(rational_non_integral_points) == 0

    for pt in integral_points:
        multiples = compute_multiple(pt, 36)
        mD = multiples[0] 
        # mD is either m*[pt-infty] or -m*[pt-infty] = iota(pt) -infty 
        # where iota is the hyperell involution
        Logpt = compute_abelian_logs(mD, p, prec, m, abelian_logs_formal_series)
        coeffsP = index_matrix*M*Matrix(2,1,Logpt)
        pol1 = algdep(coeffsP[0][0],1)
        pol2 = algdep(coeffsP[1][0],1)
        integral_coeffs.append([-pol1[0]/pol1[1], -pol2[0]/pol2[1]])

    for pt in extra_points:
        multiples = compute_multiple(pt, 36) # Uses the Kummer to avoid
        mD = multiples[0]                    # precision loss in Cantor
        # mD is either m*[pt-infty] or -m*[pt-infty] = iota(pt) -infty 
        # Since we only need to run the MW sieve for extra points modulo
        # iota, this suffices.
        #Logpt = abelian_logs(J(pt)-J(H(1,1)), p, prec, 36 , LOG = abelian_logs_formal_series)
        Logpt = compute_abelian_logs(mD, p, prec, m, abelian_logs_formal_series)
        coeffsP = index_matrix*M*Matrix(2,1,Logpt)
        extra_coeffs.append([ZZ(coeffsP[0][0]), ZZ(coeffsP[1][0])])

    return alphas, integral_points, extra_points, integral_coeffs, extra_coeffs

alphas, integral_points, extra_points, integral_coeffs, extra_coeffs = main()

# write coefficients of QC roots in terms of bas as magma input
file = open('../QCExample/coeffs_qc_mws.m','w')
file.write("extra_coeffs := ")
file.write(str(extra_coeffs))
file.write(";\n")
file.write("integral_coeffs := ")
file.write(str(integral_coeffs))
file.write(";")
file.close()
   

t2 = cputime()
print("time to run QC and find MWS-coeffs", t2-t1)

# run the Mordell-Weil sieve on `extra` roots
load("../QCExample/qc_example.m")
