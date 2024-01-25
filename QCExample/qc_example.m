load "../QCExample/mws_qc.m";
load "../QCExample/coeffs_qc_mws.m";

// Magma computations for Example 5.5
// In particular
// - check the claims about the regular model at 2
// - run the Mordell-Weil sieve to show ``extra'' points aren't rational.

f := x^5 + x^3 - 2*x + 1;
C := HyperellipticCurve(f); 
J := Jacobian(C);
p := 5;
N := 10;
gC := Genus(C);

// Check that the regular model at 2computed by magma has 4 irreducible
// components and that the four points given below reduce to distinct
// components.
C2 := RegularModel(C, 2);
assert #Components(C2) eq 4;
assert #[PointOnRegularModel(C2, P)`component_indices : P in [C![1,0,0], 
       C![-1,1], C![1,1], C![4/9, 339/9^3]]] eq 4;

ptsC := Points(C : Bound:=100);
intptsC := [P : P in ptsC | P[3] eq 1];

"Small integral points: ", intptsC;
K := pAdicField(5,8);

bas := [C![0,1]-C![1,1], C![1,-1] - C![1,0,0]];
// J(Q) = <bas>, See https://www.lmfdb.org/Genus2Curve/Q/85280/d/682240/1
assert #TorsionSubgroup(J) eq 1;
base_pt := [1,0,0];

fake_coeffs :=  [[extra_coeffs]] ;

""; "Use the Mordell-Weil sieve to show that the extra points aren't rational.";

qc_fake_coeffs_mod_M := coeffs_CRT(fake_coeffs, [5], [6]);
"number of cosets", #qc_fake_coeffs_mod_M;
qc_M := 5^6; // modulus

mws_primes := [311]; // chosen because J(F_311) = (Z/2Z)^2 + Z/(2^3*5^5)Z
printf "starting MW-sieve modulo %o\n", 311;
time done_fake := MWSieve(J, mws_primes, qc_M, bas, C!base_pt, qc_fake_coeffs_mod_M : known_rat_coeffs := integral_coeffs );
assert done_fake;
printf "No additional solutions are rational\n";
printf "Hence the set of integral points on \n%o is precisely \n %o\n", C, intptsC;
