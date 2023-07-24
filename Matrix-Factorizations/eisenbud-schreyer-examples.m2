
needs "ZZdFactorizations.m2"

quadraticMF = method()

--input: integer n, and polynomial ring in 2*(n+1) variables, e.g. x_0..y_n
--output: MF for q = \sum_{i=0}^n x_i*y_i
--note: this function really assumes the ring is presented as k[x_0..y_n]

quadraticMF(ZZ, PolynomialRing) := ZZdFactorization => (n, S) -> (
    if dim S != 2*(n+1) then error "expected polynomial ring in 2*(n+1) variables";
    phi0 := matrix{{S_0}};
    psi0 := matrix{{S_(n+1)}};
    phi := phi0;
    psi := psi0;
    for i from 1 to n do(
	phi = matrix{
	    {S_i * id_(S^(2^(i-1))), phi0},
	    {psi0, (-S_(n+1+i))*id_(S^(2^(i-1)))}};
	psi = matrix{
	    {S_(n+1+i) * id_(S^(2^(i-1))), phi0},
	    {psi0, (-S_i)*id_(S^(2^(i-1)))}};
	phi0 = phi;
	psi0 = psi;
	);
    ZZdfactorization{phi0, psi0}	
    )

--input: list {f_0..g_n} of 2(n+1) elements of the same ring (generalizing x_0..y_n)
--output: MF for \sum_{i=0}^n f_i*g_i
quadraticMF(List) := ZZdFactorization => L -> (
    S := ring L#0;
    for i to #L-1 do if not instance(L#i, S) then error "expected list of elements from same ring";
    if #L%2 !=0 then error "expected an even length list";
    n := (#L-2)//2;
    A' := matrix{{L#0}};
    B' := matrix{{L#(n+1)}};
    A := A';
    B := B';
    for i from 1 to n do(
	A = matrix{
	    {L#i * id_(S^(2^(i-1))), A'},
	    {B', (-L#(n+1+i))*id_(S^(2^(i-1)))}};
	B = matrix{
	    {L#(n+1+i) * id_(S^(2^(i-1))), A'},
	    {B', (-L#i)*id_(S^(2^(i-1)))}};
	A' = A;
	B' = B;
	);
    ZZdfactorization{A, B}	
    )

--In: skew symmetric 2(n+1) x 2(n+1) matrix with scalar entries
--In: polynomial ring in 2(n+1) x 2(n+1) variables
--Out: matrix factorization of combination of variables determined by input matrix
quadraticMF(Matrix, PolynomialRing) := ZZdFactorization => (M,S) -> (
    if (rank source M)%2 != 0 then error "expected skew symmetric matrix of even dimension";
    if transpose M + M != 0 then error "expected skew symmetric matrix";
    n := (rank source M - 2)//2;
    if (dim S)!= rank source M then error "expected polynomial ring to have same number of variables as matrix dimension";
    R := ring M;
    G := matrix{{0, id_(R^(n+1))},{id_(R^(n+1)),0}};
    d := (vars S)*(G*M);
    quadraticMF(flatten entries d)
    )

--In: List of distinct scalars d_0..d_n + coefficient ring; the function will make a polynomial ring
--Out: MF associated to the skew symmetric matrix 0 D \\ -D 0
--where D is the diagonal matrix made from the list
quadraticMF(List, Ring) := ZZdFactorization => (L,R) -> (
    if L != unique L then print "Warning: input with repeated entries may not give Ulrich module";
    D := diagonalMatrix L;
    M := matrix{{0, D}, {-D, 0}};
    n := #L-1;
    quadraticMF(M, R[x_0..y_n])
    )

---------------

end--

restart
needs "eisenbud-schreyer-examples.m2"

S = QQ[x_0..y_2]

--this will give factorization of \sum x_i*y_i
Q1 = quadraticMF(2, S)
(Q1.dd_0)*(Q1.dd_1)

--this will give factorization of \sum L_i*L_(n+1+i)
--if the list is length 2(n+1)
Q2 = quadraticMF({x_0, 2*x_1-3*x_2, -x_2, 3*y_0, 2*y_1, 3*y_2})
(Q2.dd_0)*(Q2.dd_1)

--this is closer to the full construction in Eisenbud--Schreyer
S' = QQ[x_0..y_1]
Lambda = matrix{{0,1,2,3}, {-1,0,-1,-2}, {-2,1,0,-3}, {-3,2,3,0}}
Q1 = quadraticMF(1, S')
q1 = ((Q1.dd_0)*(Q1.dd_1))_(0,0)
Q2 = quadraticMF(Lambda, S')
q2 = (Q2.dd_0)*(Q2.dd_1)

A = Q1.dd_1 | Q2.dd_1
R = S'/ideal(q1, q2)
A' = R**A
M = coker A' --this should be an Ulrich module

--diagonal entries
F = quadraticMF({2,-1,4,3}, QQ)
S = ring F
Q1 = quadraticMF(3, S)
q2 = ((F.dd_0)*(F.dd_1))_(0,0)
q1 = ((Q1.dd_0)*(Q1.dd_1))_(0,0)
R = S/ideal(q1, q2)
M = coker(sub(Q1.dd_1 | F.dd_1, R))


quadraticMF(ZZ, PolynomialRing) := ZZdFactorization => (n, S) -> (
    if dim S != 2*(n+1) then error "expected polynomial ring in 2*(n+1) variables";
    phi0 := matrix{{S_0}};
    psi0 := matrix{{S_(n+1)}};
    phi := phi0;
    psi := psi0;
    for i from 1 to n do(
	phi = matrix{
	    {S_i * id_(S^(2^(i-1))), phi0},
	    {psi0, (-S_(n+1+i))*id_(S^(2^(i-1)))}};
	psi = matrix{
	    {S_(n+1+i) * id_(S^(2^(i-1))), phi0},
	    {psi0, (-S_i)*id_(S^(2^(i-1)))}};
	phi0 = phi;
	psi0 = psi;
	);
    ZZdfactorization{phi0, psi0}	
    )

S = QQ[x_0..y_2]
phi0 = matrix{{S_0}}
psi0 = matrix{{S_(2+1)}}
n=2

--this will make everything homogeneous
for i from 1 to n do(
    phi = matrix{
    	{S_i * matrix( map(S^(2^(i-1)), S^{2^(i-1):-1}, id_(S^(2^(i-1))))), phi0},
    	{psi0, -S_(n+i+1) * matrix( map(S^(2^(i-1)), S^{2^(i-1):-1}, id_(S^(2^(i-1)))))}
    	};
    psi = matrix{
    	{S_(n+i+1) * matrix( map(S^(2^(i-1)), S^{2^(i-1):-1}, id_(S^(2^(i-1))))), phi0},
    	{psi0, -S_i * matrix( map(S^(2^(i-1)), S^{2^(i-1):-1}, id_(S^(2^(i-1)))))}
    	};
    phi0 = phi;
    psi0 = psi;
    );

phi0
psi0
degrees source psi0
isHomogeneous psi0
isHomogeneous phi0
