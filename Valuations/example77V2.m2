needsPackage "Valuations"
needsPackage "SubalgebraBases"
needsPackage "Tropical"
needsPackage "gfanInterface"
needsPackage "Binomials"

-- Given an ideal I compute the prime cones of Trop(I)
--
-- A prime cone is one whose corresponding ideal is prime and binomial
-- so we use the Binomials package to check primality (since we are
-- working over \CC)
--
primeConesOfIdeal = I -> (
    F:=tropicalVariety(I, IsHomogeneous=>true,Prime=>true);
    r:=rays(F);
    c:=maxCones(F);
    cns := for i in c list(r_i);
    inCns := for c in cns list (flatten entries( c * transpose matrix{toList(numColumns(c) : 1)}));
    L:= for i from 0 to #cns-1 list (J = gfanBuchberger(I, "w" => -1*(inCns#i));
	H := gfanInitialForms(J, -1*(inCns#i), "ideal" =>true);
	K := H_1;
	if binomialIsPrime(ideal(K)) then cns#i);
    delete(null,L)
    )
						    
-- given a set of rays of a 2D cone,
-- get two interior points of the cone that span it (as a vector space)
coneToMatrix = coneRays -> (
    independentConeRays := getMaxIndependent(coneRays);
    coeffs := matrix for i from 1 to numcols independentConeRays list for j from 1 to numcols independentConeRays list if i == j then 2 else 1;
    print(coeffs);
    print(independentConeRays);
    coeffs*(transpose independentConeRays)
    )

-- get a maximal set of independent columns of a matrix
getMaxIndependent = M -> (
    -- compute the pivot columns to obtain a maximal linearly independent subset of columns of M
    R := reducedRowEchelonForm(sub(M, QQ));
    P := for i from 1 to rank M list min (for j from 1 to numcols M list if R_(i-1,j-1) != 0 then j-1 else numcols M);
    M_P
    )

-- scale the rows of a list of matrices
-- using a positive vector of the lineality space of a tropical
-- variety f
positivity = (f, matL) -> (
    l := transpose linealitySpace(f);
    finalScaledMats := {};
    matList := for i from 0 to #matL-1 list entries matL_i;    
    for i from 0 to #matList-1 do (
	scaledRows = {};
	for j from 0 to #(matList_i)-1 do (
	    coeff := -1*floor(min apply(#(matList_i)_j, k -> (((matList_i)_j)_k)/(flatten entries l)_k));
	    scaledRows = append(scaledRows, (1/gcd(flatten entries (coeff*l + matrix{(matList_i)_j})))*(coeff*l + matrix{(matList_i)_j}));
	    );
	mat := scaledRows_0;
	for i from 1 to #scaledRows-1 do mat = mat || scaledRows_i;
	finalScaledMats = append(finalScaledMats, mat);
    	);
    finalScaledMats
    )

coneToValuation = (coneRays, I) -> (
    F := tropicalVariety(I, IsHomogeneous=>true,Prime=>true);
    M := coneToMatrix(coneRays);
    scaledM := (positivity(F, {M}))/(i -> sub(i, ZZ));
    T := QQ[e_1, e_2, e_3, y, MonomialOrder=>{Weights=>((entries scaledM_0)_0), Weights=>((entries scaledM_0)_1)}];
    val := leadTermValuation(T);
    orderedM := orderedQQn(2, {Lex});
    func := (f -> (
	    valf := val(sub(f, T));
	    if valf == infinity then infinity else (
		(gens orderedM)*(scaledM_0)*(valf)
		)
	    )
	);
    valuation(func, S, orderedM)
    )

-- construct the new valuation by taking min
valM = (T, valMTwiddle) -> (
    valMfunc = (g) -> (
	R := QQ[x_1, x_2, x_3, e_1, e_2, e_3, y, MonomialOrder => Eliminate 3];
	I := ideal{x_1 + x_2 + x_3 - e_1, x_1*x_2 + x_1*x_3 + x_2*x_3 - e_2, x_1*x_2*x_3 - e_3, (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3) - y};
	f := e_1^2*e_2^2 - 4*e_2^3 - 4*e_3*e_1^3 + 18*e_1*e_2*e_3 - 27*e_3^2 - y^2;
	S := valMTwiddle#"domain";
	m := map(S, R, matrix{{0,0,0}} | matrix {gens S});
	gTwiddle := m (sub(g, R) % I);
	maxTwiddle := gTwiddle % ideal(sub(f, S));
	use T; -- something above changes the user's ring (what could it be?) let's assume it was T
	valMTwiddle(maxTwiddle)
	);
    valuation(valMfunc, T, valMTwiddle#"codomain")
    )


end ---

restart
load "example77V2.m2"

-- setup the example
R = QQ[x_1, x_2, x_3];

A = subring {
    x_1 + x_2 + x_3,
    x_1*x_2 + x_1*x_3 + x_2*x_3,
    x_1*x_2*x_3,
    (x_1 - x_2)*(x_1 - x_3)*(x_2 - x_3)
    };

S = QQ[e_1, e_2, e_3, y];

presMap = map(R, S, gens A);
I = ker presMap

-- The primes cones of the tropical variety:
C = primeConesOfIdeal I

-- turn them into weights:
flatten (C/coneToMatrix/(i -> positivity(tropicalVariety I, {i})))

-- create weight valuations on the polynomial ring S
v0 = coneToValuation(C#0, I);
v1 = coneToValuation(C#1, I);
v2 = coneToValuation(C#2, I);
use S;

v0(e_1^2 + e_2*e_3 - y^3) -- lead term from e_1^2
v1(e_1^2 + e_2*e_3 - y^3) -- lead term from y^3
v2(e_1^2 + e_2*e_3 - y^3) -- lead term from e_2*e_3


-- create the induced valuation on the subring A
vA0 = valM(R, v0);
vA1 = valM(R, v1);
vA2 = valM(R, v2);
use R;

vA0(x_1^2 + x_2^2 + x_3^2)
vA1(x_1^2 + x_2^2 + x_3^2)
vA2(x_1^2 + x_2^2 + x_3^2)

vA0((x_1^2 - x_2^2)*(x_1^2 - x_3^2)*(x_2^2 - x_3^2))
vA1((x_1^2 - x_2^2)*(x_1^2 - x_3^2)*(x_2^2 - x_3^2))
vA2((x_1^2 - x_2^2)*(x_1^2 - x_3^2)*(x_2^2 - x_3^2))

vA0(0_R)

-- Note, for elements not in A, the valuation returns nonsense 
-- because the valuation does not come from a weight valuation
-- on R
vA0(x_2)
vA0(x_2^2)
vA0(x_2^3)
