newPackage("Valuations",
        Headline => "implementation of valuations for rings",
        Version => "1.0",
        Date => "June 5, 2023",
        Authors => {
            {Name => "Michael Burr", Email => "burr2@clemson.edu", HomePage => "https://cecas.clemson.edu/~burr2/"},
            {Name => "Colin Alstad", Email => "calstad@clemson.edu"},
            {Name => "Michael Byrd", Email => "mbyrd6@clemson.edu", HomePage => "https://michael-byrd.github.io"},
            {Name => "Ethan Partida", Email => "ethan_partida@brown.edu", HomePage => "https://ethanpartida.github.io/"},
            {Name => "Shelby Cox", Email => "spcox@umich.edu"},
            {Name => "Courtney George", Email => "courtney.george@uky.edu"},
            {Name => "Oliver Clarke", Email => "oliver.clarke@ed.ac.uk", HomePage => "oliverclarkemath.com"}},
        DebuggingMode => false,
        HomePage => "https://github.com/Macaulay2/Workshop-2023-Minneapolis/tree/valuations",
        Configuration => {},
        PackageExports => {"LocalRings", "SubalgebraBases", "InvariantRing", "Tropical", "gfanInterface", "Binomials"}
        )

-- importFrom_"LocalRings" {"LocalRing"}
-- importFrom_LocalRings {"LocalRing"}
-- importFrom_SubalgebraBases {"Subring"}

----- MOVE TO SubalgebraBases -> Subring
ring Subring := A -> ambient A
ring RingOfInvariants := A -> ambient A
ring LocalRing := A -> A
ambient LocalRing := A -> A

export{"valuation",
       "trivialValuation",
       "padicValuation",
       "leadTermValuation",
       "lowestTermValuation",
       "localRingValuation",
       "getMExponent",
       "domain",
       "codomain",
       "orderedQQn",
       "valM"
       }
   
OrderedQQn = new Type of Module
OrderedQQVector = new Type of Vector

--------------------------------------------------------------------------------
-------------------------------- Valuation Type --------------------------------
--------------------------------------------------------------------------------

Valuation = new Type of HashTable

valuation = method()

valuation Function := v -> (
    internalValuation(v, null, null)
    )

ourSources := {Ring,Subring,LocalRing,RingOfInvariants}
ourTargets := {Ring,Subring,LocalRing,RingOfInvariants,OrderedQQn}

for i in ourSources do (
    for j in ourTargets do (
        valuation (Function, i, j) := (v, S, T) -> (
            internalValuation(v, S, T)
            )
        )
    )

-- parameter class was Type, but the null
-- input is not a "Type", but a "Thing"
internalValuation = method()
internalValuation (Function, Thing, Thing) := (v, S, T) -> (
    new Valuation from{
        "function" => v,
        "domain" => S,
        "codomain" => T,
        cache => new CacheTable
        }
    )

-- Concerns with subrings and local rings, will need testing.
Valuation Thing := (v,t) -> (
    if (v#"domain" === null) or (ring t) === v#"domain" then
        v#"function" t
    else if (isMember(ring t, v#"domain".baseRings)) then
        v#"function" promote(t, v#"domain")
    )

--------------------------------------------------------------------------------
--------------------------- Ordered QQ-module Types ----------------------------
--------------------------------------------------------------------------------

-- Ordered Module based on the monomial order of a polyomial ring
--
-- given two elements a, b in QQ^n they are compared by using the
-- the monomial order of the polynomial ring:
-- 1) clear the denominators of a and b: d*a, d*b \in \ZZ^n
-- 2) clear the negative values: d*a + c, d*b + c \in \NN^n
-- 3) compare x^(d*a + c) and x^(d*b + c) in the polynomial ring

orderedQQn = method()
orderedQQn(PolynomialRing) := R -> (
    n := numgens R;
    ordMod := new OrderedQQn of OrderedQQVector from QQ^n;
    ordMod.cache.Ring = R;
    ordMod
    )
orderedQQn(ZZ, List) := (n, monOrder) -> (
    R := QQ[Variables => n, MonomialOrder => monOrder];
    ordMod := orderedQQn R;
    ordMod
    )

-- Two ordered modules are equal iff their cached rings are identitcal
OrderedQQn == OrderedQQn := (N, M) -> (
    N.cache.Ring === M.cache.Ring
    )

OrderedQQn#{Standard,AfterPrint} =
OrderedQQn#{Standard,AfterNoPrint} = M -> (
    << endl; -- double space
    << concatenate(interpreterDepth:"o") << lineNumber;
    << " : Ordered QQ^" | toString (numgens M) | " module" << endl
    );

--
-- comparison of ordered vectors 
OrderedQQVector ? OrderedQQVector := (a, b) -> (
    M := class a;
    N := class b;
    assert(M == N);
    assert(instance(M, OrderedQQn));

    d := lcm((entries a | entries b)/denominator);
    aScaled := d*a;
    bScaled := d*b;

    R := M.cache.Ring;
    c := for i from 0 to numgens R-1 list min(a_i, b_i, 0);

    aMonomial := product for i from 0 to numgens R-1 list (R_i)^(sub(aScaled_i - c_i,ZZ));
    bMonomial := product for i from 0 to numgens R-1 list (R_i)^(sub(bScaled_i - c_i,ZZ));

    if aMonomial > bMonomial then symbol <
    else if aMonomial < bMonomial then symbol >
    else symbol ==
    )

OrderedQQVector == InfiniteNumber := (a, b) -> false
InfiniteNumber ==  OrderedQQn := (a, b) -> false

--
-- monomialToOrderedQQVector
--
-- A function that takes a monomial and an ordered QQ-module and returns the
-- exponent vector of the monomial as a vector in the passed QQ-module    
monomialToOrderedQQVector = method()
monomialToOrderedQQVector (RingElement, OrderedQQn) := (monomial, orderedQQModule) -> (
    exponentVector := vector flatten exponents monomial;
    modGens := gens orderedQQModule;
    modGens*exponentVector
    )
--------------------------------------------------------------------------------
------------------------- Built-in Valuation Functions -------------------------
--------------------------------------------------------------------------------

-- Trivial Valuation
trivialValuation = valuation (x -> if x == 0 then infinity else 0)

-- p-adic Valuation
countPrimeFactor = (p, x) -> (
    -- Returns the number of times that p divides x
    numFactors := 0;
    while x % p == 0 do (
        x = x // p;
        numFactors = numFactors + 1
        );
    numFactors
    )

padicValuation = method()
padicValuation ZZ := p -> (
    if not isPrime p then error "expected a prime integer";
    func := x -> (
        if x == 0 then infinity
        else countPrimeFactor(p, numerator x_QQ) - countPrimeFactor(p, denominator x_QQ)
        );
    valuation func
    )

-- Leading Term Valuation (max convention)
leadTermValuation = method()
leadTermValuation PolynomialRing := R -> (
    monOrder := (options R).MonomialOrder;
    orderedMod := orderedQQn(R);
    valFunc := f -> (print f;print describe f;if f == 0 then infinity else monomialToOrderedQQVector(leadTerm f, orderedMod));
    internalValuation(valFunc, R, orderedMod)
    )

-- Lowest Term Valuation
lowestTermValuation = method()
lowestTermValuation PolynomialRing := R -> (
    monOrder := (options R).MonomialOrder;
    orderedMod := orderedQQn(R);
    valFunc := f -> (
        if f == 0_R then infinity
        else monomialToOrderedQQVector((sort flatten entries monomials f)_0, orderedMod)
        );
    internalValuation(valFunc, R, orderedMod)
    )

-- Local Ring Valuation
getMExponent = method()
getMExponent (Ideal, RingElement) := (m, x) -> (
    numFactors := 0;
    n := m;
    while x % n == 0 do (
        numFactors = numFactors + 1;
        n = n*m;
        );
    numFactors
    )

localRingValuation = method()
localRingValuation LocalRing := R -> (
    m := R.maxIdeal;
    S := ring m;
    func := x -> (
        if x == 0 then infinity
        else getMExponent(m, sub(x, S))
        );
    valuation(func, R, ZZ)
    )

--------------------------------------------------------------------------------
-------------------------------- example77 Valuation ---------------------------
--------------------------------------------------------------------------------

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

-- TODO need to generalize!
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

-- TODO need to generalize!
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

--------------------------------------------------------------------------------
-------------------------------- Documentation ---------------------------------
--------------------------------------------------------------------------------

beginDocumentation()

doc ///
     Key
         "trivialValuation"
     Headline
         Constructs the trivial valuation
     Usage
         v = trivialValuation
     Outputs
         v:Valuation
             the trivial valuation
     Description
       Text
           A function to construct the trivial valuation, returning infinity when the valuation input is zero and returning zero otherwise.
       Example
           v = trivialValuation;
           v (-13)
           v 100000000
           v (14/23)
           v 0
     SeeAlso
         valuation
     ///


doc ///
     Key
         padicValuation
         (padicValuation, ZZ)
     Headline
         Construct a p-adic valuation
     Usage
         v = padicValuation(p)
     Inputs
         p:ZZ
             a prime
     Outputs
         v:Valuation
             p-adic valuation using prime p
     Description
       Text
           A function to construct the p-adic valuation, returning the number of times that p divides the numerator minus the number of times it divides the denominator.
       Example
           v = padicValuation 7;
           v 98
           v (2/7)
           v 0
           v (-42)
     SeeAlso
         valuation
     ///

doc ///
     Key
        lowestTermValuation
        (lowestTermValuation, PolynomialRing)
     Headline
        The valuation which returns the lowest term of an element of an ordered ring
     Usage
         v = lowestTermValuation
     Inputs
         R:PolynomialRing
     Outputs
         v:Valuation
             the lowest term valuation
     Description
       Text
           This valuation returns the lowest (trailing) term of a polynomial with respect to the ring's term order.
	   The valuation returns a value in an ordered module. For more details see @TO "Ordered modules"@.
       Example
           R = QQ[a,b,c, MonomialOrder => Lex];
           vR = lowestTermValuation R;
           f = 13*a^2*b + a*c^3;
           vR f
           S = QQ[a,b,c, MonomialOrder => RevLex, Global => false];
           vS = lowestTermValuation S;
           f = 13*a^2*b + a*c^3;
           vS f
     SeeAlso
         valuation
     ///

doc ///
     Key
         valuation
         (valuation, Function)
	 (valuation, Function, Ring, Ring) 
	 (valuation, Function, Ring, Subring)
	 (valuation, Function, Ring, LocalRing)
	 (valuation, Function, Ring, RingOfInvariants)
	 (valuation, Function, Subring, Ring)
	 (valuation, Function, Subring, Subring)
	 (valuation, Function, Subring, LocalRing)
	 (valuation, Function, Subring, RingOfInvariants)
	 (valuation, Function, LocalRing, Ring)
	 (valuation, Function, LocalRing, Subring)
	 (valuation, Function, LocalRing, LocalRing)
	 (valuation, Function, LocalRing, RingOfInvariants)
	 (valuation, Function, RingOfInvariants, Ring)
	 (valuation, Function, RingOfInvariants, Subring)
	 (valuation, Function, RingOfInvariants, LocalRing)
	 (valuation, Function, RingOfInvariants, RingOfInvariants)
     Headline
         Constructs a user defined valuation object
     Usage
         v = valuation(f)
         v = valuation(f, S, T)
     Inputs
         f:Function
         S:{Ring,LocalRing,Subring}
         T:{Ring,LocalRing,Subring}
     Outputs
         v:Valuation
            user defined valuation function
     Description
         Text
             A function to construct a user defined valuation function.
         Example
             v = valuation(x -> if x == 0 then infinity else 0)
             v = valuation(x -> if x == 0 then infinity else 0, ZZ, ZZ)
     SeeAlso
          lowestTermValuation
	  padicValuation
	  "trivialValuation"
	  leadTermValuation
	  localRingValuation
///

doc ///
     Key
        leadTermValuation
        (leadTermValuation, PolynomialRing)
     Headline
        The valuation which returns the exponent of the lead term of an element of an ordered ring
     Usage
         v = leadTermValuation R
     Inputs
         R:PolynomialRing
     Outputs
         v:Valuation
             the lead term valuation
     Description
       Text
           This valuation returns the exponent vector of the lead term of a polynomial with respect to the ring's term order.
           The valuation returns vectors in an \textit{ordered $\QQ$-module}, which respects the monomial order of the
           given @TO "PolynomialRing"@. For more details see @TO "Ordered modules"@.
       Example
           R = QQ[a,b,c, MonomialOrder => Lex];
           v = leadTermValuation R;
           f = 13*a^2*b + a*c^3;
           g = 5*a^2*c + b^3;
           v f
           v g
           v f > v g
     SeeAlso
         valuation
///

doc ///
     Key
         localRingValuation
         (localRingValuation, LocalRing)
     Headline
         Construct a local ring valuation given a local ring.
     Usage
         v = localRingValuation(R)
     Inputs
         R:LocalRing
             a local ring
     Outputs
         v:Valuation
             local ring valuation using a local ring R
     Description
       Text
           A function to construct a local ring valuation. This function returns the largest power of the maximal ideal
	   of R that contains the input of the valuation.
       Example
           R = QQ[x,y];
	   I = ideal(x,y);
	   S = R_I
	   localVal = localRingValuation(S)
	   localVal(1 + x + y)
	   localVal(x^4 + x^2*y^2 + x^7 + y^3)
	   localVal(x^2 + x*y + y^2)
     SeeAlso
         valuation
     ///
     
doc ///
      Key
      	  Valuations
      Headline
      	  A package for constructing and using valuations.
      Description
        Text
	      A valuation is a function $v:R\rightarrow G\cup\{\infty\}$ 
	      where $R$ is a ring and $G$ is a linearly ordered group with
	      the following properties: $v(ab)=v(a)+v(b)$, 
	      $v(a+b)\geq\min\{v(a),v(b)\}$, and $v(a)=\infty$ iff $a=0$.
          As implemented in @TT "Macaulay2"@, a valuation acts like @ofClass Function@,
          perhaps with extra information.
        Example
          pval = padicValuation 3;
          pval(54)
          R = QQ[x,y];
          leadval = leadTermValuation R;
          leadval(x^3+3*x^3*y^2+2*y^4)
          lowestval = lowestTermValuation R;
          lowestval(x^3+3*x^3*y^2+2*y^4)
          pval(0)
      ///

doc ///
      Key
      	  Valuation
      Headline
      	  The type of all valuations
      Description
      	  Text
	      @TT "Valuation"@ is a type that contains the data needed 
	      to evaluate a valuation.
      SeeAlso
      	  valuation
      ///

doc ///
     Key
    	"Ordered modules"
     Headline
         Overview of the ordered module $\QQ^n$
     Description
       Text
           Many standard valuations take values in a totally ordered subgroup $\Gamma \subseteq \QQ^n$.
	   These standard valuations implement the \textit{OrderedQQn} object, whose order is based on the
	   monomial order of a given ring $R$.
	   The values in $\QQ^n$ are compared using the monomial order of $R$.
       Example
           R = QQ[x,y];
	   I = ideal(x,y);
	   v = leadTermValuation R;
	   a = v(x)
	   b = v(y)
	   c = v(x+y)
	   a > b
	   a == c
       Text
       	   In the future, this object will be implemented at a deeper level.
	   A @TO "Module"@ object does not naturally contain a monomial order. 
	   We aim to implemenet this like we see in the object @TO "Ring"@. 
     SeeAlso
     	 leadTermValuation
	 lowestTermValuation
    

///

doc ///
     Key
    	"valM"
     Headline
         Add headline!
     Description
       Text
           Add description!
       Example
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
     SeeAlso
    

///

--------------------------------------------------------------------------------
------------------------------------ Tests -------------------------------------
--------------------------------------------------------------------------------

-- Trivial Valuation Tests
TEST ///
-- Everything should have valuation 0 except 0
val = trivialValuation
assert(val 0 == infinity)
assert(val 5 == 0)
assert(val (-9/2) == 0)
///

-- p-adic Valuation Tests
TEST ///
val = padicValuation(7)
assert(val 0 == infinity)
assert(val 7 == 1)
assert(val (-7) == 1)
assert(val (-9/7) == -1)
assert(val (7/9) == 1)
///

-- Leading Term Valuation Tests
TEST///
R = QQ[x,y]
val = leadTermValuation(R)
assert(val(x) < val(y^2))
///
TEST///
R = QQ[x,y, MonomialOrder=>Lex]
val = leadTermValuation(R)
assert(val(x) > val(y^2))
///

-- Lowest Term Valuation tests
TEST///
R = QQ[x,y,MonomialOrder => Lex];
val = lowestTermValuation R;
assert(val(x^2 + x*y) > val(y^3 + x*y^4));
///

TEST///
R = QQ[x,y,z, Degrees => {1,2,3}, MonomialOrder => GLex];
val = lowestTermValuation R;
assert(val(x^2*y^2*z + x^7*y) > val(x*y*z^2 + y^3*z))
///

-- Local ring valuation tests

-- 

end--


--------------------------------------------------------------------------------
-------------------------------- Dev Functions ---------------------------------
--------------------------------------------------------------------------------

buildPackage = x ->(
    uninstallPackage("Valuations");
    installPackage("Valuations", RunExamples => false);
    )

testPackage = x -> (
    uninstallPackage("Valuations");
    installPackage("Valuations");
    check Valuations;
    )
