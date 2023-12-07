newPackage("Valuations",
        Headline => "implementation of valuations for rings",
        Version => "0.1",
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
        PackageExports => {"LocalRings", "SubalgebraBases", "InvariantRing"}
--      PackageExports => {"SubalgebraBases"}
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
       "orderedQQn"
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
    valFunc := f -> if f == 0_R then infinity else monomialToOrderedQQVector(leadTerm f, orderedMod);
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
