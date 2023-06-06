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
        PackageExports => {"LocalRings", "SubalgebraBases"}
        )

export{"function",
       "valuation",
       "trivialValuation",
       "padicValuation",
       "leadTermValuation",
       "lowestTermValuation",
       "domain",
       "codomain"
       }


Valuation = new Type of HashTable

valuation = method()

valuation Function := v -> (
    new Valuation from{
        function => v,
        domain => null,
        codomain => null,
        cache => new CacheTable
        }
--    internalValuation(v, null, null)
    )

ourRings := {Ring,Subring,LocalRing}

for i in ourRings do (
    for j in ourRings do (
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
        function => v,
        domain => S,
        codomain => T,
        cache => new CacheTable
        }
    )

Valuation Thing := (v,t) -> (
    if (v.domain === null) or (ring t) === v.domain then
    -- Concerns with comparing things like ZZ and QQ
    -- Concerns with subrings and local rings
    v.function t
    )

-- Trivial Valuation
--trivialValuation = symbol trivialValuation
trivialValuation = valuation (x -> if x == 0 then infinity else 0)


-- padic Valuation
getExponent = (p, x) -> (
    numFactors := 0;
    while x % p == 0 do (
        x = x // p;
        numFactors = numFactors + 1
        );
    numFactors
    )


padicValuation = method()
padicValuation ZZ := p -> (
    if not isPrime p
    then error "expected a prime integer";
    func := x -> (
        if x == 0 then infinity
        else getExponent(p, numerator x_QQ) - getExponent(p, denominator x_QQ)
        );
    valuation func
    )

--leading term valuation (max convention)
leadTermValuation = valuation (x -> if x == 0 then infinity else leadMonomial x)

lowestTermValuation = valuation (f -> if f == 0 then infinity else (sort flatten entries monomials f)_0 )

-- localRingValuation = valuation (f -> if f == 0 then infinity else


OrderedQQn = new Type of Module
OrderedQQVector = new Type of Vector

orderedQQn = method()
orderedQQn(ZZ, List) := (n, monOrder) -> (
    R := QQ[Variables => n, MonomialOrder => monOrder];
    ordMod := new OrderedQQn of OrderedQQVector from QQ^n;
    ordMod.cache.Ring = R;
    ordMod
    )

OrderedQQn == OrderedQQn := (N, M) -> (
    N.cache.Ring === M.cache.Ring
    )

OrderedQQn#{Standard,AfterPrint} =
OrderedQQn#{Standard,AfterNoPrint} = M -> (
    << endl; -- double space
    << concatenate(interpreterDepth:"o") << lineNumber;
    << " : Ordered QQ^" | toString (numgens M) | " module" << endl
    );

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

    if aMonomial < bMonomial then symbol <
    else if aMonomial > bMonomial then symbol >
    else symbol ==
    )

---Documentation
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
         MethodFunction
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
           Stuff goes here
       Example
           v = padicValuation 7;
           v 98
           v (2/7)
           v 0
           v (-42)
     SeeAlso
         MethodFunction
     ///

doc ///
     Key
        "lowestTermValuation"
     Headline
        The valuation which returns the lowest term of an element of an ordered ring
     Usage
         v = lowestTermValuation
     Outputs
         v:Valuation
             the lowest term valuation
     Description
       Text
           This valuation returns the lowest (trailing) term of a polynomial with respect to the ring's term order.
       Example
           R = QQ[a,b,c, MonomialOrder => Lex];
           v = lowestTermValuation;
           f = 13*a^2*b + a*c^3;
           v f
           S = QQ[a,b,c, MonomialOrder => RevLex, Global => false];
           f = 13*a^2*b + a*c^3;
           v f
     SeeAlso
         MethodFunction
     ///

     TEST ///
         assert(trivialValuation 5 == 0)
         ///
