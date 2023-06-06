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
        Configuration => {}
        )

export{"function",
       "valuation",
       "trivialValuation",
       "padicValuation",
       "leadTermValuation",
       "lowestTermValuation"
       }

Valuation = new Type of HashTable

valuation = method()
valuation Function := v -> (
    new Valuation from{
        function => v,
        source => null,
        target => null,
        cache => new CacheTable
        }
    )

Valuation Thing := (v,t) -> (v.function t)

-- Trivial Valuation
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
    assert isPrime p;
    func := x -> (
        if x == 0 then infinity
        else getExponent(p, numerator x_QQ) - getExponent(p, denominator x_QQ)
        );
    valuation func
    )

--leading term valuation (max convention)
leadTermValuation = valuation (x -> if x == 0 then infinity else leadMonomial x)

lowestTermValuation = valuation (f -> if f == 0 then infinity else (sort flatten entries monomials f)_0 )


---Documentation
beginDocumentation()

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
    	lowestTermValuation     	
     Headline
    	The valuation which returns the lowest term of an element of an ordered ring
     Usage
     	 v = lowestTermValuation
     Inputs
     	 null
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
     	 MethodFunction, leadTermValuation
     ///
