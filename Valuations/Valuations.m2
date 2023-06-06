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

localRingValuation = valuation (f -> if f == 0 then infinity else 


---Documentation
beginDocumentation()

doc ///
     Key
     	 triviaValuation
	 (triviaValuation)
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
	   v -13
	   v 100000000
	   v 14/23
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
