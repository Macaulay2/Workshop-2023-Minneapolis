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
       "valuation"}

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
