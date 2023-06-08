-- -*- coding: utf-8 -*-
newPackage(
    "Subrings",
    Version => "1.0",
    Date => "June 6, 2023",
    Authors => {
	{Name => "Casey Hill", Email => "caseybhill2@gmail.com"},
        {Name => "Trevor K. Karn", Email => "karnx018@umn.edu"},
        {Name => "Miranda Moore", Email => "moor2340@umn.edu"},
        {Name => "Christopher O'Neill", Email => "cdoneill@sdsu.edu"}},
    Headline => "a package for subrings",
    Keywords => {"Documentation"},
    DebuggingMode => true
    )


export {"Subring",
        "subring",
	"presentationRing",
	"presentationMap",
        "presentationIdeal",
        "toQuotientRing"}

Subring = new Type of HashTable

-- a method to create subrings
subring = method()
subring Matrix := genMatrix -> (
    
    -- compute presentation ring
    R := ring genMatrix;
    nGens := numgens source genMatrix;
    k := coefficientRing R;
    x := symbol x;
    P := k[x_1..x_nGens];
    
    -- compute presentation map 
    f := map(R, P, genMatrix);
    
    S := new Subring from {
	generators => genMatrix,
	ambient => R,
	
	-- presentation ring: one variable for each generator
	presentationRing => P,
      
	-- presentation map: presentation ring --> ambient ring
	presentationMap => f,
	 
	-- presentation ideal?? -- probably compute into cache 
	cache => new CacheTable
	};
    S
    )

subring List := genList -> (
    subring matrix {genList}
    )

presentationRing = method()
presentationRing Subring := S -> (
    S#presentationRing
    )

presentationMap = method()
presentationMap Subring := S -> (
    S#presentationMap
    )

presentationIdeal = method()
presentationIdeal Subring := S -> (
    P := presentationRing S;
    f := presentationMap S;
    return ker f; --kernel is cached automatically
    )

toQuotientRing = method()
toQuotientRing Subring := S -> (
    P := presentationRing S;
    I := presentationIdeal S;
    return P/I;
    )


beginDocumentation()
doc ///
 Node
  Key
   Subrings
  Headline
   a package to deal with subrings
  Description
   Text
     {\em Subrings} is a package to give basic subroutines for subrings.
    Caveat
     There are other subring flavor things out there. 
    Subnodes
     Subring
     subring
     presentationRing
     presentationMap
     presentationIdeal
     toQuotientRing
///
doc ///
    Key
        subring
        (subring, Matrix)
        (subring, List)
    Headline
         Construct a subring of a polynomial ring
    Usage
        S = subring M
        S = subring L
    Inputs
        M:Matrix
            of generators for a subring of a @ ofClass{PolynomialRing} @
        L:List
            of generators for a subring of a @ ofClass{PolynomialRing} @
    Outputs
        S:Subring
            the subring of the polynomial ring
    Description
        Text
            An easy way to specify a subring is to specify the ambient ring and the generators of a desired subring as a matrix. The ambient ring is implicit in the constructor function.
        Example
            R = QQ[x,y]
            M = matrix(R, {{x^2, x*y, y^2}})
            S = subring M
        Text
            This function also accepts a list of elements of a polynomial ring as an input.
        Example
            R = QQ[x,y]
            L = {x^2, x*y, y^2}
            S = subring L
///
doc ///
    Key
        presentationIdeal
        (presentationIdeal, Subring)
    Headline
        Compute the presentation ideal of a subring
    Usage
        I = presentationIdeal S
    Inputs
        S:Subring
            a @ ofClass{Subring} of a @ ofClass{PolynomialRing}
    Outputs
        I:Ideal
            the presentation ideal of the subring
    Description
        Text
            This function finds the presentation ideal of the subring, which is defined to be the presentation ring modulo the presentation map.
        Example
            R = QQ[x,y]
            S = subring {x^2, x*y, y^2}
            I = presentationIdeal S
///

TEST ///
    assert(true == true)
///

end--

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
