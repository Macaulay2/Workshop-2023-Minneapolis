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
generators Subring := S -> S#generators
ambient  Subring := S -> S#ambient

-- a method to create subrings
subring = method()
subring Matrix := genMatrix -> (
    
    -- compute presentation ring
    R := ring genMatrix;
    nGens := numgens source genMatrix;
    k := coefficientRing R;
    p := symbol p;
    P := k[p_1..p_nGens];
    
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

-- a quotient ring isomorphic to the image of the subring inside of the presentation ring
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
    Example
     needsPackage("Subrings")
///

doc ///
 Node
  Key
   presentationRing
  Headline
   a polynomial ring with a variable for each subring generator
  Description
   Text
    The {\tt presentationRing} of a subring {\tt S} is a polynomial ring
    with a variable for each generator of {\tt S}.
   Example
    R = QQ[x,y];
    g = {x^2, x*y, y^2};
    S = subring g
    presentationRing S
///

doc ///
 Node
  Key
   presentationMap
  Headline
   the map from the {\tt presentationRing} to a subring {\tt S}
  Description
   Text
    There is a map sending each generator of the presentation ring to
    the corresponding element of the ambient ring. This is that map.
   Example
    R = QQ[x,y];
    g = {x^2, x*y, y^2};
    S = subring g;
    f = presentationMap S
    P = presentationRing S
    p = P_0 * P_1 - P_2
    f(p)
///


TEST ///
    assert(true == true)
///

end--

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
