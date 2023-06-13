-*

   Copyright 2023, Oliver Clarke, Francesca Gandini,
   Casey Hill, Trevor Karn, Miranda Moore, Chris O'Neill
    
   You may redistribute this file under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 2 of
   the License, or any later version.
*-


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
    Key
    	Subring
    Headline
    	the class of finitely generated subrings of polynomial rings
    Description
    	A subring of a @ ofClass{PolynomialRing} @ is a ring with
	unity contained inside of another ambient ring that is closed
	under the operations of the ambient ring.
    Example
    	R = QQ[x,y];
	L = {x^2, y^2};
	S = subring L
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
            a @ ofClass{Subring} @ of a @ ofClass{PolynomialRing} @
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
