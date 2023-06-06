restart
debug needsPackage "Brackets"
-*

The goal of this exercise is to prove that 4 lines in P^3 can have either 1, 2, or infinitely many transversals. 

To do this, we fix 4 pairs of points in P^3, and look at the condition that the common transversal of the first three lines
also intersects with the other line.  

The first part of this script is to be converted as the Grassmann option for the GCExpression method. 

*-
----------------------------------------------------------------

-- grassmann method for computing GCalgebra 
G37 = Grassmannian(3,7) --ideal of Pluecker relations

G37gb = gb G37  -- Groebner basis

G37ring = ring G37 -- grassmannian ring before quotient

G = G37ring/G37 -- quotient ring is desired GCalgebra
--G = gc(a .. h, 4) -- old method that is too computationally expensive

----------------------------------------------------------------
R = (G, [u,v]) -- algebra with indeterminants u and v over GCalgebra G

l1 = a*b; --
l2 = c*d;
l3 = e*f;
--l4 = g*h;

l = (((u*a + v*b)*c*d)^(e*f)) * (u*a + v*b)

quantity1 = l*g*h

toBracketPolynomial(quantity1, G)

