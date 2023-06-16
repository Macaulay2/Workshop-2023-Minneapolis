
-- This is a method for obtaining the rank of a global algebra k[x_1..x_n] over a field, given the input of a zero-dimensional ideal (f_1,..f_n) < k[x_1..x_n]
rankGlobalAlgebra = method()
rankGlobalAlgebra (List) := (ZZ) => (Endo) -> (
    -- Get the underlying field    
    kk := coefficientRing(ring(Endo#0));    
    if isField(kk) == false then(
    	kk = toField(kk);
    	);
    
    -- Let S = k[x_1..x_n] be the ambient polynomial ring
    S:=ring(Endo#0);
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then (
	print "Error: ideal is not zero-dimensional";
	return Endo;
	);
    
    -- Get the rank of S/ideal(Endo) as a kk-vector space
    return numColumns(basis(S/ideal(Endo)));
    
    );


-- TESTS

-- SQQ = QQ[x,y]
-- LQQ = {x^2 - y, y-2}    
-- assert(rankGlobalAlgebra(LQQ)==2)

-- SGF = GF(101)[x,y]
-- LGF = {x^2 - y, y-2}    
-- assert(rankGlobalAlgebra(LGF)==2)

-- SRR = RR[x,y]
-- LRR = {x^2 - y, y-2}    
-- assert(rankGlobalAlgebra(LRR)==2)

-- SCC = CC[x,y]
-- LCC = {x^2 - y, y-2}    
-- assert(rankGlobalAlgebra(LCC)==2)



--(entries(gens ideal(LQQ)))_0
