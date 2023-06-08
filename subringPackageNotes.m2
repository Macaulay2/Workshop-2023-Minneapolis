-- Code snippets for Subrings:
-- these are methods that could make sense for any subring without needing a SAGBI basis
-- they work by using partial GB in the tensor ring wrt elimination order
-- they do not require a SAGBI bases (right?)


--Functions to include in the Subring package which should not require SAGBI bases

--------------------
--PRESENTATION METHOD
--------------------
-- Currently there is a presentation method in M2 but it is not designed for presenting a subring
-- maybe one should add a new option to this method to present a subring as a quotient ring
-- this should be a simple function that takes a subring, builds the presentations ring,
-- and finds the kernel of the presentation map so that one can present the subring as a quotient ring


-- Presentation Ring for a Subring
-- Presentation map = K[P_0 .. P_k] -> R : P_i -> f_i, where f_i are the gens of your subring
-- Presentation Ideal for a Subring = ker Presentation map -- compute at your own risk


-------------------
--MEMBERSHIP METHODS WITHOUT SAGBI
-------------------
-- groebnerMembershipTest(f, S) = (f lies in S)
-- f an element of the ambient ring of S
-- S a subring of a polynomial ring (or quotient ring)
-- Should this come before of after presentation of a subring?

-- 
groebnerMembershipTest = method() 
groebnerMembershipTest(RingElement, Subring) := (f, S) -> (
    subringGens := gens S;
    Q := ring subringGens;
    R := ambient Q;
    J := ideal Q;
    fLifedToR := lift(f, R);
    subringGensLiftedToR := lift(subringGens, R);
    -- construct the tensor ring
    tensorRingNumVars := (numgens R) + (numcols subringGens);
    -- Should the field below NOT be QQ automatically?
    tensorRing := QQ[Variables => tensorRingNumVars, MonomialOrder => {Eliminate(numgens R)}];    
    liftToTensorRing := map(tensorRing, R, (vars tensorRing)_{0 .. numgens R - 1});
    fInTensorRing := liftToTensorRing fLifedToR;
    subringGensInTensorRing := liftToTensorRing subringGensLiftedToR;
    JInTensorRing := liftToTensorRing J;
    I := ideal((vars tensorRing)_{numgens R .. tensorRingNumVars - 1} - subringGensInTensorRing);
    fNormalForm := fInTensorRing % (I + JInTensorRing);
    numcols selectInSubring(1, matrix {{fNormalForm}}) == 1
    )



-- groebnerSubductionQuotient(f, S) = h
-- Computes the subduction quotient of f with respect to the generators of S
-- setup: gens S = {g_1 .. g_s} in Q = K[x_1 .. x_n]/I and f in Q 
-- output: h in K[y1 .. ys] such that (f - f%S) = h(g_1 .. g_s)

groebnerSubductionQuotient = method() 
groebnerSubductionQuotient(RingElement, Subring) := (f, S) -> (
    local outputRing;
    subringGens := gens S;
    Q := ring subringGens;
    R := ambient Q;
    J := ideal Q;
    FF := coefficientRing R;
    fLifedToR := lift(f, R);
    subringGensLiftedToR := lift(subringGens, R);
    -- construct the tensor ring
    tensorRingNumVars := (numgens R) + (numcols subringGens);
    oldOrder := (options R).MonomialOrder;
    newOrder := prepend(Eliminate(numgens R), oldOrder);
    tensorRing := FF[Variables => tensorRingNumVars, MonomialOrder => oldOrder];
    liftToTensorRing := map(tensorRing, R, (vars tensorRing)_{0 .. numgens R - 1});
    fInTensorRing := liftToTensorRing fLifedToR;
    subringGensInTensorRing := liftToTensorRing subringGensLiftedToR;
    JInTensorRing := liftToTensorRing J;
    I := ideal((vars tensorRing)_{numgens R .. tensorRingNumVars - 1} - subringGensInTensorRing);
    fNormalForm := fInTensorRing % (I + JInTensorRing);
    -- output fNormalForm in the subductionQuotientRing
    outputRing = subductionQuotientRing S;
    outputMap := map(outputRing, tensorRing, matrix {toList((numgens R):0)} | vars outputRing);
    outputMap fNormalForm
    )



-- Nice shortcut
-- RingElement // Subring 
-- returns the subduction quotient
--

RingElement // Subring := (f, S) -> (
    groebnerSubductionQuotient(f, S)    
    )

-- Another nice shortcut
-- Matrix or RingElement % Subring
-- M % S, f % S
-- Returns the smallest r in the ambient ring of S
-- such that M or f = r + s for some s in S
-- Two modes of operation: 
-- 1) If Subring has a sagbi basis (stored in its cache) then use subduction
-- 2) If there is not a complete sagbi basis then use the 'extrinsic method' - see groebnerMembershipTest above
-- Note: we construct the tensor ring with a monomial order lifted from the ambient ring
--

Matrix % SAGBIBasis := (M, SB) -> (
    assert(ambient SB === ring M);
    subduction(SB, M)
    );

RingElement % SAGBIBasis := (f, SB) -> (
    assert(ambient SB === ring f);
    first first entries subduction(SB, matrix{{f}})
    );

Matrix % Subring := (M, S) -> (
    local result;
    assert(ring M === ambient S); 
    if (S#cache#?SAGBIBasis) and (S#cache#SAGBIBasis#SAGBIdata#"sagbiDone") then (
	-- S has a complete sagbi basis so use subduction
	SB := S#cache#SAGBIBasis;
	result = subduction(SB, M);	
	) else (
	-- extrinsic subduction
	subringGens := gens S;
	Q := ring subringGens;
    	R := ambient Q;
    	J := ideal Q;
    	FF := coefficientRing R;
    	MLiftedToR := lift(M, R);
    	subringGensLiftedToR := lift(subringGens, R);
    	-- construct the tensor ring
	tensorRingNumVars := (numgens R) + (numcols subringGens);
    	oldOrder := (options R).MonomialOrder;
    	newOrder := prepend(Eliminate(numgens R), oldOrder);
    	tensorRing := FF[Variables => tensorRingNumVars, MonomialOrder => oldOrder];
	liftToTensorRing := map(tensorRing, R, (vars tensorRing)_{0 .. numgens R - 1});
    	MInTensorRing := liftToTensorRing MLiftedToR;
    	subringGensInTensorRing := liftToTensorRing subringGensLiftedToR;
    	JInTensorRing := liftToTensorRing J;
    	I := ideal((vars tensorRing)_{numgens R .. tensorRingNumVars - 1} - subringGensInTensorRing);
    	MNormalForm := MInTensorRing % (I + JInTensorRing);
	projectToQ := map(Q, tensorRing, matrix {toList(numgens Q : 0_Q)} | subringGens);
        result = M - (projectToQ MNormalForm);
	);    
    result
    );

RingElement % Subring := (f, S) -> (
    first first entries (matrix{{f}} % S)
    );


----------------
--SUBRING TYPE
----------------
--Look up in the InvariantRing and the SAGBIBases packages how the Subring and InvariantRing types have been defined
--Check how subrings are handled by monomial algebras
--Decide what should be in the hastable of a subring type
--Create a subring type in the subring package that can be available independently of other packages





-------------------------------------

-- Design notes 6 June

-- Subring is a hash table in classes.m2 inside of SubalgebraBases
-- presentation should return a matrix

-- presentation B -- B should be a Subring object, subring of A
--   [f1, f2]
--  1      1
-- A <--- C 

-- presentationRing (c.f. PresRing in SubalgebraBases) (not subduction quotient)


-- 7 June (half-day)

-- projections

-- goal for AM: add Subring to Subring.m2, do some documentation, examples


-- questions - subring?

-- functions to:
    -- compute presentation map/presentation ring
    -- kernel of presentation map (slow)
    
 -- tensor ring should go in cache
 
 -- Subring then refactor
 
 -- SAGBI should load Subring package.
 
 -- mirror functionnality of Rings
 -- matrices in presentation ring (reduce? designn decisionn - could be sslow.)

-- uninstall Subring
-- values PackageDictionary
-- install Subrings -- plural
-- values PackageDictionary
