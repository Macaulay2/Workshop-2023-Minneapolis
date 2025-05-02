newPackage(
"NeuralIdeals",
Version => "1.00",
Date => "June 5, 2023",
Authors => {{Name => "Hugh Geller"},{Name => "Rebecca R.G."}},
Headline => "canonical forms of neural ideals",
Keywords => {"Coding Theory", "Combinatorial Commutative Algebra", "Commutative Algebra"},
DebuggingMode => false,
PackageImports => {"PrimaryDecomposition","PseudomonomialPrimaryDecomposition"},
Reload => false
)

--************************************************************************************************************
--************************************************************************************************************
--***Given a neural code, this package computes the corresponding neural ideal and canonical form of the neural ideal.
--************************************************************************************************************
--************************************************************************************************************
--***Acknowledgements:                                                                                     ***
--***Special thanks to Juliette Bruce for contributing the original code for allCodeWords and neuralCodeComplement.***
--************************************************************************************************************
--************************************************************************************************************

export{--types
    "NeuralCode",
    --methods/functions
    "neuralCode",
    "polarizedRing",
    "neuralIdeal",
    "canonicalForm",
    "codeSupport",
    "neuralIdealToCode",
    "isPseudomonomial",
    "receptiveFieldRelation",
    "polarizePseudomonomial",
    "polarizedCanonicalForm",
    "polarizedCanonicalIdeal",
    "polarizedCanonicalResolution",
    "depolarizationMap",
    "canonicalResolution",
    "isCanonical",
    --Symbols
    "codeWords",
    "Factor",
    "Iterative",
    "SharedIndex",
    "Polarized"}

--protect codeWords
protect dimension
--protect Factor
--protect Iterative
--protect SharedIndex
--protect Polarized

--creates a ring with n or 2n variables
createRing = method()

createRing(ZZ,String) := Ring => (n,z) -> (
    x := getSymbol z;
    R := (ZZ/2)(monoid[x_1..x_n])
    )

createRing(ZZ) := Ring => n -> (
    createRing(n,"x")
    )

createPolarizedRing = method()

createPolarizedRing(ZZ,String,String) := Ring => (n,z,w) -> (
    x := getSymbol z;
    y := getSymbol w;
    S := (ZZ/2)(monoid[x_1..x_n,y_1..y_n])
    )

createPolarizedRing(ZZ) := Ring => n -> (
    createPolarizedRing(n,"x","y")
    )

--type that will store the data of a neural code
NeuralCode = new Type of HashTable
NeuralCode.synonym = "neural code"

--constructs a neural code object from a list of codewords given as binary strings of the same length
neuralCode = method()

--add this in in case don't want to add polarized data, default is not to
--reason: what if want to give it only 1 ring, or 1 variable symbol? (say if not working with the polarization at all)
--create an option to add a polarization, or an option not to?
--neuralCode = method(Options => {
--	Polarized => false
--	})

neuralCode (List,Ring,Ring) := NeuralCode => (codeList,R,S) -> (
    d := #(codeList#0);
    X:=new NeuralCode from {
	symbol codeWords => codeList,
	symbol dimension => d,
	symbol cache => new CacheTable,
	};
    X.cache.ring = R;
    X.cache.polarizedRing = S;
    X
    )
--need to add a check that S contains R, has twice as many variables

neuralCode (List,String,String) := NeuralCode => (codeList,z,w) -> (
    d := #(codeList#0);
    R := createRing(d,z);
    S := createPolarizedRing(d,z,w);
    neuralCode(codeList,R,S)
    )

neuralCode List := NeuralCode => codeList -> (
    neuralCode(codeList,"x","y")
    )

--short way to get the dimension of a neural code 
dim NeuralCode := C -> C.dimension

--short way to get the ring of a neural code
ring NeuralCode := C -> C.cache.ring

--short way to get the polarized ring of a neural code
polarizedRing = method()

polarizedRing NeuralCode := Ring => C -> C.cache.polarizedRing

--checks whether a NeuralCode is well-defined
isWellDefined NeuralCode := Boolean => X -> (
    --check keys
    K:=keys X;
    expectedKeys := set{symbol codeWords, symbol dimension}; 
    if set K =!= expectedKeys then (
	if debugLevel > 0 then (
	    added := toList(K - expectedKeys);
	    missing := toList(expectedKeys - K);
	    if #added > 0 then
	    << "-- unexpected key(s): " <<toString added << endl;
	    if #missing >0 then
	    << "-- missing key(s): " <<toString missing << endl
	    );
	return false
	);
    -- check types
    if not instance(X.codeWords, List) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a list" <<endl;
	return false
	);
    if X.codeWords === {} or not all (X.codeWords, r->instance(r,String)) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a nonempty list of strings" <<endl;
	return false
	);
    if not all (X.codeWords, r->all(r,i->(value(i)==0 or value(i)==1))) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a list of strings of 0's and 1's" << endl;
	return false
	);
    codeList := codeWords X;
    d:= # (codeList#0);
    if not all (X.codeWords, r-> #r === d) then (
	if debugLevel > 0 then
	<< "-- expected 'codes' to be a list of equal length strings" << endl;
	return false
	);
    --if codeList == {} then (
	--if debugLevel > 0 then
	--<< "--expected 'codes' to be a nonempty list" <<endl;
--	return false
	--);
    if dim X != numgens ring X then (
	if debugLevel >0 then
	<< "-- expected dimension of ring to equal length of code words" << endl;
	return false
	);
    if numgens polarizedRing X != 2*(numgens ring X) then (
	if debugLevel >0 then
	<< "--expected dimension of polarized ring to be twice dimension of first ring" << endl;
	return false
	);
    true);

--given a neural code, this constructs a ring for polarizations of the neural ideal to live in
--do I want to create this when I create the neural code? probably eventually yes.
--polarizedRing = method();

--polarizedRing(NeuralCode,String,String) := Ring => (C,z,w) -> (
--    d := dim C;
--    x := getSymbol z;
--    y := getSymbol w;
--    S := (ZZ/2)(monoid[x_1..x_d,y_1..y_d]);
--    C.cache.polarizedRing = S;
--    S
--    )

--polarizedRing(NeuralCode) := Ring => C -> (
--    polarizedRing(C,"x","y")
--    )

--gives a list of all code words on a given number of neurons
--used internally in the neuralIdeal function
allCodeWords = method();
allCodeWords ZZ := List => d ->(
    L1 := apply(d+1,i->(
	    apply(i,i->1)|apply(d-i,j->0)
	    ));
    L2 := unique flatten apply(L1,i->permutations i);
    apply(L2, i-> concatenate(apply(i,j->toString j)))
    )

--given a neural code, gives the list of code words not in it
--used internally in the neuralIdeal function
neuralCodeComplement = method();
neuralCodeComplement NeuralCode := List => C ->(
    d := dim C;
    L1 := allCodeWords(d);
    L:=C.codeWords;
    sort(toList(set(L1)-set(L))) --may not need to sort
    --for i in L do L1=delete(i,L1);
    --L1
    )    

--gives the neural ideal of a neural code by the method of Curto, Itskov, et al
neuralIdeal = method();

neuralIdeal NeuralCode := Ideal => C -> (
    R := ring C;
    oppC:=neuralCodeComplement C;
    genList := for a in oppC list (
	prodList := for j to d-1 list (
	    if a#j == "1" then R_j else (1-R_j)
	    );
	product(prodList)
	);
    --genList:=for i to #oppC-1 list (
    	--prod:=1;
    	--for j to d-1 do
	    --prod=prod*(1-value((oppC#i)#j)-R_j);
	--prod
	--);
    ideal genList
    )

--computes the canonical form of a neural code by the iterative method of Petersen, Youngs, et al
--used internally as the default option for canonicalForm(NeuralCode)
iterCanonicalForm = method()

iterCanonicalForm NeuralCode := List => C -> (
    R := ring C;
    initialCodeWord := C.codeWords#0;
    currentGens := for i to d-1 list (
	R_i-value(initialCodeWord#i) --creates canonical form for single codeword as a starting list without using append
	);
    for i from 1 to #C.codeWords - 1 do (
	currentCodeWord := C.codeWords#i;
	codeCoordinates := for j to d-1 list(
	    value(currentCodeWord#j)
	    );
	factors := apply(0..(d-1),j->(R_j-codeCoordinates#j)); --or would doing vars R - codeCoordinates be better?
	substitutionMap := map(R,R,codeCoordinates);
	H := partition(gen -> substitutionMap(gen)==0,currentGens,{true,false});
	keepList := H#true;
	changeList := H#false;
	newList := flatten (
	    for elem in changeList list (
		for comp in factors list (
		    if elem%(comp-1) == 0 then continue;
		    goToNext := false;
		    g := elem*comp;
		    for mgen in keepList do (
		    	if g%mgen == 0 then (goToNext = true;
		    	    break)
		    	);
		    if goToNext then continue;
		    g
		    )
	    	)
	    );
	currentGens = join(keepList,newList);
	);
    currentGens
    )

--original algorithm for the canonical form of a pseudomonomial ideal from Curto, Itskov et al
--does everything but remove gens divisible by another gen
primaryDecompositionAlmostCanonicalForm = method()

primaryDecompositionAlmostCanonicalForm Ideal := List => I -> ( 
    decomp := primaryDecompositionPseudomonomial I;
    multipliedGens :=product(decomp);
    R := ring I;
    d := numgens R;
    booleanIdeal := ideal(apply(gens R,g -> g*(1-g)));
    --booleanIdeal := ideal(apply(d,i->(R_i*(1-R_i))));
    booleanR := R/booleanIdeal;
    reducedGens := promote(multipliedGens,booleanR);
    noZeroGens := compress gens reducedGens;
    --reducedGens := apply(first entries gens multipliedGens,i->sub(i,booleanR));
    --noZeroGens := delete(sub(0,booleanR),reducedGens);
    almostGens := unique first entries lift(noZeroGens,R)
    --almostGens := unique apply(noZeroGens,i->(sub(i,R)))
    )

--functions needed to implement Geller-R.G. algorithm for canonical form
isSharedIndex = method()

isSharedIndex (RingElement,RingElement,ZZ,Ring) := Boolean => (g,h,i,R) -> (
--    R:=ring g;
--    if ring g =!= ring h then error "Expected two elements from the same ring";
    if i > dim R then error "Expected index at most the dimension of the ring";
    if i < 1 then error "Expected index at least 1";
    x:=R_(i-1);
    (g*h)%(x*(1-x))==0
    )

isUniqueSharedIndex = method()

isUniqueSharedIndex (RingElement,RingElement,ZZ,Ring) := Boolean => (g,h,i,R) -> (
    n := numgens R;
    if isSharedIndex(g,h,i,R) then (
	onlySharedIndex := true;
	for j from 1 to n when onlySharedIndex do (
	    if j == i then continue;
	    if isSharedIndex(g,h,j,R) then (
		onlySharedIndex = false;
		break
		);
	    );
	onlySharedIndex
	)
    else false
    )

newGens = method()

newGens (List,ZZ,Ring) := List => (listGens,i,R) -> ( 
    unique flatten (for g in listGens list (
	for h in listGens list (
	    if h==g then continue;
	    if isUniqueSharedIndex(g,h,i,R) then lcm(g,h)//(R_(i-1)*(1-R_(i-1))) else continue
	    )
	)
    ) )

--produces the almost canonical form of an ideal I using the shared index method of Geller-R.G.
almostCanonicalForm = method()

almostCanonicalForm Ideal := List => I -> (
    R := ring I;
    n := numgens R;
    listGensI := first entries gens I;
    for i from 1 to n do (
	listGensI=join(listGensI,newGens(listGensI,i,R))
	);
    unique listGensI
    )

--removes generators divisible by another generator to get from almost canonical form to canonical form
removeGens = method()

removeGens List := List => almostGens -> (
    for i in almostGens list (
	isDivisible := false;
	for j in almostGens do (
	    if i%j==0 and i =!= j then (isDivisible=true; break));
	if isDivisible then continue;
	i
	)
    )

sharedIndexCanonicalForm = method()

sharedIndexCanonicalForm Ideal := List => I -> (
    removeGens(almostCanonicalForm(I))
    )

---------

--exported function to compute the canonical form of a neural ideal or neural code
--can decide to display it factored
--default method for a neural code is the iterative method
--default method for an ideal is the primary decomposition method
canonicalForm = method(
    Options => {
	Factor => false,
	Iterative => true, --made iterative the default,
	SharedIndex => false
	})

canonicalForm Ideal := List => opts -> I -> (
    if not isSquarefreePseudomonomialIdeal(I) then error "Expected a squarefree pseudomonomial ideal.";
    if opts.SharedIndex then (
	canon := sharedIndexCanonicalForm(I);
	if opts.Factor then apply(canon,factor) else canonForm
	)
    else (
	canonP := removeGens(primaryDecompositionAlmostCanonicalForm(I));
	if opts.Factor then apply(canonP,factor) else canonPForm
	)
    )
--TO DO: throw error if ideal is not pseudomonomial

canonicalForm NeuralCode := List => opts -> C -> (
    if opts.Iterative then (
	C.cache.canonicalForm = iterCanonicalForm(C);
	if opts.Factor then apply(C.cache.canonicalForm,factor) else
	C.cache.canonicalForm
	)
    else
    C.cache.canonicalForm = canonicalForm(neuralIdeal(C));
    if opts.Factor then apply(C.cache.canonicalForm,factor) else C.cache.canonicalForm
    )

--finds the support of a given neural code (list of sets of neurons that fire together)
codeSupport = method();
codeSupport NeuralCode := List => C -> (
    fullSupport := for c in C.codeWords list (
	cSupport := for i to #c-1 list (if value(c#i) == 0 then continue; i+1)
	)
    )

--outputs the neural code of a list of pseudomonomials, usually the canonical form of a neural ideal
canonicalFormToCode = method();

canonicalFormToCode List := NeuralCode => L -> (
    R := ring L#0;
    d := numgens R;
    --checks that entries in list are squarefree pseudomonomials
    if not isSquarefreePseudomonomialIdeal(ideal(L)) then error "Expected elements that generate a squarefree pseudomonomial ideal.";
    --checks that generators don't generate the unit ideal
    if ideal(L)==sub(ideal(1),R) then error "Expected generators of a non-unit ideal.";
    --checks that all elements in list are in the same ring
    for ell in L do (if ring ell =!= R then error "Expected elements of the same ring.");
    allCodes := allCodeWords(d);
    codeList := for i in allCodes list (
	validCode := true;
	for j in L do (
	    M:=matrix{apply(d,k->sub(value(i#k),R))};
	    if sub(j,M) != 0 then (
		validCode = false;
		break
		);
	    );
	if not validCode then continue else i
	);
    neuralCode codeList
    )

----The following function is an internal function from the PseudomonomialPrimaryDecomposition package by Alan Veliz-Cuba

-- determines if a polynomial is square free pseudomonomial
-- Input:
-- Polynomial P in bitwise form
-- Output:
-- true or false
isPseudomonomial = method();
isPseudomonomial RingElement := Boolean => P -> ( 
    -- check if polynomial is a unit or zero
    if P == 0 then return false;
    if isUnit P then return true;
    -- factor polynomial P and initialize the support list
    factoredP := factor P;
    allSupport := {};
    -- test if some factor is not of the form (xi-a) where a=0 or 1
    for i to #factoredP-1 do ( 
        -- evaluate ith factor
        base := value factoredP#i; 
        -- if factor is not a unit but is a constant -> not a square free pseudomonomial
        if isUnit base then continue;
        if isConstant base then return false;
        -- find if factor is equal to xi or xi-1
        suppi := support base;
        if #suppi >= 2 then return false;
        if suppi_0 =!= base and suppi_0-1 =!= base then return false;
        allSupport = append(allSupport,suppi_0);
    );
    -- find if there are factors xi, xi-1 simultaneously -> not a square free pseudomonomial
    #(support P) == #allSupport
    -- if #(support P) != #allSupport then return false;
    -- true
)

--------------------------------------

--input a pseudomonomial, return the lists corresponding to the receptive field relation, i.e. {sigma,tau}
--in other words {x's dividing pseudomonomial,(1-x)'s dividing pseudomonomial}
receptiveFieldRelation = method();

receptiveFieldRelation(RingElement) := List => P -> (
    if isPseudomonomial(P) == false then error "Expected input to be a squarefree pseudomonomial";
    R := ring P;
    d := numgens R;
    H := partition(i -> (P%R_(i-1)==0,P%(1-R_(i-1))==0),toList(1..d),{(true,true),(true,false),(false,true)});
    sigma := flatten{H#(true,true),H#(true,false)}; 
    tau := flatten{H#(true,true),H#(false,true)}; 
    {sigma,tau}
    )
	

--input a pseudomonomial, outputs the polarization.
--can specify the ring it comes from and the ring it goes to
--or the ring it goes to (recommended at least this)
--or no rings and it will create them
polarizePseudomonomial = method();


polarizePseudomonomial(RingElement,Ring) := RingElement => (P,S) -> (
    if not isPseudomonomial(P) then error "Expected input to be a Pseudomonomial";
    if (numgens S)%2 != 0 then error "Ring must have an even number of generators";
    if 2*(numgens ring P) > numgens S then error "Target ring does not have enough generators for polarization";
    d := (numgens S)//2;
    st := receptiveFieldRelation(P);
    sigma := st_0;
    tau := st_1;
    use S;
    mon := 1_S;
    for i in sigma do mon = mon*S_(i-1);
    for i in tau do mon = mon*S_(d+i-1);
    mon
    )

--want to be able to specify variable to create new ring without making ring
--currently stuck
--polarizePseudomonomial (RingElement,String) := RingElement => (P,z) -> (
--    R := ring P;
--    d := numgens R;
--    x := getSymbol z;
--    S := R[x_1..x_d]
--    )

polarizePseudomonomial RingElement := RingElement => P -> (
    R := ring P;
    d := numgens R;
    S := createPolarizedRing(d,x,y);
    polarizePseudomonomial(P,S)
    )

--polarizes the elements of a list into a particular polynomial ring
--applied internally to polarize canonical forms by work of Gunturkun, Jeffries, and Sun
polarizeList = method()

polarizeList(List,Ring) := List => (L,S) -> (
    for P in L list polarizePseudomonomial(P,S)
    )

--given a neural code, produces the polarized canonical form in S
polarizedCanonicalForm = method()

polarizedCanonicalForm NeuralCode := List => C -> (
    S := polarizedRing C;
    C.cache.polarizedCanonicalForm = polarizeList(canonicalForm(C),S);
    C.cache.polarizedCanonicalForm
    )

polarizedCanonicalForm(Ideal,Ring) := List => (I,S) -> (
    L := canonicalForm(I);
    polarizeList(L,S)
    )

--need to be able to create the polarized ring of an ideal for this to work
--see if can fix later
--polarizedCanonicalForm(Ideal) := List => I -> (
--    S :=polarizedRing(I);
--    polarizedCanonicalForm(I,S)
--    )

--add code that if an ideal is already polarized, can still get the canonical form


--given a neural code, gives the ideal of S generated by its polarized canonical form
--ask Mike if this is worth having or not
polarizedCanonicalIdeal = method()

polarizedCanonicalIdeal(NeuralCode,Ring) := Ideal => (C,S) -> (
    ideal(polarizedCanonicalForm(C,S))
    )

polarizedCanonicalIdeal(NeuralCode) := Ideal => C -> (
    S := polarizedRing C;
    ideal(polarizedCanonicalForm(C,S))
    )

polarizedCanonicalIdeal(Ideal,Ring) := Ideal => (I,S) -> (
    ideal(polarizedCanonicalForm(I,S))
    )

--also need a polarizedRing of an ideal for this to work
--polarizedCanonicalIdeal(Ideal) := Ideal => I -> (
--    S := polarizedRing I;
--    ideal(polarizedCanonicalForm(I,S))
--    )

--given a pseudomonomial (or squarefree monomial) ideal, determines whether it's in canonical form
--issue: will computing the canonical form respect order such that this will work?
--issue: would be nice to have this be faster than computing the canonical form, but currently it takes just as long
isCanonical = method(
    Options => {
	Polarized => false
	}
    );

isCanonical Ideal := Boolean => opts -> I -> (
    if opts.Polarized then (
	polarizedCanonicalForm(I) == first entries gens I
	)
    else (
	canonicalForm(I) == first entries gens I
	)
    )
--make this cache the canonical form too? or use cached one if it exists?


-------------------------------------------------

--given a neural code, computes its canonical form, polarizes it, and computes a minimal resolution
--can input the polarizedRing or not
--note that the polarized canonical form is a set of minimal generators, so res will give a minimal resolution
polarizedCanonicalResolution = method();

polarizedCanonicalResolution (NeuralCode) := Resolution => C -> (
    S := polarizedRing C;
    L := polarizedCanonicalIdeal(C,S);
    res L
    )

--polarizedCanonicalResolution (NeuralCode,Ring) := Resolution => (C,S) -> (
--    L := polarizedCanonicalIdeal(C,S);
--    res L
--    )

--polarizedCanonicalResolution(NeuralCode) := Resolution => C -> (
--    d := dim C;
--    x := getSymbol "x";
--    y := getSymbol "y";
--    S := (ZZ/2)(monoid[x_1..x_d,y_1..y_d]);
--    polarizedCanonicalResolution(C,S)
--    )

--sets up a depolarization map from the polarized ring to a polynomial ring in half the variables
depolarizationMap = method();

depolarizationMap(Ring,Ring) := (R,S) -> ( ----Target ring followed by source ring
    if 2*(numgens R) < numgens S then error "Target ring must have at least half the number of generators of the source";
    if (numgens S)%2 != 0 then error "Source ring must have an even number of generators";
    d := (numgens S)//2;
    maintain := for i to d-1 list R_i;
    change := for i to d-1 list 1+R_i;
    depolarizationList := maintain|change;
    dePolMap := map(R,S,depolarizationList)
    )

--uses the depolarization map and polarizedCanonicalResolution to create the canonical resolution of a neural code
canonicalResolution = method();

canonicalResolution NeuralCode := Resolution => C -> (
    polarRes := polarizedCanonicalResolution(C);
    depolarMap := depolarizationMap(ring C,polarizedRing C);
    depolarMap(polarRes)
    )

--canonicalResolution (NeuralCode,Ring,Ring) := (C,R,S) -> (
--    polarRes := polarizedCanonicalResolution(C,S);
--    depolarMap := depolarizationMap(R,S);
--    depolarMap(polarRes)
    --d := dim C;
    --quotientIdeal := ideal(for i to d-1 list (S_i+S_(d+i)-1));
    --R := S/quotientIdeal;
    --polarRes ** R
--    )

--canonicalResolution(NeuralCode,Ring) := (C,R) -> (
--    S := polarizedRing(C);
--    canonicalResolution(C,R,S)
--    )

--canonicalResolution(NeuralCode) := C -> (
--    S := polarizedRing(C);
--    R := ring C;
--    canonicalResolution(C,R,S)
--    )
    
beginDocumentation()

document{
  Key => NeuralIdeals,
  Headline => "neural ideals",
  EM "NeuralIdeals", " is a package that allows computation of a neural ideal or its canonical form from a neural code",
  Caveat => "In progress"
  }

document{
    Key => {NeuralCode},
    Headline => "NeuralCode -- type of a neural code",
    Description => "Stores the data of a neural code, specifically its code words, dimension, ring, and polarized ring."
    }

document{
  Key => {neuralCode},
  Headline => "neuralCode -- creates a NeuralCode object",
  Usage => "neuralCode(L) or neuralCode(L,s,t) or neuralCode(L,R,S)",
  Inputs => {"L, a list of binary strings of the same length like 000 and 101","s and t, strings to names the variables in the ring of C and the polarized ring of C","R and S, the ring of C and polarized ring of C"},
  Outputs => {"a NeuralCode"},
  Description => "Create a NeuralCode from a list of binary strings. By default, Macaulay2 will choose the ring and polarized ring of the code, but these can be specified by giving strings for the variable names or inputting rings.",
  EXAMPLE lines ///
  C=neuralCode({"000","001","101"});
  ring C
  polarizedRing C
  ///,
  EXAMPLE lines ///
  C=neuralCode({"00","10"},"z","w");
  ring C
  polarizedRing C
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1,x_2];
  S=ZZ/2[x_1,x_2,y_1,y_2];
  C=neuralCode({"00","10","11"},R,S);
  ring C
  polarized ring C
  ///  
  }

document{
    Key => {polarizedRing},
    Headline => "Polarized Ring",
    TEX "Gives the ring of the polarized form of the neural ideal.",
    Usage => "polarizedRing(neuralCode(code))",
    Inputs => {"A NeuralCode or ideal in a polynomial ring (typically a neural ideal)"},
    Outputs => {"A polynomial ring over ZZ/2 in 2n variables, where n is the number of neurons."},
    TEX "We give some examples",
    EXAMPLE lines ///
    polarizedRing(neuralCode("000","001","101")
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    I=ideal(x_1*x_2,x_1*x_3);
    polarizedRing(I)
    ///
    }

document{
  Key => {neuralIdeal, (neuralIdeal,NeuralCode),(neuralIdeal,NeuralCode,Ring)},
  Headline => "Neural ideal.",
  TEX "A method which computes the neural ideal for a given neural code (not necessarily in canonical form) by the method of Curto, Itskov, et al in The Neural Ring.",
  Usage => "neuralIdeal(neuralCode(code)) or neuralIdeal(neuralCode(code),Ring)",
  Inputs => {"neuralCode or neuralCode,Ring"},
  Outputs => {"The neural ideal corresponding to the given neural code, in the given Ring or in ZZ/2[x_1..x_d] where d is the dimension of the neural code"},
  TEX "We compute an example",
  EXAMPLE lines ///
  C=neuralCode("000","001");
  neuralIdeal(C)
  ///,
  EXAMPLE lines ///
  C=neuralCode("000","001");
  R=ZZ/2[x_1..x_3];
  neuralIdeal(C,R)
  ///,
}

--make sure to talk about Factor and Iterative
document{
  Key => {canonicalForm, (canonicalForm,Ideal,Ring),(canonicalForm,Ideal),(canonicalForm,NeuralCode,Ring),(canonicalForm,NeuralCode),[canonicalForm,Factor],[canonicalForm,SharedIndex],[canonicalForm,Iterative]},
  Headline => "Canonical Form",
  TEX "A method which computes the canonical form of a given squarefree pseudomonomial ideal or neural code. If entering a neural code, it is recommended that you also specify the ring where the elements of the canonical form will live. The option Factor returns the canonical form with every element factored. If starting from an Ideal, the default option SharedIndex=>false will compute the canonical form by the method of The Neural Ring. If SharedIndex=>true is specified, the canonical form will be computed using the method of Geller and R.G. instead. If starting from a NeuralCode, the default option Iterative=> true will compute the canonical form of a neural code using the newer method from Neural Ideals in SageMath, but you can use the older method from The Neural Ring by selecting Iterative=>false.",
  Usage => "canonicalForm(Ideal,Ring) or canonicalForm(Ideal) or canonicalForm(Ideal,Ring,SharedIndex=>true) or canonicalForm(Ideal,SharedIndex=>true) or canonicalForm(Ideal,Ring,Factor=>true) or canonicalForm(Ideal,Factor=>true) or canonicalForm(NeuralCode,Ring) or canonicalForm(NeuralCode) or canonicalForm(NeuralCode,Ring,Iterative=>false) or canonicalForm(NeuralCode,Iterative=>false) or canonicalForm(NeuralCode,Ring,Factor=>true) or canonicalForm(NeuralCode,Factor=>true)",
  Inputs => {"Squarefree pseudomonomial ideal or NeuralCode (recommend specifying a ring for the latter)"},
  Outputs => {"The canonical form as a list of elements of the ring of the ideal, the specified ring, or ZZ/2[x_1..x_d] where d is the dimension of the neural code."},
  TEX "We compute an example",
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  I=ideal(x_1*x_3,x_2*(1-x_1));
  canonicalForm(I)
  ///,
    EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  I=ideal(x_1*x_3,x_2*(1-x_1));
  canonicalForm(I,Factored=>true)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  I=ideal(x_1*x_3,x_2*(1-x_1));
  canonicalForm(I,SharedIndex=>true)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  C=neuralCode({"000","001"},R);
  canonicalForm(C)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  C=neuralCode("000","001");
  canonicalForm(C,R,Iterative=>false)
  ///,
}

document{
    Key => {codeSupport, (codeSupport,NeuralCode)},
    Headline => "Support of a NeuralCode",
    TEX "A method which returns a list of the sets of neurons that fire together.",
    Usage => "codeSupport(NeuralCode)",
    Inputs => {"a NeuralCode"},
    Outputs => {"a List of lists of neurons that fire together."},
    TEX "We compute an example",
    EXAMPLE lines ///
    C=neuralCode("000","100","101","001","101");
    codeSupport(C)
    ///
    }

document{
  Key => {neuralIdealToCode, (neuralIdealToCode,List)},
  Headline => "Neural Ideal To Code",
  TEX "A method that computes the neural code corresponding to a list of pseudomonomial generators (generally expected to be in canonical form).",
  Usage => "neuralIdealToCode(List)",
  Inputs => {"List of squarefree pseudomonomials in a single polynomial ring which do not generate the unit ideal"},
  Outputs => {"The corresponding neural code"},
  TEX "We compute some examples",
  EXAMPLE lines ///
  R=ZZ/2[x_1,x_2];
  L={x_1*x_2};
  neuralIdealToCode(L)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  M=ideal(x_1*x_2,x_3*(1-x_1),x_2*x_3);
  neuralIdealToCode(M)
  ///,
}

document{
    Key => {isPseudomonomial,(isPseudomonomial,RingElement)},
    Headline => "isPseudomonomial",
    TEX "A method which determines whether an element of a polynomial ring is a squarefree pseudomonomial. This function was written by Alan Veliz-Cuba for a package on primary decomposition of squarefree pseudomonomial ideals.",
    Usage => "isPseudomonomial(RingElement)",
    Inputs => {"An element of a polynomial ring"},
    Outputs => {"Boolean"},
    TEX "We compute some examples",
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    f=x_1*(1-x_2);
    isPseudomonomial(f)
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    f=1;
    isPseudomonomial(f)
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    f=x_1^2
    ///
    }

document{
    Key => {receptiveFieldRelation,(receptiveFieldRelation,RingElement)},
    Headline => "Receptive field relation corresponding to a pseudomonomial",
    TEX "A method which returns the receptive field relation corresponding to a pseudomonomial.",
    Usage => "receptiveFieldRelation(RingElement)",
    Inputs => {"A squarefree pseudomonomial"},
    Outputs => {"A list containing two lists, such that the intersection of the firing fields corresponding to the elements of the first list is contained in the union of the firing fields corresponding to the elements of the second list."},
    TEX "We compute an example",
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    f=x_1*(1-x_2)*(1-x_3);
    receptiveFieldRelation(f)
    ///
    }

document{
    Key => {polarizePseudomonomial,(polarizePseudomonomial,RingElement,Ring),(polarizePseudomonomial,RingElement)},
    Headline => "Polarize a pseudomonomial",
    TEX "A method which takes a pseudomonomial in a polynomial ring and replaces every instance of (1-var) with a new variable. It is recommended that you specify the ring in which the new monomial will live, see examples below.",
    Usage => "polarizePseudomonomial(RingElement,Ring) or polarizePseudomonomial(RingElement)",
    Inputs => {"A squarefree pseudomonomial, or a squarefree pseudomonomial and the ring the polarization will live in"},
    Outputs => {"A squarefree monomial in a larger polynomial ring."},
    TEX "We compute an example",
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    f=x_1*(1-x_2);
    S=ZZ/2[x_1..x_3,y_1..y_3];
    polarizePseudomonomial(f,S)
    ///
    }

document{
    Key => {polarizedCanonicalForm,(polarizedCanonicalForm,NeuralCode,Ring),(polarizedCanonicalForm,NeuralCode),(polarizedCanonicalForm,Ideal,Ring)},
    Headline => "Polarized canonical form",
    TEX "A method that returns the polarized canonical form of a neural code or neural ideal",
    Usage => "polarizedCanonicalForm(NeuralCode,Ring) or polarizedCanonicalForm(NeuralCode) or polarizedCanonicalForm(Ideal,Ring) or polarizedCanonicalForm(Ideal)",
    Inputs => {"A neural code and a target ring, or just a neural code, or a neural ideal and a target ring."},
    Outputs => {"The polarized canonical form, in the target ring if one is given."},
    TEX "We compute some examples",
    EXAMPLE lines ///
    S=ZZ/2[x_1..x_3,y_1..y_3];
    C=neuralCode("000","100","110","101","001");
    polarizedCanonicalForm(C,S)
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    I=ideal(x_1,x_2*(1-x_1));
    S=ZZ/2[x_1..x_3,y_1..y_3];
    polarizedCanonicalForm(I,S)
    ///
    }

document{
    Key => {polarizedCanonicalIdeal,(polarizedCanonicalIdeal,NeuralCode,Ring),(polarizedCanonicalIdeal,NeuralCode),(polarizedCanonicalIdeal,Ideal,Ring),(polarizedCanonicalIdeal,Ideal)},
    Headline => "Polarized canonical ideal",
    TEX "A method that returns the ideal generated by the polarized canonical form of a neural code or neural ideal.",
    Usage => "polarizedCanonicalIdeal(NeuralCode,Ring) or polarizedCanonicalIdeal(NeuralCode) or polarizedCanonicalIdeal(Ideal,Ring) or polarizedCanonicalIdeal(Ideal)",
    Inputs => {"A neural code or neural ideal, possibly with the ring the polarized canonical form will live in."},
    Outputs => {"The ideal generated by the polarized canonical form, in the target ring if this is specified."},
    TEX "We compute some examples",
    EXAMPLE lines ///
    C=neuralCode("00","10");
    S=ZZ/2[x_1,x_2,y_1,y_2];
    polarizedCanonicalIdeal(C,S)
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1,x_2];
    I=ideal(x_1*x_2);
    S=ZZ/2[x_1,x_2,y_1,y_2];
    polarizedCanonicalIdeal(I,S)
    ///
    }

document{
    Key => {isCanonical,(isCanonical,Ideal,Ring),(isCanonical,Ideal)},
    Headline => "Is Canonical",
    TEX "A method that determines whether a neural ideal is in canonical form. Currently, it computes the canonical form and compares, so it may not be accurate and is not quick.",
    Usage => "isCanonical(Ideal,Ring) or isCanonical(Ideal)",
    Inputs => {"A neural ideal, possibly with its ambient ring"},
    Outputs => {"Boolean"},
    TEX "We compute some examples",
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    I=ideal(x_1,x_2*(1-x_1));
    isCanonical(I,R)
    ///,
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    I=ideal(x_1*x_2,x_1*x_3);
    isCanonical(I)
    ///
    }

document{
    Key => {polarizedCanonicalResolution,(polarizedCanonicalResolution,NeuralCode,Ring),(polarizedCanonicalResolution,NeuralCode)},
    Headline => "Polarized Canonical Resolution",
    TEX "A method that computes the polarized canonical resolution of a neural code, i.e. the minimal free resolution of the ideal generated by its polarized canonical form. When depolarized, this gives the canonical resolution of Gunturkun, Jeffries, and Sun.",
    Usage => "polarizedCanonicalResolution(NeuralCode,Ring) or polarizedCanonicalResolution(NeuralCode)",
    Inputs => {"A neural code and a target ring for the resolution, or just a neural code"},
    Outputs => {"The polarized canonical resolution of the neural code"},
    TEX "We compute an example",
    EXAMPLE lines ///
    C=neuralCode("00","10","01");
    S=ZZ/2[x_1,x_2,y_1,y_2];
    polarizedCanonicalResolution(C,S)
    ///
    }

document{
    Key => {depolarizationMap,(depolarizationMap,Ring,Ring)},
    Headline => "Depolarization Map",
    TEX "A method that returns the depolarization map of Gunturkun, Jeffries, and Sun.",
    Usage => "depolarizationMap(Ring,Ring)",
    Inputs => {"Two rings, first the target polynomial ring and then the source polynomial ring over twice as many variables"},
    Outputs => {"The depolarization map of rings."},
    TEX "We compute an example",
    EXAMPLE lines ///
    R=ZZ/2[x_1..x_3];
    S=ZZ/2[x_1..x_3,y_1..y_3];
    depolarizationMap(R,S)
    ///
    }

document{
    Key => {canonicalResolution,(canonicalResolution,NeuralCode,Ring,Ring),(canonicalResolution,NeuralCode,Ring),(canonicalResolution,NeuralCode)},
    Headline => "Canonical Resolution",
    TEX "A method that returns the canonical resolution of a neural code.",
    Usage => "canonicalResolution(NeuralCode,Ring,Ring) or canonicalResolution(NeuralCode,Ring) or canonicalResolution(NeuralCode)",
    Inputs => {"A neural code, then either both a ring serving as the polarized ring of the neural code and a ring serving as the ring of the neural code (over which the resolution will be computed), or just the ring of the neural code, or neither."},
    Outputs => {"A free resolution over the target ring"},
    TEX "We compute some examples",
    EXAMPLE lines ///
    C=neuralCode("00","10","01");
    S=ZZ/2[x_1..x_3,y_1..y_3];
    R=ZZ/2[x_1..x_3];
    canonicalResolution(C,R,S)
    ///,
    EXAMPLE lines ///
    C=neuralCode("00","10","01");
    R=ZZ/2[x_1..x_3];
    canonicalResolution(C,R)
    ///,
    }



-- **TEST0**
TEST ///
  C=neuralCode("100","010","110","101","011","111");
  I=neuralIdeal(C);
  assert(I == ideal((1-x_1)*(1-x_2)*(1-x_3),(1-x_1)*(1-x_2)*x_3))
///

-- **TEST1**
TEST ///
    C=neuralCode("00","10");
    I=neuralIdeal(C);
    assert(I==ideal((1-x_1)*x_2,x_1*x_2))
///
    
-- **TEST2**
TEST ///
    C=neuralCode("00","10");
    R=ZZ/2[x_1,x_2];
    I=neuralIdeal(C,R);
    cI=canonicalForm(I);
    cC=canonicalForm(C,R);
    cCIter=canonicalForm(C,R,Iterative=>true);
    L={x_2};
    assert((cI==cC) and (cI==L) and (cCIter==cC))
///
    
-- **TEST3**
TEST ///
    R=ZZ/2[x_1,x_2];
    L={x_1*x_2};
    assert(neuralIdealToCode(L)==neuralCode("00","10","01"))
///

--need tests for code support, allCodeWords, isPseudomonomial, sigmaTau, polarizePseudomonomial, polarizedCanonicalResolution, depolarizationMap


end

--***Changelog***---

--1.01, someday


