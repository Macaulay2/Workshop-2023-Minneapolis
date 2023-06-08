newPackage(
"NeuralIdeals",
Version => "1.00",
Date => "June 5, 2023",
Authors => {{Name => "Hugh Geller"},{Name => "Rebecca R.G."}},
Headline => "canonical forms of neural ideals",
Keywords => {"Commutative Algebra", "Squarefree monomial ideals"},
DebuggingMode => false,
PackageImports => {"PrimaryDecomposition","PseudomonomialPrimaryDecomposition"},
Reload => false
)

--************************************************************************************************************
--************************************************************************************************************
--***Given a neural code, this package computes the corresponding neural idealand canonical form of the neural ideal.
--************************************************************************************************************
--************************************************************************************************************
--***Acknowledgements:                                                                                     ***
--***Special thanks to Juliette Bruce for contributing the original code for allCodeWords and neuralCodeOpp.***
--************************************************************************************************************
--************************************************************************************************************

export{--types
    "NeuralCode",
    --methods/functions
    "neuralCode",
    "neuralIdeal",
    "canonicalForm",
    "iterCanonicalForm",
    "codeSupport",
    "canonicalCode",
    "isPseudomonomial",
    "sigmaTau",
    "polarizePseudomonomial",
    --Symbols
    "codes"}

protect codes
protect dimension
protect Factor

NeuralCode = new Type of HashTable
NeuralCode.synonym = "neural code"

neuralCode=method(List):= codeList -> (
    d := #(codeList#0);
    X:=new NeuralCode from {
	symbol codes => codeList,
	symbol dimension => d,
	symbol cache => new CacheTable
	};
    X
    )
 
dim NeuralCode := C -> C.dimension


isWellDefined NeuralCode := Boolean => X -> (
    --check keys
    K:=keys X;
    expectedKeys := set{symbol codes, symbol dimension}; 
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
    if not instance(X.codes, List) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a list" <<endl;
	return false
	);
    if X.codes === {} or not all (X.codes, r->instance(r,String)) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a nonempty list of strings" <<endl;
	return false
	);
    if not all (X.codes, r->all(r,i->(value(i)==0 or value(i)==1))) then (
	if debugLevel >0 then
	<< "-- expected 'codes' to be a list of strings of 0's and 1's" << endl;
	return false
	);
    codeList := codes X;
    d:= # (codeList#0);
    if not all (X.codes, r-> #r === d) then (
	if debugLevel > 0 then
	<< "-- expected 'codes' to be a list of equal length strings" << endl;
	return false
	);
    --if codeList == {} then (
	--if debugLevel > 0 then
	--<< "--expected 'codes' to be a nonempty list" <<endl;
--	return false
	--);
    true);

allCodeWords = method();
allCodeWords (ZZ) := (n) ->(
    N := n+1;
    L1 := apply(N,i->(
	    apply(i,i->1)|apply(n-i,j->0)
	    ));
    L2 := unique flatten apply(L1,i->permutations i);
    apply(L2, i-> concatenate(apply(i,j->toString j)))
    )

neuralCodeOpp = method();
neuralCodeOpp  (NeuralCode) := NeuralCode => C ->(
    d := dim C;
    L1 := allCodeWords(d);
    L:=C.codes;
    for i in L do L1=delete(i,L1);
    L1
    )    

neuralIdeal = method();

neuralIdeal(NeuralCode,Ring) := Ideal => (C,R) -> (
    d:=dim C;
    if numgens R =!= d then error "Expected ring of the same dimension as the neuralCode";
    if coefficientRing R =!= ZZ/2 then error "Expected coefficientRing of ring to be ZZ/2";
    if instance(R,PolynomialRing)==false then error "Expected ring to be a PolynomialRing";
    oppC:=neuralCodeOpp(C);
    genList:=for i to #oppC-1 list (
    	prod:=1;
    	for j to d-1 do
	    prod=prod*(1-value((oppC#i)#j)-R_j);
	prod);
    ideal(genList)
    )

neuralIdeal NeuralCode := Ideal => C -> (
    d:=dim C;
    x:=getSymbol "x"; 
    R:=(ZZ/2)(monoid[x_1..x_d]); 
    neuralIdeal(C,R)
    )

canonicalForm = method(Options => {Factor => false});

canonicalForm Ideal := List => opts -> I -> (
    --if isSquarefreePseudomonomialIdeal I==true then (
    decomp := primaryDecompositionPseudomonomial I;
    multipliedGens :=product(decomp, i->i);
    R := ring I;
    d :=numgens R;
    booleanIdeal :=ideal(apply(d,i->(R_i*(1-R_i))));
    booleanR :=R/booleanIdeal;
    reducedGens :=apply(first entries gens multipliedGens,i->sub(i,booleanR));
    noDuplicateGens :=delete(sub(0,booleanR),reducedGens);
    almostGens :=unique apply(noDuplicateGens,i->(sub(i,R)));
    actualGens := for i in almostGens list (
	isDivisible := false;
	for j in almostGens do (
	    if i%j==0 and i =!= j then (isDivisible=true; break));
	if isDivisible==true then continue; 
	if opts.Factor == true then factor(i) else i
	)
    --)
    --else error "Input must be a squarefree pseudomonomial ideal"
    )

canonicalForm NeuralCode := List => opts -> C -> (
    canonicalForm(neuralIdeal(C),Factor => opts.Factor)
    )

canonicalForm(NeuralCode,Ring) := List => opts -> (C,R) -> (
    canonicalForm(neuralIdeal(C,R),Factor => opts.Factor)
    )



iterCanonicalForm = method();
iterCanonicalForm(NeuralCode,Ring) := List => (C,R) -> (
    d := dim C;
    if numgens R =!= d then error "Expected ring of the same dimension as the neuralCode";
    initCode := C.codes#0;
    canonical := {};
    for i to d-1 do (
	canonical = append(canonical,R_i - value(initCode#i))
	);
    for i from 1 to #C.codes - 1 do (
	current := C.codes#i;
	codeCoordinate := {};
	factors := {};
	subs := {};
	for j to #current - 1 do (
	    c :=  value(current#j);
	    codeCoordinate = append(codeCoordinate,c);
	    factors = append(factors,R_j  - c);
	    subs = append(subs,R_j => c);
	    );
	currentGens := canonical;
	M := {};
	N := {};
	L := {};
	for gen in currentGens do (
	    if sub(gen,subs) == 0
	    then M = append(M,gen)
	    else N = append(N,gen);
	    );
	for ngen in N do (
	    for fac in factors do (
		goToNext := false;
		g := ngen*fac;
		if ngen%(fac - 1) == 0 then continue;
		for mgen in M do (
		    if g%mgen == 0 then goToNext = true;
		    break;
		    );
		if goToNext then continue;
		L = append(L,g);
		);
	    );
	canonical = join(M,L);
	);
--    C.cache#iCF = canonical;
    canonical
    )

iterCanonicalForm(NeuralCode) := List => C -> (
    d:=dim C;
    x:=getSymbol "x"; 
    R:=(ZZ/2)(monoid[x_1..x_d]); 
    iterCanonicalForm(C,R)
    )

codeSupport = method();
codeSupport(NeuralCode) := C -> (
    fullSupport := {};
    L := C.codes;
    for c in L do (
	cSupport := for i to #c-1 list (if value(c#i) == 0 then continue; i+1);
	fullSupport = append(fullSupport, cSupport);
	);
    fullSupport
    )

canonicalCode = method();

canonicalCode List := NeuralCode => L -> (
    --want to check that elements in L are in same ring, return error otherwise
    R := ring L#0;
    d := numgens R;
    allCodeList := allCodeWords(d);
    codeList := for i in allCodeList list (
	validCode := true;
	for j in L do (
	    if sub(j,matrix{apply(d,k->sub(value(i#k),R))})=!= 0 then (
		validCode = false;
		break
		);
	    );
	if validCode == false then continue else i
	);
    neuralCode(codeList) --getting error array index 0 out of bounds 0 .. -1 with current state
    )

----The following function is an internal function from the PseudomonomialPrimaryDecomposition package by Alan Veliz-Cuba

-- determines if a polynomial is square free pseudomonomial
-- Input:
-- Polynomial P in bitwise form
-- Output:
-- true or false
isPseudomonomial = method();
isPseudomonomial(RingElement) : = P -> ( 
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

sigmaTau = method();

sigmaTau(RingElement) := P -> (
    if isPseudomonomial(P) == false then error "Expected input to be a Pseudomonomial";
    R := ring P;
    d := numgens R;
    sigma := {};
    tau := {};
    )
	

polarizePseudomonomial = method();

polarizePseudomonomial(RingElement) := P -> (
    if isPseudomonomial(P) == false then error "Expected input to be a Pseudomonomial";
    )


beginDocumentation()

document{
  Key => NeuralIdeals,
  Headline => "neural ideals",
  EM "NeuralIdeals", " is a package that allows computation of a neural ideal or its canonical form from a neural code",
  Caveat => "In progress"
  }

document{
  Key => {neuralCode},
  Headline => "The type neuralCode",
  TEX "Turns a list of binary strings of the same length into a NeuralCode type.",
  Usage => "neuralCode(code)",
  Inputs => {"Binary strings of the same length like 000"},
  Outputs => {"The neural code consisting of the given codes."},
  TEX "We demonstrate how to enter a neural code as a list of binary strings of the same length.",
  EXAMPLE lines ///
  neuralCode("000","001","101")
  ///
  }

document{
  Key => {neuralIdeal, (neuralIdeal,NeuralCode)},
  Headline => "Neural ideal.",
  TEX "A method which computes the neural ideal for a given neural code.",
  Usage => "neuralIdeal(neuralCode(code))",
  Inputs => {"neuralCode"},
  Outputs => {"The neural ideal corresponding to the given neural code"},
  TEX "We compute an example",
  EXAMPLE lines ///
  C=neuralCode("000","001");
  neuralIdeal(C)
  ///,
}

document{
  Key => {canonicalForm, (canonicalForm,Ideal),(canonicalForm,NeuralCode)},
  Headline => "Canonical Form",
  TEX "A method which computes the canonical form of a given squarefree pseudomonomial ideal or neural code.",
  Usage => "canonicalForm(Ideal) or canonicalForm(NeuralCode)",
  Inputs => {"Squarefree pseudomonomial ideal or NeuralCode"},
  Outputs => {"The canonical form"},
  TEX "We compute an example",
  EXAMPLE lines ///
  C=neuralCode("000","001");
  canonicalForm(C)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  I=ideal(x_1*x_3,x_2*(1-x_1));
  canonicalForm(I)
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
    
-- **TEST2**
TEST ///
    C=neuralCode("00","10");
    R=ZZ/2[x_1,x_2];
    I=neuralIdeal(C,R);
    cI=canonicalForm(I);
    cC=canonicalForm(C,R);
    L={x_2};
    assert((cI==cC) and (cI==L))

-- **TEST1**  This makes a pinch point.  We check that it has one minimal prime, that it has 3 variables, and that the singular locus is dimension 1 while the ambient object is dimension 2.  We also check that the ring we construct is a subring of A.
--TEST ///
--  A = QQ[x,y];
--  I = ideal(x);
--  B = A/I;
--  C = QQ[u];
--  f = map(B, A);
--  g = map(B, C, {y^2}); 
--  l1 = pullback(f,g);
--  vlist = first entries vars (l1#0);
--  assert ( (#(vlist) == 3) and (dim l1#0 == 2) and ((#minimalPrimes (ideal l1#0)) == 1) and (1 == dim singularLocus (l1#0)) and (ker (l1#1) == ideal(sub(0,l1#0))) )
--///
  


end

--***Changelog***---

--1.01, someday


