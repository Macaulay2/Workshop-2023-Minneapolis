newPackage(
"NeuralIdeals",
Version => "1.00",
Date => "June 5, 2023",
Authors => {{Name => "Hugh Geller"},{Name => "Rebecca R.G."}},
Headline => "canonical forms of neural ideals",
Keywords => {"Commutative Algebra", "Squarefree monomial ideals"},
DebuggingMode => true,
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
    "codeSupport",
    "canonicalCode",
    "isPseudomonomial",
    "sigmaTau",
    "polarizePseudomonomial",
    "allCodeWords",
    "polarizedCanonicalRes",
    --Symbols
    "codes",
    "Factor",
    "Iterative"}

protect codes
protect dimension
protect Factor
protect Iterative

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

ring NeuralCode := C -> (
    d := dim C;
    x := getSymbol "x";
    R := (ZZ/2)(monoid[x_1..x_d])
    )


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

neuralCodeComplement = method();
neuralCodeComplement  (NeuralCode) := NeuralCode => C ->(
    d := dim C;
    L1 := allCodeWords(d);
    L:=C.codes;
    for i in L do L1=delete(i,L1);
    L1
    )    

--input: a NeuralCode, outputs a neural ideal by the method of The Neural Ring, not necessarily in canonical form
neuralIdeal = method();

neuralIdeal(NeuralCode,Ring) := Ideal => (C,R) -> (
    d:=dim C;
    if numgens R =!= d then error "Expected ring of the same dimension as the neuralCode";
    if coefficientRing R =!= ZZ/2 then error "Expected coefficientRing of ring to be ZZ/2";
    if instance(R,PolynomialRing)==false then error "Expected ring to be a PolynomialRing";
    oppC:=neuralCodeComplement(C);
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

--not exported, used as an option for canonicalForm below it, iterative method from NeuralIdeals in SageMath paper
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


--inputs a squarefree pseudomonomial ideal or a NeuralCode, outputs the canonical form of its neural ideal
canonicalForm = method(
    Options => {
	Factor => false,
	Iterative => false
	});

canonicalForm Ideal := List => opts -> I -> (
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
    )

canonicalForm(NeuralCode,Ring) := List => opts -> (C,R) -> (
    if opts.Iterative == true then (
	if opts.Factor == true then apply(iterCanonicalForm(C),factor) else
	iterCanonicalForm(C,R)
	)
    else
    canonicalForm(neuralIdeal(C,R),Factor => opts.Factor)
    )

canonicalForm NeuralCode := List => opts -> C -> (
    if opts.Iterative == true then (
	if opts.Factor == true then apply(iterCanonicalForm(C),factor) else
	iterCanonicalForm(C)
	)
    else 
    canonicalForm(neuralIdeal(C),Factor => opts.Factor)
    )

--finds the support of a given neural code
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

--given a non-unit squarefree pseudomonomial ideal, preferably in canonical form, and outputs the corresponding NeuralCode
canonicalCode = method();

canonicalCode List := NeuralCode => L -> (
    --want to check that elements in L are in same ring, return error otherwise
    R := ring L#0;
    d := numgens R;
    if isSquarefreePseudomonomialIdeal(ideal(L))==false then error "Expected elements that generate a squarefree pseudomonomial ideal.";
    if ideal(L)==sub(ideal(1),R) then error "Expected generators of a non-unit ideal.";
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
	if validCode == false then continue else i
	);
    neuralCode(codeList)
    )

----The following function is an internal function from the PseudomonomialPrimaryDecomposition package by Alan Veliz-Cuba

-- determines if a polynomial is square free pseudomonomial
-- Input:
-- Polynomial P in bitwise form
-- Output:
-- true or false
isPseudomonomial = method();
isPseudomonomial(RingElement) := P -> ( 
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
    for i to d-1 do (
	if P%R_i == 0 then sigma = append(sigma,i+1);
	if P%(1-R_i) == 0 then tau = append(tau,i+1);
	);
    {sigma,tau}
    )
	

polarizePseudomonomial = method();

polarizePseudomonomial(RingElement,Ring,Ring) := RingElement => (P,R,S) -> (  -----ring element, the ring you want the element to currently live in, the ring you want the polarization to live in----
    if isPseudomonomial(P) == false then error "Expected input to be a Pseudomonomial";
    if (numgens S)%2 != 0 then error "Second ring must have an even number of generators";
    d := (numgens S)//2;
    st := sigmaTau(P);
    sigma := st_0;
    tau := st_1;
    use S;
    mon := 1_S;
    for i in sigma do mon = mon*S_(i-1);
    for i in tau do mon = mon*S_(d+i-1);
    mon
    )

polarizePseudomonomial(RingElement,Ring) := RingElement => (P,R) -> (
     d := numgens R;
     y := getSymbol "y";
     x := getSymbol "x";
     S := (ZZ/2)(monoid[x_1..x_d,y_1..y_d]);
     polarizePseudomonomial(P,R,S)
     )

polarizePseudomonomial(RingElement) := RingElement => P -> (
    R := ring P;
    polarizePseudomonomial(P,R)
    )

polarizedCanonicalRes = method();

polarizedCanonicalRes(NeuralCode,Ring,Ring) := Resolution => (C,R,S) -> (
    L := canonicalForm(C,R);
    polarL := for P in L list polarizePseudomonomial(P,R,S);
    res ideal(polarL)
    )

polarizedCanonicalRes(NeuralCode,Ring) := Resolution => (C,S) -> (
    if (numgens S)%2 != 0 then error "Second ring must have an even number of generators";
    d := (numgens S)//2;
    x := getSymbol "x";
    R := (ZZ/2)(monoid[x_1..x_d]);
    polarizedCanonicalRes(C,R,S)
    )

polarizedCanonicalRes(NeuralCode) := Resolution => C -> (
    d := dim C;
    x := getSymbol "x";
    y := getSymbol "y";
    S := (ZZ/2)(monoid[x_1..x_d,y_1..y_d]);
    polarizedCanonicalRes(C,S)
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
  Key => {canonicalForm, (canonicalForm,Ideal),(canonicalForm,NeuralCode,Ring),(canonicalForm,NeuralCode),[canonicalForm,Factor],[canonicalForm,Iterative]},
  Headline => "Canonical Form",
  TEX "A method which computes the canonical form of a given squarefree pseudomonomial ideal or neural code. If entering a neural code, you can also specify the ring where the elements of the canonical form will live. The option factor factors the elements of the list. The option Iterative=> true will compute the canonical form of a neural code using the newer method from NeuralIdeals in SageMath.",
  Usage => "canonicalForm(Ideal) or canonicalForm(NeuralCode,Ring) or canonicalForm(NeuralCode) or canonicalForm(NeuralCode,Ring,Iterative=true) or canonicalForm(NeuralCode,Iterative=true)",
  Inputs => {"Squarefree pseudomonomial ideal or NeuralCode (recommend specifying a ring for the latter), Iterative or not (if entering a neural code)"},
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
  C=neuralCode({"000","001"},R);
  canonicalForm(C)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  C=neuralCode({"000","001"},R,Iterative=>true)
  ///,
}

document{
  Key => {canonicalCode, (canonicalCode,List)},
  Headline => "Canonical Code",
  TEX "A method which computes the neural code corresponding to a list of pseudomonomial generators (generally expected to be in canonical form).",
  Usage => "canonicalCode(List)",
  Inputs => {"List of squarefree pseudomonomials in a single polynomial ring which do not generate the unit ideal"},
  Outputs => {"The corresponding neural code"},
  TEX "We compute an example",
  EXAMPLE lines ///
  R=ZZ/2[x_1,x_2];
  L={x_1*x_2};
  canonicalCode(L)
  ///,
  EXAMPLE lines ///
  R=ZZ/2[x_1..x_3];
  M=ideal(x_1*x_2,x_3*(1-x_1),x_2*x_3);
  canonicalCode(M)
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
    assert(canonicalCode(L)==neuralCode("00","10","01"))
///


end

--***Changelog***---

--1.01, someday

