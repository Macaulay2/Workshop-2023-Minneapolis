newPackage(
"NeuralIdeals",
Version => "1.00",
Date => "June 5, 2023",
Authors => {{Name => "Hugh Geller"},{Name => "Rebecca R.G."}},
Headline => "neural ideals",
Keywords => {"Commutative Algebra"},
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
    "factoredCanonicalForm"}

protect codes
protect dimension
protect ambientRing

NeuralCode = new Type of HashTable
NeuralCode.synonym = "neural code"

neuralCode=method(List):= codeList -> (
    d := #(codeList#0);
    x := getSymbol x;
    X:=new NeuralCode from {
	symbol codes => codeList,
	symbol dimension => d,
	symbol ambientRing => ZZ/2[x_1..x_d],
	symbol cache => new CacheTable
	};
    X
    )
 
 dim NeuralCode := C -> C.dimension
 
 ring NeuralCode :=  C -> C.ambientRing


isWellDefined NeuralCode := --Boolean =>
 X -> (
    --check keys
    K:=keys X;
    expectedKeys := set{symbol codes, symbol dimension, symbol ambientRing};
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

neuralIdeal(NeuralCode) := C -> (
    d:=dim C;
    R:=ring C;
    oppC:=neuralCodeOpp(C);
    genList:=for i to #oppC-1 list (
    	prod:=1;
    	for j to d-1 do
	    prod=prod*(1-value((oppC#i)#j)-R_j);
	prod);
    ideal(genList)
    )

canonicalForm = method();

canonicalForm(Ideal) := I -> (
    --this is from our old code, doesn't work yet, may not need pseudomonomial stuff actually
    if isSquarefreePseudomonomialIdeal I==true then (
    decomp := primaryDecomposition I;
    multipliedGens :=product(decomp, i->i);
    R :=ring I;
    d :=numgens R;
    booleanIdeal :=ideal(apply(d,i->(R_i*(1-R_i))));
    booleanR :=R/booleanIdeal;
    reducedGens :=apply(first entries gens multipliedGens,i->sub(i,booleanR));
    almostGens :=unique apply(delete(sub(0,booleanR),reducedGens),i->(sub(i,R)));
    actualGens := for i in almostGens list (
	isDivisible := false;
	for j in almostGens do (
	    if i%j==0 and i =!= j then (isDivisible=true; break));
	if isDivisible==true then continue; i))
    else print "Input must be a squarefree pseudomonomial ideal"
    --delete(,actualGens)
    --to factor:
    --apply(J_*,factor)
    --not working right now
    )

canonicalForm(NeuralCode) := C -> (canonicalForm(neuralIdeal(C)))

factoredCanonicalForm = method()

factoredCanonicalForm(Ideal):= I -> (
    apply(canonicalForm(I),factor)
    )

factoredCanonicalForm(NeuralCode) := C -> (apply(canonicalForm(neuralIdeal(C)),factor))



iterCanonicalForm = method();
iterCanonicalForm(NeuralCode) := C -> (
    d := dim C;
    R := ring C;
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
    canonical
    )

codeSupport = method();
codeSupport(NeuralCode) := C -> (
    fullSupport := {};
    for c in C.codes do (
	cSupport = for i to #c-1 when c#i == 1 list i;
	fullSupport = append(fullSupport, cSupport);
	);
    fullSupport
    )


beginDocumentation()

document{
  Key => NeuralIdeals,
  Headline => "neural ideals",
  EM "NeuralIdeals", " is a package that allows computation of a neural ideal or its canonical form from a neural code",
  Caveat => "In progress"
  }

--document{
--  Key => {neuralCode},
--  Headline => "The type neuralCode",
--  TEX "Turns a list of binary strings of the same length into a NeuralCode type.",
--  Usage => "neuralCode(code)",
--  Inputs => {"Binary strings of the same length"},
--  Outputs => {"The neural code consisting of the given codes."},
--  TEX "We demonstrate how to enter a neural code as a list of binary strings of the same length.",
--  EXAMPLE lines ///
--  neuralCode("000","001","101")
--  ///
--  }

--document{
--  Key => {neuralIdeal(NeuralCode)},
--  Headline => "Neural ideal.",
--  TEX "A method which computes the neural ideal for a given neural code.",
--  Usage => "neuralIdeal(neuralCode(code))",
--  Inputs => {"neuralCode"},
--  Outputs => {"The neural ideal corresponding to the given neural code"},
--  TEX "We compute an example",
--  EXAMPLE lines ///
--  C=neuralCode("000","001");
--  neuralIdeal(C)
--  ///,
--}


-- **TEST0**
TEST ///
  C=neuralCode{"100","010","110","101","011","111"};
  I=neuralIdeal(C);
  assert(I == ideal((1-x_1)*(1-x_2)*(1-x_3),(1-x_1)*(1-x_2)*x_3))
///

-- **TEST1**
TEST ///
    C=neuralCode{"00","10"};
    I=neuralIdeal(C);
    assert(I==ideal((1-x_1)*x_2,x_1*x_2))

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


