needsPackage "Depth"; -- for the Gorenstein case 
needsPackage "TestIdeals";

expDecomp = (R,p,e,f) -> apply(exponents(f),
      exponent->{coefficient(R_exponent,f)*R_(exponent //p^e),exponent%p^e});
--Gets the exponent vectors of each monomial X^u of the polynomial f, and associates to u the two-element list whose
        --first entry is cX^v and second entry is w, where c is the coefficient of X^u in f and u = p^e*v + w.

simpExpDecomp = (e,f) -> (
    -- returns a MutableHashTable whose keys are (sequences representing) exponents of monomial basis vectors
    R:= ring f; 
    p:= char R; 
    dec := expDecomp(R,p,e,f);
    hashT := new MutableHashTable;
    for t in dec do (ex:=toSequence t#1; coe:=t#0;
        if hashT#?ex then 
        hashT#ex = hashT#ex+coe else hashT#ex = coe);
    hashT
    )        

AKmatrix = (R,e,f) -> (
    time(p:=char R;
    d:=dim R;
    n:=p^e;
    I:= frobenius^e(ideal(vars R));
    monBasis:=(entries basis(R^1/I))#0;
    expBasis:=flatten(apply(monBasis,exponents));
    L:=for i from 0 to p^(e*d)-1 list f*(monBasis#i);
    -- T#i is a hashtable whose keys are exponents & values are coeefs
    --     appearing in f*(ith monomial)
    T:=(for i from 0 to p^(e*d)-1 list simpExpDecomp(e,L#i)););
    return time transpose time matrix time( 
        for t in T list( 
        for exponent in expBasis list (
            if t#?exponent then t#exponent else 0)
	 ))

)

-- returns number which is the sequence as a base p expansion
getIndex = (seq,p,e) -> sum for i from 0 to #seq-1 list (seq#i*p^(e*i));


-- Sparse Version
AKmatrix = (R,e,f) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    I:=frobenius^e(ideal(vars R));
    monBasis:=(entries basis(R^1/I))#0;
    expBasis:=flatten(apply(monBasis,exponents));
    L:=for i from 0 to p^(e*d)-1 list f*(monBasis#i);
    -- OLD VERSION
    -- T#i is a hashtable whose keys are exponents & values are coeefs
    --     appearing in f*(ith monomial)
    T:=flatten (for i from 0 to p^(e*d)-1 list ( expList := simpExpDecomp(e,L#i);
       for key in keys expList list (getIndex(key,p,e), i)=>expList#key));
    map(R^(n^d),R^(n^d),T)
)


-- encode all of our matrices as pairs (rHT, cHT), 
-- where rHT & cHT are lists of hashTables
-- rHT#i#j = cHT#j#i = entry (i,j) [when entry is non-zero]


AKhashTable = (R,e,f) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    I:=frobenius^e(ideal(vars R));
    monBasis:=(entries basis(R^1/I))#0;
    expBasis:=flatten(apply(monBasis,exponents));
    L:=for i from 0 to p^(e*d)-1 list f*(monBasis#i);
    -- T#i is a hashtable whose keys are exponents & values are coeefs
    --     appearing in f*(ith monomial)
    rHT := for i from 0 to p^(e*d)-1 list (
        new MutableHashTable
    );
    cHT := for j from 0 to p^(e*d)-1 list ( 
        rowList := new MutableHashTable;
        expList := simpExpDecomp(e,L#j);
        for key in keys expList do (
            i = getIndex(key,p,e);
            coe = expList#key;
            rowList#i = coe;
            rHT#i#j = coe;
        )
    );
    (rHT,cHT)
)

swapRows = (rHT,cHT,r,s) -> (
    -- swaps rows r & s
    rRow = rHT#r;
    rHT#r = rHT#s;
    rHT#s = rRow;
    for jthCol in cHT do (
        if jthCol#?r then (
            if jthCol#?s then (
                temp = jthCol#r;
                jthCol#r = jthCol#s;
                jthCol#s = temp;
            ) else (
                jthCol#s = jthCol#r; remove(jthCol,r)
            )
        ) else (
            if jthCol#?s then (
                jthCol#r = jthCol#s; remove(jthCol,s)
            ) else ();
        )
    );
    (rHT,cHT)
);

AKmatrixIdeal = (R,e,I) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    M:=fold((i,j) -> (i|j), for g in I_* list AKmatrix(R,e,g));
    return M
    )
    

fSplittingNumber = (R,e) -> (
    p:=char R;
    n:=p^e;
    d:=dim(R);
    J:=systemOfParameters(R);
    I:=ideal(vars(R/J));
    delta:=sub(ann(I),R);
    return degree((frobenius^e(J)):delta^n)
    )
--Assumes R is Gorenstein

fSplittingNumberNonGor = (R,e,f) -> (
    p=char R;
    n=p^e;
    S:=R/ideal(f);
    M:=coker(sub(AKmatrix(R,e,f),S));
    I:=ideal(vars S);
    phi:=inducedMap(S^1,module(I));
    return numgens source basis(coker(Hom(M,phi)))
)

fSplittingNumberNonGorIdeal = (R,e,J) -> (
    p=char R;
    n=p^e;
    S:=R/J;
    M:=coker(sub(AKmatrixIdeal(R,e,J), S));
    I:=ideal(vars S);
    phi:=inducedMap(S^1,module(I));
    return numgens source basis(coker(Hom(M,phi)))
)

-- ----- AKmatrix Tests -------------
    R = ZZ/2[x]; 
    f=1;
    -- assert(AKmatrix(R,1,f) == matrix{{1,0},{0,1}})


-- -----AKmatrxIdeal Tests --------------



-- fSplittingNumber Tests --------------
    R = ZZ/2[x,y,z]; f = z^2-x*y;
    aNum={2, 8, 32};
    assert(fSplittingNumber(R/ideal f,1) == aNum#0)
    assert(fSplittingNumber(R/ideal f,2) == aNum#1)
    assert(fSplittingNumber(R/ideal f,3) == aNum#2)
    assert(fSplittingNumberNonGor(R,1,f) == aNum#0)
    assert(fSplittingNumberNonGor(R,2,f) == aNum#1)
    assert(fSplittingNumberNonGor(R,3,f) == aNum#2)

    R = ZZ/3[x,y,z]; f = z^2-x*y;
    aNum = {5, 41};
    assert(fSplittingNumber(R/ideal f,1) == aNum#0)
    assert(fSplittingNumberNonGor(R,1,f) == aNum#0)
    assert(fSplittingNumber(R/ideal f,2) == aNum#1)
    assert(fSplittingNumberNonGor(R,2,f) == aNum#1)

    R = ZZ/2[a,b,c,d,e,symbol f]; 
    I = ideal(c*e+b*f, c*d+a*f, b*d+a*e);
    aNum = {10, 135, 2010, 31076, 488840, 7755280};
    assert(fSplittingNumberNonGorIdeal(R,1,I) == aNum#0)
    -- even e=2 takes too long...  
    -- assert(fSplittingNumberNonGorIdeal(R,2,I) == aNum#1)


-- Our running example:
needsPackage "InvariantRing";

S = ZZ/2[x,y,z]; f= z^2+x*y;
A = AKmatrix(S,1,f);

EMat = (i,j,n) -> matrix for k from 0 to n-1 list for ell from 0 to n-1 list (
    if (k==i) and (j==ell) then 1 else 0
)
iMat = id_(S^8);
M02 = iMat-z*EMat(0,2,8);
M73 = iMat-z*EMat(7,3,8);
P13=permutationMatrix(8,{[1,3]});
P24=permutationMatrix(8,{[2,4]});
P28=permutationMatrix(8,{[2,8]});
P48=permutationMatrix(8,{[4,8]});
P67=permutationMatrix(8,{[6,7]});
P67*P48*P24*P13*M73*M02*A*M02*M73*P28*P67

