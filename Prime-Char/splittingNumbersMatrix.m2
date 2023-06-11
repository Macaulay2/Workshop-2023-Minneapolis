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

AKmatrix = method(TypicalValue=>Matrix)

AKmatrixHelper = (R,e,f) -> (
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

AKmatrix (Ring, ZZ, RingElement) := (R,e,f) -> AKmatrixHelper(R,e,f)

-- why is this necessary? TODO figure out how to avoid
AKmatrix (Ring, ZZ, ZZ) := (R,e,f) -> AKmatrixHelper(R,e,f)

AKmatrix (Ring, ZZ, Ideal) := (R,e,I) -> (
    p:=char R;
    d:=dim R;
    n:=p^e;
    M:=fold((i,j) -> (i|j), for g in I_* list AKmatrixHelper(R,e,g));
    return M
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
        );
        rowList
    );
    (rHT,cHT)
)

swapRows = (rHT,cHT,r,s) -> (
    -- swaps rows r & s
    for jthCol in cHT do (
        if jthCol#?r then (
            if jthCol#?s then (
                temp := jthCol#r;
                jthCol#r = jthCol#s;
                jthCol#s = temp;
            ) else (
                jthCol#s = jthCol#r; remove(jthCol,r)
            )
        ) else (
            if jthCol#?s then (
                jthCol#r = jthCol#s; remove(jthCol,s)
            );
        )
    );
    (switch(r,s,rHT),cHT)
);

swapCols = (rHT,cHT,c,d) -> (
    -- swap cols c & d
    for ithRow in rHT do(
        if ithRow#?c then (
            if ithRow#?d then (
                temp := ithRow#c;
                ithRow#c = ithRow#d;
                ithRow#d = temp; 
            ) else (
                ithRow#d = ithRow#c; 
                remove(ithRow,c);
            )
        ) else (
            if ithRow#?d then (
                ithRow#c = ithRow#d;
                remove(ithRow,d);
            );
        )
    );
    (rHT,  switch(c,d,cHT))
)

hashToMatrix = (R,rHT,cHT) -> (
    map(R^(#cHT),R^(#rHT), 
        flatten for i from 0 to #rHT-1 list (
            for j in keys rHT#i list ((i,j)=>rHT#i#j)
        )
    )
)

blockDiagRecursive = (rHT, cHT, startRow, startCol) -> (
    -- construct a preliminary attempt at a block using 1st row & col
    numRows := startRow; -- # rows in the top left block, including recursed
    numCols := startCol; -- # cols in the top left block, including recursed
    eligibleRows := sort select(keys cHT#startCol, t -> t>=startRow);
    for i in eligibleRows do (
        if i!=numRows then (
            (rHT,cHT) = swapRows(rHT,cHT, numRows,i);
        );
        numRows = numRows + 1;
    );
    eligibleCols := sort select(keys rHT#startRow, t -> t>= startCol);
    for j in eligibleCols do (
        if j!=numCols then (
            (rHT,cHT) = swapCols(rHT,cHT, numCols,j);
        );
        numCols = numCols+1;
    );
    -- check if this swapping already gave a block, or if we need to increase
    needsCheck := true;
    while needsCheck do (
        needsCheck = false;
        for i from startRow to numRows-1 do (
            eligibleCols = sort select(keys rHT#i, 
                            t-> t>=startCol+numCols);
            for j in eligibleCols do (
                if j!=startCol+numCols then (
                    (rHT,cHT) = swapCols(rHT,cHT, numCols,j);
                );
                numCols = numCols+1;
                needsCheck=true;
            );
        );
        for j from startCol to startCol+numCols-1 do (
            eligibleRows = sort select(keys cHT#j,
                            t -> t>= startRow+numRows);
            for i in eligibleRows do (
                if i!=startRow+numRows then (
                    (rHT,cHT) = swapRows(rHT,cHT,numRo)
                )
            )

        )
    );
);

blockDiagonalize = (rHT,cHT) -> (
    fixedRows := 0; -- number of already "fixed" rows 
    fixedCols :=0; -- number of already "fixed" cols 
    topRow := 0;
    topCol := 0;
    doRows := true; -- flag for whether to work on rows or cols
    -- account for not square blocks! may need that on the rectangular matrix
    print(hashToMatrix(R,rHT,cHT));
    while (fixedRows < #rHT and fixedCols < #cHT) do (
        if doRows then (
        topRow = fixedRows;
        eligableRows := sort select(keys cHT#fixedCols, t -> t>=fixedRows);
        print("eligable rows are ");
        print(eligableRows);
        for i in eligableRows do(
            if i==topRow then (
                topRow = i + 1;
            ) else (
                (rHT,cHT) = swapRows(rHT,cHT,topRow,i);
                topRow = topRow+1;
            );
        );
        fixedRows = topRow;
        print("Finished increasing rows; number is");
        print(fixedRows);
        print(hashToMatrix(R,rHT,cHT));
        doRows = not doRows;
        ) else (
        topCol = fixedCols;
        eligableCols := sort select (keys rHT#fixedRows, t -> t>=fixedCols);
        print("eligableCols are");
        print(eligableCols);
        for j in eligableCols do (
            if j == topCol then (
                topCol = j+1;
            ) else (
                (rHT,cHT) = swapCols(rHT,cHT,topCol,i);
                topCol = topCol+1;
            )
        );
        fixedCols = topCol;
        doRows = not doRows;
        print("Finished increasing cols; number is");
        print(fixedCols);
        print(hashToMatrix(R,rHT,cHT));
        )
    );
    (rHT,cHT)
)

---------------------------------------------

fSplittingNumber = method(
    TypicalValue=>ZZ)
    -- Options=>{isGorenstein=>false})
-- TODO unify gor & nongor with an option

fSplittingNumberHelper = (R,e) -> (
    p:=char R;
    n:=p^e;
    J:=systemOfParameters(R);
    I:=ideal(vars(R/J));
    delta:=sub(ann(I),R);
    return degree((frobenius^e(J)):delta^n)
    )
--Assumes R is Gorenstein; computes number for R

fSplittingNumber (Ring, ZZ) := (R,e) -> fSplittingNumberHelper(R,e)
fSplittingNumber (Ring, ZZ, Ideal) := (R,e,I) -> fSplittingNumberHelper(R/I,e)
fSplittingNumber (Ring, ZZ, RingElement) := (R,e,f) -> fSplittingNumberHelper(R/ideal f, e)

fSplittingNumberNonGor = method(TypicalValue=>ZZ)

fSplittingNumberNonGorHelper = (R,e,f) -> (
    p=char R;
    n=p^e;
    S:=R/f; -- will be right if f is ring elem OR an ideal
    M:=coker(sub(AKmatrix(R,e,f),S));
    I:=ideal(vars S);
    phi:=inducedMap(S^1,module(I));
    return numgens source basis(coker(Hom(M,phi)))
)
-- Computes number for R/f

fSplittingNumberNonGor (Ring, ZZ, RingElement) := (R,e,f) -> fSplittingNumberNonGorHelper(R,e,f)

-- again, this feels unnecessary...
fSplittingNumberNonGor (Ring, ZZ, ZZ) := (R,e,f) -> fSplittingNumberNonGorHelper(R,e,f)

fSplittingNumberNonGor (Ring, ZZ, Ideal) := (R,e,J) -> fSplittingNumberNonGorHelper(R,e,J)
