-- File with tests, currently disorganized


-- ----- AKmatrix Tests -------------
    R = ZZ/2[x]; 
    f=1;
    assert(AKmatrix(R,1,f) == sub(matrix{{1,0},{0,1}},R) )


-- -----AKhashTable Tests --------------
    R = ZZ/2[x]; f = 1;
    (rHT,cHT) = AKhashTable(R,1,f);
    assert(keys (rHT#0) == {0});
    assert(keys (rHT#1) == {1});
    assert(rHT#0#0 == 1);
    assert(rHT#1#1 == 1);
    assert(keys (cHT#0) == {0});
    assert(keys (cHT#1) == {1});
    assert(cHT#0#0 == 1);
    assert(cHT#1#1 == 1);

    (rHT,cHT) = swapRows(rHT,cHT,0,1);
    assert(pairs rHT#0 == {(1,1)});
    assert(pairs rHT#1 == {(0,1)});
    assert(pairs cHT#0 == {(1,1)});
    assert(pairs cHT#1 == {(0,1)});

    (rHT,cHT) = swapCols(rHT,cHT,0,1);
    assert(pairs rHT#0 == {(0,1)});
    assert(pairs rHT#1 == {(1,1)});
    assert(pairs cHT#0 == {(0,1)});
    assert(pairs cHT#1 == {(1,1)});
    

    -- note "pairs" has no guaranteed order; need to sort!
    R = ZZ/3[x]; f=1+2*x;
    (rHT,cHT) = AKhashTable(R,1,f);
    assert(sort pairs rHT#0 == sort {(0,1),(2,2*x)});
    assert(sort pairs rHT#1 == sort {(0,2),(1,1)});
    assert(sort pairs rHT#2 == sort {(1,2),(2,1)});
    assert(sort pairs cHT#0 == sort {(0,1),(1,2)});
    assert(sort pairs cHT#1 == sort {(1,1),(2,2)});
    assert(sort pairs cHT#2 == sort {(0,2*x),(2,1)});

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
    assert(fSplittingNumberNonGor(R,1,I) == aNum#0)
    -- even e=2 takes too long...  
    -- assert(fSplittingNumberNonGor(R,2,I) == aNum#1)


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

print("End Tests")