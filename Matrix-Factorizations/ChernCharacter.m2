chernCharacter = method();
chernCharacter(ZZdFactorization) := F -> (L := reverse ((ring F)_*);
    genNo := length L;
    if odd(genNo) then return sub(0,ring F);
    dFeven := product(1..genNo,i->diff(L_(i-1),F.dd_i));
    dFodd := product(1..genNo,i->diff(L_(i-1),F.dd_(i-1)));
    trace(dFeven)-trace(dFodd)
    )


boundaryBulk = method();
boundaryBulk(ZZdFactorization,ZZdFactorizationMap) := (F,phi) -> (

