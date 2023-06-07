--takes in a rational number or integer and outputs the smallest magnitude integer in its square class
squarefreePart = method()

squarefreePart (QQ) := (ZZ) => (n) -> (
    if n==0 then (
        return 0
        );
    if n >= 0 then (
        H=hashTable(factor(numerator(n)*denominator(n)));
        return product(apply(keys(H),p->p^(H#p%2)))
        );
    if n < 0 then (
        H=hashTable(factor(numerator(-n)*denominator(-n)));
        return -product(apply(keys(H),p->p^(H#p%2)))
        );
    )

squarefreePart (ZZ) := (ZZ) => (n) -> (
    if n==0 then (
        return 0
        );
    if n >= 0 then (
        H=hashTable(factor(n));
        return product(apply(keys(H),p->p^(H#p%2)))
        );
    if n < 0 then (
        H=hashTable(factor(-n));
        return -product(apply(keys(H),p->p^(H#p%2)))
        );
    )

