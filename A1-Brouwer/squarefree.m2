squarefree = method()

squarefree (QQ) := (ZZ) => (n) -> (
    H=hashTable(factor(numerator(n)*denominator(n)));
    product(apply(keys(H),p->p^(H#p%2)))
    );

