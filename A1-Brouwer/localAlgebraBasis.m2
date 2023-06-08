-- load "./GW-type.m2"


-- local A1
-- test



-- Input: L is list of functions {f1,...,fn} over same ring R and p is prime ideal of an isolated zero

-- Output: list of basis elements of local k-algebra Q_p(f) where f = (f1,...,fn):A^n --> A^n

localAlgebraBasis = method()
localAlgebraBasis (List, Ideal) := (List) => (L,p) -> (
    if isPrime(p) == false then (
        error "Error: ideal is not prime"
        );
    R := ring L#0;
    I := ideal(L);
    if dim I > 0  then (
        error "Error: morphism does not have isolated zeroes"
        );
    if isSubset(I,p) == false then (
        error "Error: prime is not a zero of function"
        );
    J := I:saturate(I,p);
    A := R/J;
    B := basis(A);
    return flatten(entries(B))
    )