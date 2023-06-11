
-- Utility routines --

--------------------------------
--auxiliary function for generating a generic matrix of z variables
--INPUT: integers n,m
--OUTPUT: an n by m generic matrix eith entries z_(i,j)
--NOTE: the ring automatically comes equipped with the antidiagonal term order
--TODO: allow user to input the field they want as an option
-----------------------------------
-*
genMat = (n,m) -> (
    k := QQ;
    zEntries := flatten table (n,m,(i,j) -> (i,j));
    z := getSymbol "z";
    degs := apply(zEntries,i-> i_1-i_0 + m); --are there better ways to make the antidiagonal weights? prob
    Q := k(monoid[z_(1,1)..z_(n,m)]);
    Mmut := mutableMatrix(Q,n,m);
    for box in zEntries do (
        Mmut_(box) = Q_(n*(box_0) + box_1);
    );
    matrix Mmut
)
*-


genMat = method(
    Options => {
	CoefficientRing => QQ,
	Variable => getSymbol "z"
	}
    )
genMat (ZZ,ZZ) := o -> (n,m) -> (
    k := o.CoefficientRing;
    zEntries := flatten table (n,m,(i,j) -> (i,j));
    z := o.Variable;
    degs := apply(zEntries,i-> i_1-i_0 + m); --are there better ways to make the antidiagonal weights? prob
    Q := k(monoid[z_(1,1)..z_(n,m)]);
    Mmut := mutableMatrix(Q,n,m);
    for box in zEntries do (
        Mmut_(box) = Q_(m*(box_0) + box_1);
    );
    matrix Mmut
    )


--------------------------------
--auxiliary function for getting the index of a variable in a ring
--INPUT: an indexed variable
--OUTPUT: the index of the variable
--TODO: add docs
--SUGGESTION: (from Anton) `indexOfVariable = v -> ( i:= index v; last toList R.generatorSymbols#i )`  -- need `debug Core` to use `R.generatorSymbols`
--SUGGESTION: (from Ayah) `(expression(x_1))#1`
--SUGGESTION: (from Mahrud) `last baseName x_(1,2)`
-----------------------------------
indexOfVariable = method()
indexOfVariable RingElement := Sequence => (elem) -> (
    last baseName elem
)
indexOfVariable RingElement := List => (elem) -> (
    last baseName elem
)
indexOfVariable RingElement := ZZ => (elem) -> (
    last baseName elem
)
