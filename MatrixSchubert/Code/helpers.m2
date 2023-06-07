
-- Utility routines --

--------------------------------
--auxiliary function for generating a generic matrix of z variables
--INPUT: integers n,m
--OUTPUT: an n by m generic matrix eith entries z_(i,j)
--NOTE: the ring automatically comes equipped with the antidiagonal term order
--TODO: allow user to input the field they want as an option
-----------------------------------
genMat = (n,m) -> (
    k := QQ;
    zEntries := flatten table (n,m,(i,j) -> (i,j));
    z := local z;
    degs := apply(zEntries,i-> i_1-i_0 + m); --are there better ways to make the antidiagonal weights? prob
    Q:=k[z_(1,1)..z_(n,m),MonomialOrder=>Weights=>degs];
    Mmut:=mutableMatrix(Q,n,m);
    for box in zEntries do (
        Mmut_(box) = z_(box_0+1, box_1+1);
        );
    matrix Mmut
)


-*
findIndices(List, List) := (L, Bs) -> (for ell in L list position(Bs, b -> (b == ell)))
*-

--------------------------------
--auxiliary function for getting the index of a variable in a ring
--INPUT: an indexed variable
--OUTPUT: the index of the variable
--TODO: add docs and tests
-----------------------------------
variableIndex = method()
variableIndex RingElement := Sequence => (elem) -> (
    --convert to string, parse, and select index, convert back
    value replace(".*_+", "", toString elem)
)
-- indexOfVariable = v -> ( i:= index v; last toList R.generatorSymbols#i )  -- need `debug Core` to use `R.generatorSymbols`

