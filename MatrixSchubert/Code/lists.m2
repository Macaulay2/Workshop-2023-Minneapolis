------------------------------
-- INPUT: an integer n between 1 and 7
-- OUTPUT: a list of matrices containing all ASMs of size n
------------------------------
ASMFullList = method()
ASMFullList ZZ := List => (n) -> (
    if (n < 1 or n > 7) then error("There is no available list for this n.");
    filename := concatenate("./MatrixSchubert/ASMData/full/", toString n, ".txt");
    listOfMatrices := apply(lines get filename, i -> matrix value i);
    listOfMatrices
);

------------------------------
-- INPUT: an integer n between 1 and 6
-- OUTPUT: a list of matrices containing all CM ASMs of size n
------------------------------
cohenMacaulayASMsList = method()
cohenMacaulayASMsList ZZ := List => (n) -> (
    if (n < 1 or n > 6) then error("There is no available list for this n.");
    filename := concatenate("./MatrixSchubert/ASMData/CM/good", toString n, ".m2");
    listOfMatrices := apply(lines get filename, i -> matrix value i);
    listOfMatrices
);

------------------------------
-- INPUT: an integer n between 1 and 6
-- OUTPUT: a list of matrices containing all non-CM ASMs of size n
------------------------------
nonCohenMacaulayASMsList = method()
nonCohenMacaulayASMsList ZZ := List => (n) -> (
    if (n < 1 or n > 6) then error("There is no available list for this n.");
    filename := concatenate("./MatrixSchubert/ASMData/notCM/bad", toString n, ".m2");
    listOfMatrices := apply(lines get filename, i -> matrix value i);
    listOfMatrices
);

------------------------------
-- INPUT: an integer n between 3 and 6
-- OUTPUT: a list of matrices containing all antidiagonal initial ideals of size n
-- TODO: test and docs
------------------------------
initialIdealsList = method()
initialIdealsList ZZ := List => (n) -> (
    if (n < 3 or n > 6) then error("There is no available list for this n.");
    filename := concatenate("./MatrixSchubert/ASMData/antiDiagIniIdeal/ideals", toString n, ".txt");
    z := getSymbol "z";
    S := QQ(monoid[z_(1,1)..z_(n,n)]);
    listOfIdeals := apply(lines get filename, i -> sub(value i,S));
    listOfIdeals
);