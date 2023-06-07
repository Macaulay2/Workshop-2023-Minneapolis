------------------------------
-- INPUT: an integer n between 1 and 7
-- OUTPUT: a list of matrices containing all ASMs of size n
------------------------------

ASMFullList = method()
ASMFullList ZZ := List => (n) -> (
    if (n < 1 or n > 7) then error("There is no available list for this n.");
    filename := concatenate("./MatrixSchubert/ASMData/full/", toString n, ".m2");
    listOfMatrices := value get filename / matrix;
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
    listOfMatrices := value get filename / matrix;
    listOfMatrices
);