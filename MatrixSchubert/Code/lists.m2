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
