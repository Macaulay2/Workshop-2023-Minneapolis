loadPackage "Complexes"
k =ZZ/101;
-- Koszul complex on one generator
R = k[x];
C = koszulComplex vars R

-- Koszul complex on two generators
R = k[x,y];

C = koszulComplex vars R
