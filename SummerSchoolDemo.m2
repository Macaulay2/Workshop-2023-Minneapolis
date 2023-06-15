----------------------------------
----------------------------------
--**DEMO FILE**--
--Schubert Summer School
--16 June 2023
----------------------------------
----------------------------------

restart
uninstallPackage "MatrixSchubert"
restart
installPackage "MatrixSchubert"
restart
needsPackage "MatrixSchubert"
elapsedTime check "MatrixSchubert"
viewHelp "MatrixSchubert"

--------------------------------------------
--Matrix Schubert Varieties & ASM Varieties
--------------------------------------------

--Example from  worksheet
v = {2,7,1,6,3,5,4}
I = schubertDetIdeal v;
netList I_*
netList augmentedRotheDiagram v
essentialBoxes v
betti res trim I
betti res antiDiagInit v


A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}}
isPartialASM A
J = schubertDetIdeal A;
netList J_*


--we can see which Schubert Determintanl Ideals are the prime components of an ASM variety
primePerms = schubertDecomposition J

--we can more effectivey compute the regularity of a matrix Schubert variety using antiDiagInit
--(Knutson-Miller+Knutson or Weigandt+Conca-Varbaro)
--and do even better for a permutation (Pechenik-Speyer-Weigandt)
time regularity comodule I
time matrixSchubertReg v
time regularity comodule (schubertDetIdeal A)
time matrixSchubertReg A

--------------------------------
--operations for permutations
--------------------------------

v
permLength v
schubertPoly v
doubleSchubertPoly v
rajCode v

w
isPatternAvoiding(w,{4,1,2,3})
isCDG w --is Conca-De Negri-Gorla
isVexillary w -- 2143 avoiding
isVexillary v -- not 2143 avoiding
isCartwrightSturmfels w

--Gorenstein example
--takes some time
--u = {3,7,1,4,8,2,6,5}
--betti res trim schubertDetIdeal u
--betti res antiDiagInit u
----------------

--Regularity demo--
--Adam

--Two computations are the same.
for n from 1 to 5 do (
    apply(permutations(toList(1..n)),w->assert(matrixSchubertReg(w)==matrixSchubertRegADI(w)));
);

--Comparing reg with raj implementation and ADI implementation. --

setRandomSeed 50;

apply(1..10,i-> elapsedTime matrixSchubertRegADI(random toList (1..7)));

setRandomSeed 50;

apply(1..10,i-> elapsedTime matrixSchubertReg(random toList (1..7)));


--The speed on larger permutations. --
setRandomSeed 1001;

apply(1..10,i-> elapsedTime matrixSchubertReg(random toList (1..50)));


setRandomSeed 1001;

apply(1..5,i-> elapsedTime matrixSchubertReg(random toList (1..100)));

