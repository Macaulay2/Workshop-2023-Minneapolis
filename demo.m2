----------------------------------
----------------------------------
--**DEMO FILE**--
--Macaulay2 workshop Minneapolis
--9 June 2023
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

--Example from weekend worksheet
v = {2,1,6,3,5,4}
I = schubertDetIdeal v;
netList I_*
netList augmentedRotheDiagram v
essentialBoxes v
betti res trim I
betti res antiDiagInit v

--speed of regularity function
time regularity comodule I
time matrixSchubertReg v

--Example from weekend worksheet
A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}}
isPartialASM A
J = schubertDetIdeal A;
netList J_*

--we can see which Schubert Determintanl Ideals are the prime components of an ASM variety
primePerms = schubertDecomposition J
primes = apply(primePerms, i -> schubertDetIdeal i)
netList transpose{(intersect(primes_0,sub(primes_1,(vars ring primes_0))))_*,(trim J)_*} --they agree!


--we can more effectivey compute the regularity of a matrix Schubert variety using antiDiagInit
time regularity comodule (schubertDetIdeal A)
time matrixSchubertReg A

--Another example
w = {1,2,3,9,8,4,5,6,7}

--speed test: which way of computing the initial ideal is faster?
time I = schubertDetIdeal w;
time inI = antiDiagInit w;
time ideal leadTerm (schubertDetIdeal w);


--checking they agree
netList transpose {sort((inI)_*), sort((ideal leadTerm I)_*)}

--------------------------------
--operations for permutations
--------------------------------
v
rajCode v
permLength v
schubertPoly v
doubleSchubertPoly v

w
isPatternAvoiding(w,{4,1,2,3})
isCDG w --is Conca-De Negri-Gorla
isVexillary w -- 2143 avoiding
isVexillary v -- not 2143 avoiding
isCartwrightSturmfels w

--Gorenstein example
--takes some time
u = {3,7,1,4,8,2,6,5}
betti res trim schubertDetIdeal u

betti res antiDiagInit u
----------------

--Regularity demo--

--Two computations are the same.--
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
