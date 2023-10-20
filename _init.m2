restart
installPackage("A1BrouwerDegrees", RerunExamples=>true)
viewHelp A1BrouwerDegrees


M = matrix(GF(13),{{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = gwClass(M);
    sumDecompositionString(beta)


restart
installPackage("A1BrouwerDegrees")
check A1BrouwerDegrees


alpha = diagonalClass(RR,(1,-1,4,5))
WittIndex(alpha)
sumDecomposition(alpha)
sumDecompositionString(alpha)


beta = diagonalClass(QQ,(2,3))
WittIndex(beta)
sumDecomposition(beta)
sumDecompositionString(beta)


anisotropicDimension(beta)
isotropicDimension(beta)


WittIndex(beta)
anisotropicDimensionQp(beta,3)
anisotropicDimensionQp(beta,2)


alpha = diagonalClass(QQ,(1,-1,4,5))
anisotropicDimension(alpha)
isotropicDimension(alpha)

    n := numRows(alpha.matrix);
n - anisotropicDimension(alpha)

anisotropicDimension(alpha)
isotropicDimension(alpha)


WittIndex(alpha)

anisotropicPart(beta)




alpha2 = diagonalClass(GF(7),(1,-1,4,5))
anisotropicPart(alpha2)


anisotropicDimension(beta)
anisotropicDimensionQp(beta,2)
anisotropicDimensionQp(beta,3)

QQanisotropicDimension2(beta)


viewHelp A1BrouwerDegrees

installPackage("A1BrouwerDegrees", RerunExamples=>true)
check A1BrouwerDegrees

T1 = QQ[z_1..z_2];
f1 = {(z_1-1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = globalA1Degree(f1);
q=ideal {z_1,z_2};
r=ideal {z_1-1,z_2^2-(9/85)};
f1LDq= localA1Degree(f1,q)
f1LDr= localA1Degree(f1,r)
f1LDsum = gwAdd(f1LDq, f1LDr)

alpha = f1GD
beta = f1LDsum

numRows(alpha.matrix) == numRows(beta.matrix)
    
    -- Check if the signatures (Hasse-Witt invariants at RR) agree
signature(alpha) == signature(beta)
    
    -- Check if the discriminants agree
integralDiscriminant(alpha) == integralDiscriminant(beta)
    
    -- Check if all the Hasse-Witt invariants agree
unique(relevantPrimes(alpha) | relevantPrimes(beta))

diagonalForm(alpha)

HilbertSymbol(

diagonalEntries(alpha)

hasseWittInvariant(alpha,2)

isIsomorphicFormQ(f1GD,f1LDsum)

-- HERE"S THE ERROR
HilbertSymbol(-102,15,2)
-- ----

installPackage("A1BrouwerDegrees")
a= -102
b = 15
p = 2
alpha := PadicValuation(a,p)
beta := PadicValuation(b,p)
u := sub(a/p^alpha,ZZ)
v := sub(b/p^beta, ZZ)
d := ((u-1)/2)*((v-1)/2) + alpha*((v^2-1)/8) + beta*((u^2-1)/8)

(-1)^(-154)

(u^2-1)/8q


B1=matrix(QQ, {{1/1, -2/1, 4/1}, {-2/1, 2/1, 0}, {4/1, 0, -7/1}});
B2=matrix(QQ, {{-17198/4225, -166126/975, -71771/1560}, {-166126/975, -27758641/4050, -251077/135}, {-71771/1560, -251077/135, -290407/576}});
assert(isIsomorphicFormQ(B1, B2)===true);
assert(isIsomorphicFormQ(gwClass(B1), gwClass(B2))===true);

B3=matrix(QQ, {{-38/1, -50/1, 23/1}, {-50/1, -62/1, 41/1}, {23/1, 41/1, 29/1}});
assert(isIsomorphicFormQ(B1, B3)===true);
assert(isIsomorphicFormQ(gwClass(B1), gwClass(B3))===true);


relevantPrimes(gwClass(B1))
HasseWittInvariant(gwClass(B1),2)
B2=gwClass(matrix(QQ, {{-17198/4225, -166126/975, -71771/1560}, {-166126/975, -27758641/4050, -251077/135}, {-71771/1560, -251077/135, -290407/576}}))

relevantPrimes(B2)
HasseWittInvariant(B2,2)

L1 = diagonalEntries(integralDiagonalRep(B2))
D2= diagonalClass(QQ,toSequence(diagonalEntries(B2)))

HasseWittInvariant(integralDiagonalRep(B2),2)

L1 = diagonalEntries(B2)

HilbertSymbol(L1_0,L1_2,2)

class L1_0

liftable(-17198/4225,ZZ)


u = 31

u = (u % 2)

u


alpha = diagonalClass(QQ,(2,3))
beta = diagonalClass(QQ,(1,6))

relevantPrimes(alpha)
relevantPrimes(beta)

HasseWittInvariant(alpha,2)
HasseWittInvariant(alpha,3)
HasseWittInvariant(beta,2)
HasseWittInvariant(beta,3)

integralDiscriminant alpha
