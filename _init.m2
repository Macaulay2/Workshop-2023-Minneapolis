installPackage("A1BrouwerDegrees")
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
