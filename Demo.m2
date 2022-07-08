installPackage "InvariantRing" -- runs all checks
viewHelp "InvariantRing" -- opens documentation in browser

--no invariants example
--SL2 acting on C^2
restart
needsPackage "InvariantRing"
B = QQ[a,b,c,d]
A = ideal(a*d - b*c - 1)
SL2std = matrix{{a,b},{c,d}}
R = QQ[x_1..x_2]
V = linearlyReductiveAction(A,SL2std,R) 
invariants V
elapsedTime hilbertIdeal V

-- abelian group example C3xC3 acting on polynomial ring in two vars
restart
needsPackage "InvariantRing"
R = QQ[x_1..x_3]
W = matrix{{1,0,1},{0,1,1}}
L = {3,3}
T = diagonalAction(W,L,R)
S = R^T
invariantRing T
-- resolution of Hilbert ideal
I = definingIdeal S 
Q = ring I
F = res I
-- Hilbert series of ring of invariants
hilbertSeries S
-- equivariant Hilbert series of polynomial ring
equivariantHilbertSeries T


-- S_2 as a linearly reductive action
restart
needsPackage "InvariantRing"
S = QQ[z]
A = ideal(z^2 - 1)
M = matrix{{(1+z)/2, (1-z)/2},{(1-z)/2,(1+z)/2}}
R = QQ[a,b]
X = linearlyReductiveAction(A,M,R)
isInvariant(a,X)
invariants X

-- invariants of binary quadrics
restart
needsPackage "InvariantRing"
S = QQ[a,b,c,d]
I = ideal(a*d - b*c - 1)
A = S[u,v]
M = transpose (map(S,A)) last coefficients sub(basis(2,A),{u=>a*u+b*v,v=>c*u+d*v})
R = QQ[x_1..x_3]
L = linearlyReductiveAction(I,M,R)
hilbertIdeal L
invariants L
invariants(L,4)
invariants(L,5)


-- invariants of binary quartics
restart
needsPackage "InvariantRing"
S = QQ[a,b,c,d]
I = ideal(a*d - b*c - 1)
A = S[u,v]
M4 = M = transpose (map(S,A)) last coefficients sub(basis(4,A),{u=>a*u+b*v,v=>c*u+d*v})
R4 = QQ[x_1..x_5]
L4 = linearlyReductiveAction(I,M4,R4)
elapsedTime hilbertIdeal L4
elapsedTime X = invariants L4 
-- Other invariants can be obtained from the generators above
g2 = X_0/12
g3 = -X_1/216
256*(g2^3 - 27*g3^2) -- discriminant of the quartic can be expressed in terms of invariants
1728*(g2^3)/(g2^3 - 27*g3^2) -- j-invariant (symmetrized cross-ratio) of quartic also in terms of invariants

-- invariant of 2x2 matrices of binary linear forms with SL_2 action
restart
needsPackage "InvariantRing"
S = QQ[a_(1,1)..a_(2,2),b_(1,1)..b_(2,2),c_(1,1)..c_(2,2)]
I = ideal((det genericMatrix(S,a_(1,1),2,2))-1,
    (det genericMatrix(S,b_(1,1),2,2))-1,
    (det genericMatrix(S,c_(1,1),2,2))-1);
G1 = transpose genericMatrix(S,2,2);
G2 = transpose genericMatrix(S,b_(1,1),2,2);
G3 = transpose genericMatrix(S,c_(1,1),2,2);
R = QQ[x_(1,1,1)..x_(2,2,2)]
L=linearlyReductiveAction(I,G1**G2**G3,R)
elapsedTime inv=invariants L
elapsedTime invariants(L,2)
-- elapsedTime invariants(L,4) takes too long, linear algebra
-- method is very slow with this many variables
-- however hilbertIdeal is fast!
elapsedTime hilbertIdeal(L)
J = hilbertIdeal(L)
isInvariant((J.gens)_(0,0),L)

-- invariants of S_4 using King's algorithm
-- and with the linear algebra method
restart
needsPackage "InvariantRing"
R = QQ[x_1..x_4]
L = apply({[2,1,3,4],[2,3,4,1]},permutationMatrix);
S4 = finiteAction(L,R)
elapsedTime invariants S4
elapsedTime invariants(S4,Strategy=>"LinearAlgebra")
elapsedTime p=primaryInvariants S4
elapsedTime secondaryInvariants(p,S4)
elapsedTime hironakaDecomposition(S4)


-- invariant of 2x2 matrices of ternary linear forms
-- takes a bit of time but computes on Fred's computer
restart
needsPackage "InvariantRing"
S = QQ[a_(1,1)..a_(3,3),b_(1,1)..b_(2,2),c_(1,1)..c_(2,2)]
I = ideal((det genericMatrix(S,a_(1,1),3,3))-1,
    (det genericMatrix(S,b_(1,1),2,2))-1,
    (det genericMatrix(S,c_(1,1),2,2))-1)
G1 = transpose genericMatrix(S,a_(1,1),3,3)
G2 = transpose genericMatrix(S,b_(1,1),2,2)
G3 = transpose genericMatrix(S,c_(1,1),2,2)
R = QQ[x_(1,1,1)..x_(3,2,2)]
L=linearlyReductiveAction(I,G1**G2**G3,R)
elapsedTime H=hilbertIdeal(L,SubringLimit=>1);
needsPackage "Resultants"
A=R[u,v,w]
M=transpose (
    u*promote(genericMatrix(R,x_(1,1,1),2,2),A)+
    v*promote(genericMatrix(R,x_(2,1,1),2,2),A)+
    w*promote(genericMatrix(R,x_(3,1,1),2,2),A)
    )
d=discriminant det M
ideal d==H

-- invariant of 3x3 matrices of binary linear forms
-- this is currently too much for Fred's computer
restart
needsPackage "InvariantRing"
S = QQ[a_(1,1)..a_(2,2),b_(1,1)..b_(3,3),c_(1,1)..c_(3,3)]
I = ideal((det genericMatrix(S,a_(1,1),2,2))-1,
    (det genericMatrix(S,b_(1,1),3,3))-1,
    (det genericMatrix(S,c_(1,1),3,3))-1)
G1 = transpose genericMatrix(S,a_(1,1),2,2)
G2 = transpose genericMatrix(S,b_(1,1),3,3)
G3 = transpose genericMatrix(S,c_(1,1),3,3)
R = QQ[x_(1,1,1)..x_(2,3,3)]
L=linearlyReductiveAction(I,G1**G2**G3,R);
gbTrace=3
elapsedTime H=hilbertIdeal(L,SubringLimit=>1);
-- compare with the construction as discriminant of the
-- determinant of 2x2 generic matrix of linear forms
needsPackage "Resultants"
A=R[u,v]
M=transpose (
    u*promote(genericMatrix(R,x_(1,1,1),2,2),A)+
    v*promote(genericMatrix(R,x_(2,1,1),2,2),A)
    )
d=discriminant det M
ideal d==ideal first inv


-- 2x2 conjugation invariants
restart
needsPackage "InvariantRing"
S = QQ[g_(1,1)..g_(2,2),t]
I = ideal((det genericMatrix(S,2,2))*t-1)
Q = S/I
A = Q[y_(1,1)..y_(2,2)]
Y = transpose genericMatrix(A,2,2)
-- generic group element
g = promote(genericMatrix(S,2,2),A)
-- act by conjugation on a 2x2 generic matrix
-- get corresponding action of 1x4 matrix of variables
G = reshape(A^1,A^4,g*Y*inverse(g)) // (vars A)
G = lift(map(A^4,A^4,G),S)

R = QQ[x_(1,1)..x_(2,2)]
L=linearlyReductiveAction(I,G,R)
elapsedTime H=hilbertIdeal(L)
elapsedTime invariants L


-- 3x3 conjugation invariants
restart
needsPackage "InvariantRing"
S = QQ[g_(1,1)..g_(3,3),t]
I = ideal((det genericMatrix(S,3,3))*t-1)
Q = S/I
A = Q[y_(1,1)..y_(3,3)]
Y = transpose genericMatrix(A,3,3)
-- generic group element
g = promote(genericMatrix(S,3,3),A)
-- act by conjugation on a 2x2 generic matrix
-- get corresponding action of 1x4 matrix of variables
G = reshape(A^1,A^9,g*Y*inverse(g)) // (vars A)
G = lift(map(A^9,A^9,G),S)
R = QQ[x_(1,1)..x_(3,3)]
L=linearlyReductiveAction(I,G,R)
elapsedTime H=hilbertIdeal(L)
elapsedTime invariants(L,1)
elapsedTime invariants(L,2)
elapsedTime invariants(L,3)
