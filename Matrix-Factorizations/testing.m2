restart
needs "functionsMF_new.m2"

S = QQ[x,y,z]

--Make a factorization using just ring elements
X = ZZdfactorization{x,y,y^2}
isdFactorization X
period X
X.period

--Make a factorization using matrices
Y = ZZdfactorization{matrix{{x,y}, {-y,x}}, matrix{{x,-y}, {y, x}}}
isdFactorization Y

--Fake factorization
Y' = ZZdfactorization{matrix{{x,y}, {y,y}}, matrix{{x,-y}, {y, x}}}
isdFactorization Y'

--Things that work for d=2 only
shift Y
isdFactorization shift Y --factors same polynomial


--Hom
X' = ZZdfactorization{x, y, z}
Q = adjoinRoot(3, ring X')
H = Hom(X'**Q, X'**Q, t) --Note: base change to Q necessary!!
isdFactorization H


--Koszul factorizations
S = QQ[x,y,z,w]
Q = adjoinRoot(5, S)
K = koszulMFf({x,y,z,w}, x^5+y^5+w^5, 4, t)
dd^K
isdFactorization K


--
--messing around with cycles in Hom
Z = ZZdfactorization{x,y}
H = Hom(Z**Q, Z**Q)
isdFactorization H
C = chainComplex(H.dd_0, H.dd_1)
C.dd
ker(C.dd_1)
ker(C.dd_2)
C.dd_2
