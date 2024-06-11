restart
needsPackage "Brackets"
check(1,"Brackets")
check "Brackets"

restart
needsPackage "Brackets"
G = gc(toList(a..h),4,CoefficientRing=>QQ[l,u])
ell1 = (a*b)_G
ell2 = (c*d)_G
ell3 = (e*f)_G
ell4 = (g*h)_G
pt = (l*a + u * b)_G
ell = ((pt * ell2) ^ ell3) * pt
formula = ell * ell4
(m, c) = coefficients formula
disc = c_(2,0) * c_(0,0) - 4 * c_(1,0)
