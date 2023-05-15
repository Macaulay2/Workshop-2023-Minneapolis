-- load necessary packages
needsPackage "Elimination"

-- test code and assertions
R = ZZ/101[a..d]
I = monomialCurveIdeal(R, {1, 3, 4})
assert(eliminate(I, {b}) == ideal(c^4-a*d^3))
