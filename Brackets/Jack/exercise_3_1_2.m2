restart
needsPackage "Brackets"

B = bracketRing(6,3)

M = matrix(B, for i from 4 to 6 list {[1 2 i]_B, [1 i 3]_B, [i 2 3]_B})

normalForm (det M)
