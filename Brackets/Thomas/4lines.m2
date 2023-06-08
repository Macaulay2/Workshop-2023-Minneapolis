restart
loadPackage("Brackets",FileName=>"../Brackets.m2")
G = gc(a .. h,4,CoefficientRing => QQ[u,v])
B = bracketRing G
--R = G [u,v] --remove this?

l1 = (a * b)_G; --
l2 = (c * d)_G;
l3 = (e * f)_G;
l4 = (g * h)_G;
t = (u*a)_G + (v*b)_G;

exp1 = (t*(c*d)_G)
exp2 = (e*f)_G
exp3 = (exp1)^(exp2)
l = exp3 * t

quantity1 = l*g*h




