restart
loadPackage("Brackets",FileName=>"../Brackets.m2")
--G = gc(a .. h,4,CoefficientRing => QQ[u,v])
G = gc(a..h,4,CoefficientRing => QQ[u,v],Strategy => Grassmannian)
B = bracketRing G

l1 = (a * b)_G;
l2 = (c * d)_G;
l3 = (e * f)_G;
l4 = (g * h)_G;
t = (u*a)_G + (v*b)_G;

exp1 = (t*(c*d)_G)
exp2 = (e*f)_G
exp3 = (exp1)^(exp2)
l = exp3 * t

q = l*(g*h)_G
qe = q#RingElement
p = flatten entries last coefficients(qe,Variables => {u_(ring qe),v_(ring qe)})
q1 = (p#1^2 - 4*p#0*p#2)_B

normalForm(q1)


