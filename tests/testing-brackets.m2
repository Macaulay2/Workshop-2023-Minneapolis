needsPackage("Brackets", FileName => "/home/macaulay/Brackets/Brackets.m2")

B = bracketRing(6, 3)
T = [1 4 5]_B * [1 5 6]_B * [2 3 4]_B
n = normalForm T 
assert(net n == "[256]*[145]*[134]-[356]*[145]*[124]+[456]*[145]*[123]")

G = gc(a..f, 3)
A = (a * d)_G
B = (b * e)_G
AB = A ^ B 
C = (c * f)_G
D = AB ^ C -- Output "2*[bde]*[acf]-2*[cdf]*[abe]" is consistent with the book's answer up to sorting and sign.
assert(net D == "2*[bde]*[acf]-2*[cdf]*[abe]")

G = gc(a..f,3)
abLine = (a * b)_G
deLine = (d * e)_G
bcLine = (b * c)_G
efLine = (e * f)_G
acLine = (a * c)_G
dfLine = (d * f)_G
pt1 = abLine ^ deLine
pt2 = bcLine ^ efLine
pt3 = acLine ^ dfLine
linePerspective = pt1 * pt2 * pt3
adLine = (a * d)_G
beLine = (b * e)_G
cfLine = (c * f)_G
pointPerspective =  adLine ^ beLine ^ cfLine
assert(net pointPerspective == "2*[bde]*[acf]-2*[cdf]*[abe]")
(n1, n2) = (normalForm pointPerspective, normalForm linePerspective);
(f1, f2) = (factor n1, factor n2)
assert(net f1#0 == "[bdf]*[ace]-[bef]*[acd]-[cdf]*[abe]-[def]*[abc]")
assert(net f2#2 == "[bdf]*[ace]-[bef]*[acd]-[cdf]*[abe]-[def]*[abc]")
