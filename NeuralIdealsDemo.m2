--Rebecca does place cells demo
needsPackage "NeuralIdeals"
R=ZZ/2[x_1..x_3]
C=neuralCode("000","100","110","010","001")
I=neuralIdeal(C,R)
canonicalForm(I,Factor=>true)
canonicalForm(C,R,Factor=>true,Iterative=>true)

J=ideal(x_1*(1-x_2),x_2*(1-x_3))
canonicalForm(J,Factor=>true)
f=x_1*(1-x_3)
f%J
D=canonicalCode(canonicalForm(J))

--Hugh: add polarization demo

S = ZZ/2[x_1..x_3,y_1..y_3]

polarizedCanonicalResolution(C,R,S)

polarizedCanonicalResolution(C,S)

polarizedCanonicalResolution(C)
