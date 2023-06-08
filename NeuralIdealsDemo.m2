needsPackage "NeuralIdeals"
R=ZZ/2[x_1..x_3]
C=neuralCode("000","100","110","010","001")
I=neuralIdeal(C,R)
canonicalForm(I,R,Factored=>true)
canonicalForm(C,R,Factored=>true,Iterative=>true)

J=ideal(x_1*(1-x_2),x_2*(1-x_3))
canonicalForm(J,R,Factored=>true)
f=x_1*(1-x_3)
f%J

