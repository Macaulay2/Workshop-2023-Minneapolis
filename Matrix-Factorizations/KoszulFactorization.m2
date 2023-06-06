KoszulMF = method()
KoszulMF(List) := List -> (

X := ZZdfactorization{matrix{{L#0}}, matrix{{L#1}}};
for i from 2 to #List do X := tensorMF(X, ZZdfactorization{L#i,L#(i+1)});

X
)

---test
KoszulMF {x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100}
