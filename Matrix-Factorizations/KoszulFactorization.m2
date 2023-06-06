KoszulMF = method()
KoszulMF(List) := L -> (

X := ZZdfactorization{map(R^1, R^1, matrix{{L#0}}), map(R^1, R^1, matrix{{L#1}})};
for i from 1 to (#L)//2-1 do X = tensorMF(X, ZZdfactorization{ map(R^1, R^1, L#(2*i)), map(R^1, R^1, L#(2*i+1))});

X
)

---test
K = KoszulMF({x^2*y^2,z^2*w^2, x*y, y^10, x*y, z^100})

KMF = KoszulMF({x,x,y,y,z,z})

(dd^KMF)_1*(dd^KMF)_2
