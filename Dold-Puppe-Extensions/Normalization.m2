altSumFace = method();
altSumFace(SimplicialModule,ZZ) := (S,n) -> (sum(0..n,i->(-1)^i*(S.dd)_(n,i)))

altSumFace(Complex,ZZ) := (C,n) -> (altSumFace(simplicialModule(C,n),n))

naiveNorm = method();
naiveNorm(SimplicialModule,ZZ) := (S,n) -> (complex for i from 1 to n list altSumFace(S,i))

naiveNorm(Complex,ZZ) := (C,n) -> (naiveNorm(simplicialModule(C,n),n))



--we assume that S is Gamma(C) for some complex C here
--sym = method();
--sym(ZZ,Matrix) := (d,M) -> (

sym(ZZ,SimplicialModule) := (d,S) -> (tdeg := topDegree S;
    C := S.complex;
    L := hashTable for i to tdeg list i => symmetricPower(d,combineSFactors(S,i));
    print(L);
    H1 := hashTable for i in keys (S.dd.map) list i => symmetricPower(d,((S.dd)_i));
    print(H1);
    H2 := hashTable for i in keys (S.ss) list i => symmetricPower(d,((S.ss)_i));
    print(H2);
    simplicialModule(L,H1,H2,tdeg)
    )

ext = method()
ext(ZZ,SimplicialModule) := (d,S) -> (tdeg := topDegree S;
    C := S.complex;
    L = hashTable for i to tdeg list i => exteriorPower(d,combineSFactors(S,i));
    H1 = hashTable for i in keys (S.dd.map) list i => map(exteriorPower(d,target (S.dd)_i),exteriorPower(d,source (S.dd)_i),exteriorPower(d,((S.dd)_i)));
    H2 = hashTable for i in keys (S.ss.map) list i => exteriorPower(d,((S.ss)_i));
    print("we made it");
    simplicialModule(L,H1,H2,tdeg)
    )

degenMorphisms = method()
degenMorphisms(SimplicialModule,ZZ,ZZ) := (S,n,k) -> (Cn := naiveNorm(S,n);
    map(Cn,Cn,i->(S.ss)_(i,k),Degree => 1)
    )

degenMorphisms(Complex,ZZ,ZZ) := (C,n,k) -> (degenMorphisms(simplicialModule(C,n),n,k))

faceMorphisms = method()
faceMorphisms(SimplicialModule,ZZ,ZZ) := (S,n,k) -> (Cn := naiveNorm(S,n);
    map(Cn,Cn,i->(S.dd)_(i,k),Degree => -1)
    )

faceMorphisms(Complex,ZZ,ZZ) := (C,n,k) -> (faceMorphisms(simplicialModule(C,n),n,k))

normalize = method();
