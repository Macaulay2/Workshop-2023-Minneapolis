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

ext = method(Options => {Degeneracy=>false,TopDegree => null})
ext(ZZ,SimplicialModule) := SimplicialModule => opts -> (d,S) -> (tdeg := topDegree S;
    --C := S.complex;
    L = hashTable for i to tdeg list i => exteriorPower(d,combineSFactors(S,i));
    H1 = hashTable for i in keys (S.dd.map) list i => map(exteriorPower(d,target (S.dd)_i),exteriorPower(d,source (S.dd)_i),exteriorPower(d,((S.dd)_i)));
    if opts.Degeneracy==true then H2 = hashTable for i in keys (S.ss.map) list i => exteriorPower(d,((S.ss)_i));
    --print("we made it");
    if opts.Degeneracy==true then return simplicialModule(L,H1,H2,tdeg);
    simplicialModule(L,H1,tdeg)
    )

ext(ZZ,SimplicialModuleMap) := SimplicialModuleMap => opts -> (d,phi) -> (
    S1 := source phi;
    S2 := target phi;
    if instance((keys phi.map)#0,Sequence) then error "Expected SimplicialModuleMap to have singly graded indices";
    map(ext(d,S2),ext(d,S1),new HashTable from for i to max(topDegree S1,topDegree S2) list i => exteriorPower(d,phi_i),Degree => degree phi)
    )

ext(ZZ,Complex) := Complex => opts -> (d,C) -> (if not(opts.TopDegree === null) then S = simplicialModule(C,opts.TopDegree)
    else S = simplicialModule(C,d*length(C));
    Sn := ext(d,S);
    normalize(Sn)
    )

ext(ZZ,ComplexMap) := ComplexMap => opts -> (d,phi) -> (if not(opts.TopDegree === null) then phin = simplicialModuleMap(phi,opts.TopDegree)
    else phin = simplicialModuleMap(phi,d*(max(length(source phi),length target phi)));
    phik := ext(d,phin);
    normalize(phik)
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

makeNormMap = method();
makeNormMap(SimplicialModule,ZZ) := (S,d) -> (
    K2 := intersect(for i from 1 to d list ker S.dd_(d,i));
    if d>1 then K1 = intersect(for i from 1 to d-1 list ker S.dd_(d-1,i));
    if d ==1 then K1 = S_0;
    prune inducedMap(K1,K2,S.dd_(d,0))
    )

makeNormMap(SimplicialModuleMap,ZZ) := (phi,d) -> (
    S1 := source phi;
    S2 := target phi;
    if d==0 then return phi_0;
    K1 := intersect( for i from 1 to d list ker S1.dd_(d,i));
    K2 := intersect( for i from 1 to d list ker S2.dd_(d,i));
    prune inducedMap(K2,K1,phi_d)
    ) 

normalize = method();
normalize(SimplicialModule,ZZ) := (S,d) -> (complex for i from 1 to d list makeNormMap(S,i))

normalize(SimplicialModule) := S -> (if any(keys S,i->i==symbol complex) then return normalize(S,S.complexLength);
    normalize(S,S.topDegree)
    )

normalize(SimplicialModuleMap,ZZ) := (phi,d) -> (C1 := normalize(source phi,d);
    C2 := normalize(target phi,d);
    map(C2,C1,new HashTable from for i to max(length C1,length C2) list i => makeNormMap(phi,i))
    )

normalize(SimplicialModuleMap) := phi -> (d := max(topDegree source phi,topDegree target phi);
    normalize(phi,d)
    )

