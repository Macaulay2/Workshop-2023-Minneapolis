 needsPackage "Complexes"
 
 SimplicialModule = new Type of MutableHashTable --
  -- note: we make this mutable in order to construct the
  --   differential as a morphism of ZZdFactorizations (in the style of Complexes)
  -- BUT: after construction, it is should be IMMUTABLE!!
  -- at some point, we might want to allow lazy determination of the modules and maps
  -- but for now, we insist that all modules and maps are explicit.
  -- key:
  --  ring
  --  modules: hash table: ZZ => Module
  --  differential: ComplexMap from C --> C, degree -1 (mod d)
  --  period: an integer d, tells us for which d the factorization C is ZZ/d-graded
  --  cache: a CacheTable
  
 SimplicialModuleMap = new Type of HashTable
  -- keys:
  --   degree: ZZ
  --   source: Simplicial module over a ring R
  --   target: simplicial module over the same ring R
  --   maps themselves (HashTable of Matrices), keys lying in the concentration period of the source.
  --    not all of the keys maps#i, need be present.
  --    missing ones are presumed to be zero maps.
  --   cache: a CacheTable
  --    cache.isCommutative: whether this map commutes with the differentials
  --      not set until needed.  unset means we have not checked yet, 
  --          and the user hasn't declared it to be true/false yet.

SimplicialModule.synonym = "Simplicial Module"
SimplicialModuleMap.synonym = "Map of Simplicial Modules"

topDegree = method();
topDegree SimplicialModule := ZZ => S -> S.topDegree    
   
ring SimplicialModule := Ring => S -> S.ring

moduleMaker = (C,d) -> (
    moduleList := new MutableHashTable;
    	maxK := min (d, length C);
	for k to maxK do (
	moduleList#(d,k) = directSum toList(
	    binomial(d,k):(C_k));
	    );    
    for i in (sort keys moduleList) list (i,moduleList#i)
    )

mapMaker = (phi,d) -> (
    mapList := new MutableHashTable;
    	maxK := min (d, max(length source phi,length target phi));
	for k to maxK do (
	mapList#(d,k) = directSum toList(
	    binomial(d,k):(phi_k));
	    );  
    for i in (sort keys mapList) list (i,mapList#i)
    )

combineSFactors = method(Options => {RememberSummands => true});
combineSFactors(SimplicialModule,ZZ) := Module => opts -> (S,d) -> (modwComps := directSum for j in components S list directSum for i to min(d,j.complexLength) list (j.module)#(d,i);
    if opts.RememberSummands then modwComps.cache.components = flatten flatten for j in components S list for i to min(d,S.complexLength) list components (S.module#(d,i));
    modwComps
    )
   
   
--WARNING: it seems like there are issues when the complex is 0 in homological degree 0    
--H1 is the face maps, H2 is the degeneracy maps
simplicialModule = method(Options => {Base=>0,Degeneracy => false})
simplicialModule(Complex,HashTable,HashTable,ZZ) := SimplicialModule => opts -> (C,H1,H2,d) -> (
    spots := sort keys H1;
    if #spots === 0 then
      error "expected at least one map";
    R := ring C;
    moduleList := new MutableHashTable;
    for b to d do (
    	maxK := min (b, length C);
	for k to maxK do (
	    modwComps := directSum toList(binomial(b,k):(C_k));
	    modwComps.cache.components = flatten toList(binomial(b,k):(components (C_k)));
	    moduleList#(b,k) = modwComps
	);
	);
    S := new SimplicialModule from {
	symbol ring => R,
	symbol topDegree => d,
	symbol module => new HashTable from moduleList,
	symbol cache => new CacheTable,
    	symbol complexLength => length C,
	symbol complex => C
	};
    S.dd = map(S,S,H1,Degree=>-1);
    S.ss = map(S,S,H2,Degree=>1);
    S
    )

--H1 is the face maps
simplicialModule(Complex,HashTable,ZZ) := SimplicialModule => opts -> (C,H1,d) -> (--print("made it here!");
    spots := sort keys H1;
    if #spots === 0 then
      error "expected at least one map";
    R := ring C;
    moduleList := new MutableHashTable;
    for b to d do (
    	maxK = min (b, length C);
	for k to maxK do (
	modwComps := directSum toList(binomial(b,k):(C_k));
	    modwComps.cache.components = flatten toList(binomial(b,k):(components (C_k)));
	    moduleList#(b,k) = modwComps
	);
	);
    S := new SimplicialModule from {
	symbol ring => R,
	symbol topDegree => d,
	symbol module => new HashTable from moduleList,
	symbol cache => new CacheTable,
    	symbol complexLength => length C,
	symbol complex => C
	};
    S.dd = map(S,S,H1,Degree=>-1);
    S
    )

simplicialModule(HashTable,HashTable,HashTable,ZZ) := SimplicialModule => opts -> (L,H1,H2,d) -> (--print("we got started");
    R := ring (L#((keys L)#0));
    S := new SimplicialModule from {
	symbol ring => R,
	symbol topDegree => d,
	symbol module => L,
	symbol cache => new CacheTable,
	};
    S.dd = map(S,S,H1,Degree=>-1);
    S.ss = map(S,S,H2,Degree=>1);
    S
    )

simplicialModule(HashTable,HashTable,ZZ) := SimplicialModule => opts -> (L,H1,d) -> (--print("we got started");
    R := ring (L#((keys L)#0));
    S := new SimplicialModule from {
	symbol ring => R,
	symbol topDegree => d,
	symbol module => L,
	symbol cache => new CacheTable,
	};
    S.dd = map(S,S,H1,Degree=>-1);
    S
    )


simplicialModule(Complex,ZZ) := SimplicialModule => opts -> (C,d) -> (
    --C is a chain complex, output is the Dold-Kan image of C in the category of simplicial modules
     if not instance(opts.Base, ZZ) then
      error "expected Base to be an integer";
     if instance(C,Complex) then (
	 if opts.Degeneracy == true then (degenmapHash := hashTable flatten for n from 1 to d-1 list (
	    for i from 0 to n list (
		(n,i) => degenMapi(n,i,C)
		)
	     ););
	 facemapHash := hashTable flatten for n from 1 to d list (
	     for i from 0 to n list (
	 	 (n,i) => faceMapi(n,i,C)
	 	 )
	     );
	 --print("mde it here first");
	 if opts.Degeneracy == true then break return simplicialModule(C,facemapHash,degenmapHash,d);
	 --print("made it here");
	 return simplicialModule(C,facemapHash,d)
	 );
     )
 
 simplicialModule(Complex) := SimplicialModule => opts -> C -> (simplicialModule(C,length C,Degeneracy => opts.Degeneracy))
 

 simplicialModule(Module,ZZ) := SimplicialModule => opts -> (M,d) -> (simplicialModule(complex M,d,Degeneracy => opts.Degeneracy))
 
 simplicialModule(Ring,ZZ) := SimplicialModule => opts -> (R,d) -> (simplicialModule(R^1,d,Degeneracy => opts.Degneracy))
 
 forgetComplex = method(Options => {RememberSummands => true});
 forgetComplex(SimplicialModule) := SimplicialModule => opts -> S -> (
     if not any(keys S,i->i==symbol complex) then return S;
     L := new HashTable from for i to S.topDegree list i => combineSFactors(S,i,RememberSummands => opts.RememberSummands);
     faces := new HashTable from for i in keys S.dd.map list i => S.dd.map#i;
     if any(keys S,i->i==symbol ss) then (
	 degens := new HashTable from for i in keys S.ss.map list i => S.ss.map#i;
	 return simplicialModule(L,faces,degens,S.topDegree);
	 );
     D := simplicialModule(L,faces,S.topDegree);
     D.cache.components = components S;
     D
     )
 
 forgetDegeneracy = method();
 forgetDegeneracy(SimplicialModule) := S -> (
     if not any(keys S,i->i==symbol ss) then return S;
     if any(keys S,i->i==symbol complex) then return simplicialModule(S.complex,S.dd.map,S.topDegree);
     simplicialModule(S.module,S.dd.map,S.topDegree)
     )
 
 
SimplicialModule _ Sequence := Module => (S,p) -> (
    if #p =!= 2 then
    	error ("Expected a pair of integer indices");
    if S.module#?(p#0,p#1) then S.module#(p#0,p#1) else (ring S)^0
    )


SimplicialModule _ ZZ := Module => (S,n) -> (if S.module#?n then S.module#n else combineSFactors(S,n)) 



net SimplicialModule := S -> (
     (lo,hi) := (0,topDegree S);
     if lo > hi then 
         error "In a complex, lo <= hi should always hold in the concentration"
         --"0"
     else if lo == hi and C_lo === 0 then 
         "0"
     else if any(keys S,i->i==symbol complex) then
         (horizontalJoin between(" <-- ", 
             for i from lo to hi list
                 stack (net directSum(for k from 0 to min(i,S.complexLength) list (S.module)#(i,k)), " ", net i)))
     else 
         horizontalJoin between(" <-- ", 
             for i from lo to hi list
                 stack (net ((S.module)#i), " ", net i))  
     )
 
 
 Symbol ^ SimplicialModule := SimplicialModuleMap => (sym, C) -> (
    if sym === dd then C.dd;
    if sym === ss then C.ss
    else error "expected symbol to be 'dd' or 'ss'"
    )
 
lineOnTop := (s) -> concatenate(width s : "-") || s

source SimplicialModuleMap := SimplicialModule => f -> f.source
target SimplicialModuleMap := SimplicialModule => f -> f.target
ring SimplicialModuleMap := SimplicialModule => f -> ring source f
degree SimplicialModuleMap := ZZ => f -> f.degree

isHomogeneous SimplicialModuleMap := (f) -> all(values f.map, isHomogeneous)

simplicialModuleMap = method(Options => {Degeneracy => false});
simplicialModuleMap(ComplexMap,ZZ) := SimplicialModuleMap => opts -> (phi,d) -> (
    src := simplicialModule(source phi,d,opts);
    trg := simplicialModule(target phi,d,opts);
    result := map(trg,src,new HashTable from for i to d list i => directSum apply(mapMaker(phi,i),j->j_1),Degree => degree phi);
    result.cache.complexMap  = phi;
    result
    )

simplicialModuleMap(ComplexMap) := SimplicialModuleMap => opts -> phi -> (tDeg := max((source phi).concentration_1,(target phi).concentration_1);
    simplicialModuleMap(phi,tDeg)
    )
    



map(SimplicialModule, SimplicialModule, HashTable) := SimplicialModuleMap => opts -> (tar, src, maps) -> (
    if not(topDegree tar == topDegree src) then error "expected source and target to have the same top degree";
    R := ring tar;
    if ring src =!= R or any(values maps, f -> ring f =!= R) then
        error "expected source, target and maps to be over the same ring";
    deg := if opts.Degree === null 
           then 0 
           else if instance(opts.Degree, ZZ) then 
             opts.Degree
           else
             error "expected integer degree";
    (lo,hi) := (0,topDegree tar);
    maps' := hashTable for k in keys maps list (
    	local f;
        if instance(k, Sequence) then (
        f = maps#k;
        -- note: we use != instead of =!= in the next 2 tests,
        -- since we want to ignore any term order differences
	--print(k);
	--print(source f);
	--print(src_(first k));
        if rank source f != rank src_(first k) then (
            error ("map with index "|toString(k)|" has inconsistent source");
	);
        if rank target f !=  rank tar_(first(k)+deg) then
            error ("map with index "|toString(k)|" has inconsistent target");
        if first k < lo or first k > hi then continue else (k,f)
	)
       else (
	    f = maps#k;
        -- note: we use != instead of =!= in the next 2 tests,
        -- since we want to ignore any term order differences
	--print(k);
	--print(source f);
	--print(src_(first k));
        if source f !=  src_(k) then (
            error ("map with index "|toString(k)|" has inconsistent source");
	);
        if target f !=   tar_(k+deg) then
            error ("map with index "|toString(k)|" has inconsistent target");
        if k < lo or k > hi then continue else (k,f)
        ));
    new SimplicialModuleMap from {
        symbol source => src,
        symbol target => tar,
        symbol degree => deg,
        symbol map => maps',
        symbol cache => new CacheTable
        }
    )



SimplicialModuleMap _ ZZ := Matrix => (f,i) -> (
    if f.map#?i then f.map#i else map((target f)_(i + degree f), (source f)_i, 0))

SimplicialModuleMap _ Sequence := Matrix => (f,s) -> (
    if f.map#?s then f.map#s else map((target f)_(s#0 + degree f), (source f)_(s#0), 0))




expression SimplicialModuleMap := Expression => f -> (
    d := degree f;
    s := sort keys f.map;
    if #s === 0 then 
        new ZeroExpression from {0}
    else if instance(s_0,Sequence) then new VerticalList from for i in s list
        RowExpression {(i#0+d,i#1), ":", MapExpression { target f_i, source f_i, f_i }, ":", i}
    else if instance(s_0,ZZ) then return new VerticalList from for i in s list
        RowExpression {i+d, ":", MapExpression { target f_i, source f_i, f_i }, ":", i}
    )


 
 net SimplicialModuleMap := Net => f -> (
     v := between("",
            for i in sort keys f.map list (
                if instance(i,Sequence) then (horizontalJoin(
		            net ((i#0+f.degree,i#1)), " : ", net target f_i, " <--",
		            lineOnTop net f_i,
		            "-- ", net source f_i, " : ", net i
                    ))
	        else  (horizontalJoin(
		            net (i+f.degree), " : ", net target f_i, " <--",
		            lineOnTop net f_i,
		            "-- ", net source f_i, " : ", net i
			    ))
                ));
     if # v === 0 then net "0"
     else stack v
     )
 



SimplicialModule == ZZ := (C,n) -> (
    if n =!= 0 then error "cannot compare Simplicial module to non-zero integer";
    (lo,hi) := (0,C.topDegree);
    for i from lo to hi do if C_i != 0 then return false;
    true
    )
ZZ == SimplicialModule := (n,C) -> C == n
 
--as written, the code assumes one only takes direct sums of simplicial modules with same top degree
SimplicialModule.directSum = args -> (
    assert(#args > 0);
    R := ring args#0;
    if not all(args, C -> ring C === R) then error "expected all simplicial modules to be over the same ring";
    concentrations := for C in args list (0,C.topDegree);
    --checking if they all are Dold-Kan images
    allComplexes := all(args,i->any(keys i,j->j==symbol complex));
    allDegens := all(args,i->any(keys i, j->j== symbol ss));
    if not allComplexes then args = apply(args,i->forgetComplex(i));
    lo := concentrations/first//min;
    hi := concentrations/last//max;
    if not(all(args,i->i.topDegree==hi)) then error "all objects should have the same top degree";
    S := first args;
    LM := new HashTable from for i in keys S.module list i => directSum for j in args list if j.module#?i then j_i else continue;
    faceHashM := new HashTable from for i in keys S.dd.map list i => directSum for j in args list if j.dd.map#?i then j.dd_i;
    if allDegens then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => directSum for j in args list if j.ss.map#?i then j.ss_i;
	D := if allComplexes then simplicialModule(directSum for i in args list i.complex,faceHashM,degenMapHashM,hi)
	else  simplicialModule(LM,faceHashM,degenMapHashM,hi);
	D.cache.components = toList args;
	return D;
	);
    local D;
    D = if allComplexes then  simplicialModule(directSum for i in args list i.complex,faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree);
    D.cache.components = toList args;
    D    
    )
SimplicialModule ++ SimplicialModule := SimplicialModule => (C,D) -> directSum(C,D)
directSum SimplicialModule := C -> directSum(1 : C)

components SimplicialModule := C -> if C.cache.?components then C.cache.components else {C}

SimplicialModule#id = (C) -> (
    (lo,hi) := (0,C.topDegree);
    maps := hashTable for i from lo to hi list i => id_(C_i);
    result := map(C,C,maps);
    result.cache.isCommutative = true;
    result
    )


SimplicialModuleMap ^ ZZ := SimplicialModuleMap => (f,n) -> (
    tDeg := (source f).topDegree;
    df := degree f;
    if n === -1 then (
        maps := hashTable for i from 0 to tDeg list (i+df) => (
            f_i^(-1)
            );
        result := map(source f, target f, maps, Degree=>-df);
        if f.cache.?isCommutative then result.cache.isCommutative = f.cache.isCommutative;
        result
	    )
    else if n < 0 then (f^-1)^(-n)
    else if n === 0 then id_(source f)
    else if n === 1 then f
    else (
      if source f != target f then error "expected source and target to be the same";
      maps = hashTable for i from 0 to tDeg list i => (
          s := f_i;
          j := 1;
          while j < n do (
              s = f_(i+j*df) * s;
              j = j+1;
              );
          if s == 0 then continue else s
          );
      result = map(source f, source f, maps, Degree=> n * df);
      if f.cache.?isCommutative then result.cache.isCommutative = f.cache.isCommutative;
      result
      )
  )


SimplicialModule ** SimplicialModule := SimplicialModule => (C,D) -> simplicialTensor(C,D)

--it seemed more efficient to directly define the tensor product with a module
--as opposed to converting the module into a simplicial object then using simplicialTensor
Module ** SimplicialModule := SimplicialModule => (M,S) -> (
    LM := new HashTable from for i in keys S.module list i => M**(S.module)#i;
    faceHashM := new HashTable from for i in keys S.dd.map list i => M**(S.dd.map)#i;
    if any(keys S,i->i==symbol ss) then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => M**(S.ss.map)#i;
	if any(keys S,i->i==symbol complex) then return simplicialModule(M**(S.complex),faceHashM,degenMapHashM,S.topDegree)
	else return simplicialModule(LM,faceHashM,degenMapHashM,S.topDegree);
	);
    if any(keys S,i->i==symbol complex) then return simplicialModule(M**(S.complex),faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree)
    )
    
SimplicialModule ** Module := SimplicialModule => (S,M) -> (
    LM := new HashTable from for i in keys S.module list i => (S.module)#i**M;
    faceHashM := new HashTable from for i in keys S.dd.map list i => (S.dd.map)#i**M;
    if any(keys S,i->i==symbol ss) then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => (S.ss.map)#i**M;
	if any(keys S,i->i==symbol complex) then return simplicialModule((S.complex)**M,faceHashM,degenMapHashM,S.topDegree)
	else return simplicialModule(LM,faceHashM,degenMapHashM,S.topDegree);
	);
    if any(keys S,i->i==symbol complex) then return simplicialModule((S.complex)**M,faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree)
    )

SimplicialModule ** Matrix := SimplicialModuleMap => (S, f) -> (
    if ring S =!= ring f then error "expected Simplicial module and Matrix over the same ring";
    src := S ** source f;
    tar := S ** target f;
    map(tar, src, new HashTable from for i to S.topDegree list i => map(tar_i, src_i, S_i ** f))
    )

Matrix ** SimplicialModule := SimplicialModuleMap => (f, C) -> (
    if ring S =!= ring f then error "expected Simplicial module and Matrix over the same ring";
    src := source f ** S;
    tar := target f ** S;
    map(tar, src, new HashTable from for i to S.topDegree list i => map(tar_i, src_i, f ** S_i))
    )

SimplicialModule ** Ring := SimplicialModule => (S,R) -> (
    LM := new HashTable from for i in keys S.module list i => R**(S.module)#i;
    faceHashM := new HashTable from for i in keys S.dd.map list i => R**(S.dd.map)#i;
    if any(keys S,i->i==symbol ss) then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => R**(S.ss.map)#i;
	if any(keys S,i->i==symbol complex) then return simplicialModule(R**(S.complex),faceHashM,degenMapHashM,S.topDegree)
	else return simplicialModule(LM,faceHashM,degenMapHashM,S.topDegree);
	);
    if any(keys S,i->i==symbol complex) then return simplicialModule(R**(S.complex),faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree)
    )

Ring ** SimplicialModule := SimplicialModule => (R,S) -> S ** R

RingMap SimplicialModule := SimplicialModule => (phi,S) -> (
    LM := new HashTable from for i in keys S.module list i => phi((S.module)#i);
    faceHashM := new HashTable from for i in keys S.dd.map list i => phi((S.dd.map)#i);
    if any(keys S,i->i==symbol ss) then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => phi((S.ss.map)#i);
	if any(keys S,i->i==symbol complex) then return simplicialModule(phi(S.complex),faceHashM,degenMapHashM,S.topDegree)
	else return simplicialModule(LM,faceHashM,degenMapHashM,S.topDegree);
	);
    if any(keys S,i->i==symbol complex) then return simplicialModule(phi(S.complex),faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree)
    )

tensor(RingMap, SimplicialModule) := SimplicialModule => {} >> opts -> (phi, S) -> (
    if source phi =!= ring S then error "expected the source of the ring map to be the ring of the simplicial module";
    LM := new HashTable from for i in keys S.module list i => tensor(phi,(S.module)#i);
    faceHashM := new HashTable from for i in keys S.dd.map list i => tensor(phi,(S.dd.map)#i);
    if any(keys S,i->i==symbol ss) then (
	degenMapHashM := new HashTable from for i in keys S.ss.map list i => tensor(phi,(S.ss.map)#i);
	if any(keys S,i->i==symbol complex) then return simplicialModule(tensor(phi,(S.complex)),faceHashM,degenMapHashM,S.topDegree)
	else return simplicialModule(LM,faceHashM,degenMapHashM,S.topDegree);
	);
    if any(keys S,i->i==symbol complex) then return simplicialModule(tensor(phi,(S.complex)),faceHashM,S.topDegree)
    else simplicialModule(LM,faceHashM,S.topDegree)
    )
tensor(SimplicialModule, RingMap) := SimplicialModule => {} >> opts -> (S, phi) -> tensor(phi, S)

RingMap ** SimplicialModule := SimplicialModule => (phi, S) -> tensor(phi, S)
SimplicialModule ** RingMap := SimplicialModule => (S, phi) -> tensor(phi, S)


SimplicialModuleMap | SimplicialModuleMap := SimplicialModuleMap => (f,g) -> (
    if target f =!= target g then error "expected targets to be the same";
    if (source f).topDegree =!= (source g).topDegree then error "expected sources to have same top degree";
    deg := degree f;
    if deg =!= degree g then error "expected maps with the same degree";
    map(target f, source f ++ source g, new HashTable from for i to (source f).topDegree list i => (f_i|g_i), Degree=>deg)
    )

SimplicialModuleMap || SimplicialModuleMap := SimplicialModuleMap => (f,g) -> (
    if source f =!= source g then error "expected sources to be the same";
    if (target f).topDegree =!= (target g).topDegree then error "expected targets to have same top degree";
    deg := degree f;
    if deg =!= degree g then error "expected maps with the same degree";
    map(target f ++ target g, source f, new HashTable from for i to (source f).topDegree list i => (f_i||g_i), Degree=>deg)
    )

SimplicialModule == SimplicialModule := (C,D) -> (
    if C === D then return true;
    if keys C =!= keys D then return false;
    if topDegree C =!= topDegree D then return false;
    if ring C =!= ring D then return false;
    for i from 0 to C.topDegree do (
        if C_i != D_i then return false;
        );
    for i in keys C.dd.map do (
	if C.dd.map#i != D.dd.map#i then return false;
	);
    if any(keys C,i->i==symbol ss) then for i in keys C.ss.map do (
	if C.ss.map#i != D.ss.map#i then return false;
	);
    true    
    )

SimplicialModule == ZZ := (C,n) -> (
    if n =!= 0 then error "cannot compare Complex to non-zero integer";
    for i from 0 to C.topDegree do if C_i != 0 then return false;
    true
    )
ZZ == SimplicialModule := (n,C) -> C == n

transs := (C,v) -> (
    if C.cache.?indexComponents then (
	    Ci := C.cache.indexComponents;
	    apply(v, i -> if Ci#?i then Ci#i else error "expected an index of a component of the direct sum"))
    else (
        if not C.cache.?components then error "expected a direct sum of simplicialmodules";
	    Cc := C.cache.components;
	    apply(v, i -> if not Cc#?i then error "expected an index of a component of the direct sum");
	    v)
    )


SimplicialModule _ Array := SimplicialModuleMap => (C,v) -> (
    v = transs(C,v);
    D := directSum apply(toList v, j -> C.cache.components#j);
    Cc := if any(keys C,i->i==symbol complex) then forgetComplex(C,RememberSummands => false) else C;
    maps := hashTable for i from 0 to Cc.topDegree list i => map(C_i,D_i,Cc_i_v);
    result := map(C,D,maps);
    result.cache.isCommutative = true;
    result
    )

SimplicialModule ^ Array := SimplicialModuleMap => (C,v) -> (
    v = transs(C,v);
    D := directSum apply(toList v, j -> C.cache.components#j);
    Cc := if any(keys C,i->i==symbol complex) then forgetComplex(C,RememberSummands => false) else C;
    maps := hashTable for i from 0 to Cc.topDegree list i => map(D_i,C_i,Cc_i^v);
    result := map(D,C,maps);
    result.cache.isCommutative = true;
    result
    )


SimplicialModuleMap == SimplicialModuleMap := (f,g) -> (
    if f === g then return true;    
    if source f != source g or target f != target g 
      then return false;
    for i from 0 to (source f).topDegree do (
        if f_i != g_i then return false;
        );
    true    
    )
SimplicialModuleMap == ZZ := Boolean => (f,n) -> (
    if n === 0 then 
        all(keys f.map, k -> f.map#k == 0)
    else if n === 1 then (
        if source f != target f then return false;
        if degree f =!= 0 then return false;
        (lo,hi) := (0,(source f).topDegree);
        for i from lo to hi do
            if f_i != 1 then return false;
        f.cache.isCommutative = true;  -- this is the identity, after all!        
        true
        )
    else 
        error "cannot compare ComplexMap to integer other than 0 or 1"
    )
ZZ == SimplicialModuleMap := Boolean => (n,f) -> f == n

RingElement * SimplicialModuleMap := (r,f) -> (
    df := degree f;
    maps := hashTable for i to (source f).topDegree list i => (
        h := r * f_i;
        if h == 0 then continue else h
        );
    result := map(target f, source f, maps, Degree=>df);
    result
    )

SimplicialModuleMap * RingElement := (f,r) -> (r*f)

Number * SimplicialModuleMap := (r,f) -> (
    try r = promote(r,ring f) else error "can't promote scalar to ring of complex homomorphism";
    r * f
    )

SimplicialModuleMap * Number := (f,r) -> (
    try r = promote(r,ring f) else error "can't promote scalar to ring of complex homomorphism";
    f * r
    )

- SimplicialModuleMap := (f) -> (
    result := (-1)*f;
    if isCommutativeCached f then
        result.cache.isCommutative = true;
    result
    )

SimplicialModuleMap + SimplicialModuleMap := (f,g) -> (
    df := degree f;
    dg := degree g;
    if source f != source g then error "expected simplicial module homomorphisms with the same source";
    if target f != target g then error "expected simplicial homomorphisms with the same target";
    if df =!= dg then error "expected complex homomorphisms with the same degree";
    maps := hashTable for i from 0 to (source f).topDegree list i => (
        h := f_i + g_i;
        if h == 0 then continue else h
        );
    result := map(target f, source f, maps, Degree=>df);
    result
    )
SimplicialModuleMap + Number :=
SimplicialModuleMap + RingElement := SimplicialModuleMap => (f,r) -> (
    if r == 0 then f
    else (
        if source f != target f
        then error "expected same source and target"
        else f + r*id_(target f))
    )
Number + SimplicialModuleMap :=
RingElement + SimplicialModuleMap := SimplicialModuleMap => (r,f) -> f + r

SimplicialModuleMap - Number :=
SimplicialModuleMap - RingElement :=
SimplicialModuleMap - SimplicialModuleMap := SimplicialModuleMap => (f,g) -> f + (-1)*g

Number - SimplicialModuleMap :=
RingElement - SimplicialModuleMap := SimplicialModuleMap => (r,f) -> -f + r

SimplicialModuleMap * SimplicialModuleMap := (f,g) -> (
    df := degree f;
    dg := degree g;
    maps := hashTable for i from 0 to (source g).topDegree list i => (
        h := f_(dg + i) * g_i;
        if h == 0 then continue else h
        );
    result := map(target f, source g, maps, Degree=>df+dg);
    result
    )

--need to complete this
SimplicialModuleMap.directSum = args -> (
    -- args: sequence of SimplicialModuleMap's
    -- args: f_i : C_i --> D_i, having same degree deg
    -- result : sum(C_i) --> sum(D_i)
    R := ring args#0;
    deg := degree args#0;
    if not all(args, f -> ring f === R) then 
        error "expected maps all over the same ring";
    if not all(args, f -> degree f === deg) then
        error "expected maps to all have the same degree";
    -- WARNING: we call simplicialModule.directSum directly rather than using
    -- just directSum to avoid getting a cached copy of the direct
    -- sum.  Otherwise the labels of the cached copies might get
    -- changed (in Options.directSum).
    src := SimplicialModule.directSum (args/source);
    tar := SimplicialModule.directSum (args/target);
    -- only keep matrices in the homomorphism that are non-zero
    spots := unique flatten(args/(f -> keys f.map));
    maps := hashTable for i in spots list i => directSum(args/(f -> f_i));
    result := map(tar,src,maps,Degree=>deg);
    result.cache.components = toList args;
    result
    )

SimplicialModuleMap ++ SimplicialModuleMap := SimplicialModuleMap => (f,g) -> directSum(f,g)
directSum SimplicialModuleMap := f -> directSum(1 : f)
components SimplicialModuleMap := f -> if f.cache.?components then f.cache.components else {f}
SimplicialModuleMap ^ Array := SimplicialModuleMap => (f,v) -> (target f)^v * f
SimplicialModuleMap _ Array := SimplicialModuleMap => (f,v) -> f * (source f)_v



isCommutative SimplicialModuleMap := Boolean => f -> (
    if debugLevel == 0 and f.cache.?isCommutative then 
       return f.cache.isCommutative;
    C := source f;
    D := target f;
    deg := degree f;
    (loC,hiC) := C.concentration;
    (loD,hiD) := D.concentration;
    for i from loC to hiC do (
        if i+deg-1 >= loD and i+deg-1 <= hiD then (
            if not (dd^D_(i+deg) * f_i == (-1)^deg * (f_(i-1) * dd^C_i))
            then (
                if debugLevel > 0 then (
                    << "-- block " << (i,i-1) << " fails to commute" << endl;
                    );
                f.cache.isCommutative = false;
                return false;
                )
            )
        );
    f.cache.isCommutative = true;
    true
    )

-- the following method is not exported:
isCommutativeCached = method()
isCommutativeCached ComplexMap := Boolean => f -> f.cache.?isCommutative and f.cache.isCommutative

-*
isComplexMorphism = method(TypicalValue => Boolean)
isComplexMorphism ComplexMap := (f) -> (
    if debugLevel > 0 and degree f =!= 0 then (
        << "-- the complex map has non-zero degree" << endl;
        return false;
        );
    degree f === 0 and isCommutative f
    )
*-
--------------------------------------------------------------------
-- tensor products of simplicial maps ------------------------------
--------------------------------------------------------------------
tensor(SimplicialModuleMap, SimplicialModuleMap) := SimplicialModuleMap => {} >> opts -> (f,g) -> (
    -- f : C1 --> C2, g : D1 --> D2
    -- f**g : C1**D1 --> C2**D2
    -- (f**g)_i : sum_j(C1_j ** D1_(i-j) --> C2_(j+df) ** D2_(i-j+dg))
    df := degree f;
    dg := degree g;
    src := (source f) ** (source g);
    tar := (target f) ** (target g);
    -- for the i-th matrix src_i --> tar_(i+df+dg)
    -- we make a table of matrices, and create a block matrix from that using "matrix" and "map"
    (lo,hi) := (0,src.topDegree);
    maps := hashTable for i from lo to hi list i => f_i**g_i;
    result := map(tar, src, maps, Degree=>df+dg);
    result    
    )
SimplicialModuleMap ** SimplicialModuleMap := SimplicialModuleMap => (f,g) -> tensor(f,g)
SimplicialModule ** SimplicialModuleMap := SimplicialModuleMap => (C,g) -> id_C ** g
SimplicialModuleMap ** SimplicialModule := SimplicialModuleMap => (f,D) -> f ** id_D
Module ** SimplicialModuleMap := SimplicialModuleMap => (M,g) -> (
    map(M**(target g),M**(source g),new HashTable from for i in keys (g.map) list i => M**(g.map#i))
    )
SimplicialModuleMap ** Module := SimplicialModuleMap => (f,N) -> (
    map((target g)**M,(source g)**M,new HashTable from for i in keys (g.map) list i => (g.map#i)**M)
    )
    

SimplicialModuleMap ** Ring := SimplicialModuleMap => (f,R) -> (
    map((target g)**R,(source g)**R,new HashTable from for i in keys (g.map) list i => (g.map#i)**R)
    )
Ring ** SimplicialModuleMap := SimplicialModuleMap => (R,f) -> f ** R

RingMap SimplicialModuleMap := SimplicialModuleMap => (phi,f) -> (
    map(phi target f, phi source f, new HashTable from for i in keys (f.map) list i=> phi((f.map)#i))
    )

tensor(RingMap, SimplicialModuleMap) := SimplicialModuleMap => {} >> opts -> (phi, f) -> (
    if source phi =!= ring f then error "expected the source of the ring map to be the ring of the complex map";
    map(tensor(phi, target f), tensor(phi, source f), new HashTable from for i in keys (f.map) list i=> tensor(phi,(f.map)#i))
    )
tensor(SimplicialModuleMap, RingMap) := SimplicialModuleMap => {} >> opts -> (f, phi) -> tensor(phi, f)

RingMap ** SimplicialModuleMap := SimplicialModuleMap => (phi, f) -> tensor(phi, f)
SimplicialModuleMap ** RingMap := SimplicialModuleMap => (f, phi) -> tensor(phi, f)

----------------------------------------------------------------------------------------
------------- some functionality for complexes acting on simplicial modules ------------

Complex ** SimplicialModule := SimplicialModule => (C,S) -> (simplicialModule(C,S.topDegree)**S)

SimplicialModule ** Complex := SimplicialModule => (S,C) -> (S**simplicialModule(C,S.topDegree))

ComplexMap ** SimplicialModuleMap := SimplicialModuleMap => (f,g) -> (tDeg := max((source g).topDegree,(target g).topDegree);
    simplicialModuleMap(f,tDeg)**g
    )

SimplicialModuleMap ** ComplexMap := SimplicialModuleMap => (g,f) -> (tDeg := max((source g).topDegree,(target g).topDegree);
    g**simplicialModuleMap(f,tDeg)
    )

Complex ** SimplicialModuleMap := SimplicialModuleMap => (C,f) -> (id_C**f)
SimplicialModuleMap ** Complex := SimplicialModuleMap => (f,C) -> (f**id_C)



kernel SimplicialModuleMap := SimplicialModule => opts -> f -> (
    -- Initialize output as a local
    local result;
    -- f : B --> C
    B := source f;
    modules := hashTable for i from 0 to B.topDegree list i => kernel f_i;
    inducedMaps = hashTable for i from 0 to B.topDegree list i => inducedMap(B_i, modules#i);
    facemaps = hashTable for i in keys (B.dd.map)  list i => (
	        b1 :=B.dd_i * inducedMaps#(i_0);
		b2 := map(target b1,source inducedMaps#(i_0-1),inducedMaps#(i_0-1));
                (b1) // b2
                );
    if any(keys B,i->i==symbol ss) then (
	degenmaps = hashTable for i in keys (B.ss.map) list i => (
	    b1 := B.ss_i * inducedMaps#(i_0);
	    b2 := map(target b1,source inducedMaps#(i_0+1),inducedMaps#(i_0+1));
	    (b1) // b2
	    );
	result =  simplicialModule(modules,facemaps,degenmaps,B.topDegree);
	result.cache.kernel = f;
	return result;
	);
    result = simplicialModule(modules,facemaps,B.topDegree);
    result.cache.kernel = f;
    result
    )
cokernel SimplicialModuleMap := SimplicialModule => f -> (
    -- Initialize output as local
    local result;
    -- f : B --> C
    C := target f;
    deg := degree f;
    modules = hashTable for i from 0 to C.topDegree list i => cokernel f_(i-deg);
    facemaps = hashTable for i in keys (C.dd.map) list i => (
                map(modules#(i_0-1), modules#(i_0), matrix C.dd_i)
                );
    if any(keys C,i->i==symbol ss) then (
	degenmaps = hashTable for i in keys (C.ss.map) list i => (
	    map(modules#(i_0+1), modules#(i_0), matrix C.ss_i)
                );
	    result =  simplicialModule(modules,facemaps,degenmaps,C.topDegree);
	    result.cache.cokernel = f;
	    return result;
	    );
    result = simplicialModule(modules,facemaps,C.topDegree);
    result.cache.cokernel = f;
    result
    )

image SimplicialModuleMap := SimplicialModule => f -> (
    -- f : B --> C
    local result;
    B := source f;
    C := target f;
    deg := degree f;
    modules := hashTable for i from 0 to C.topDegree list i => image f_(i-deg);
    inducedMaps = hashTable for i from 0 to B.topDegree list i => inducedMap(C_i, modules#i);
    facemaps = hashTable for i in keys (C.dd.map)  list i => (
	        b1 :=C.dd_i * inducedMaps#(i_0);
		b2 := map(target b1,source inducedMaps#(i_0-1),inducedMaps#(i_0-1));
                (b1) // b2
                );
    if any(keys C,i->i==symbol ss) then (
	degenmaps = hashTable for i in keys (C.ss.map) list i => (
	    b1 := C.ss_i * inducedMaps#(i_0);
	    b2 := map(target b1,source inducedMaps#(i_0+1),inducedMaps#(i_0+1));
	    (b1) // b2
	    );
	result =  simplicialModule(modules,facemaps,degenmaps,C.topDegree);
	result.cache.image = f;
	return result;
	);
    result = simplicialModule(modules,facemaps,C.topDegree);
    result.cache.image = f;
    result
    )

coimage SimplicialModuleMap := SimplicialModule => f -> (
    -- f : B: --> C
    local result;
    B := source f;
    modules = hashTable for i from 0 to B.topDegree list i => coimage f_(i);
    facemaps = hashTable for i in keys (B.dd.map) list i => (
                map(modules#(i_0-1), modules#(i_0), matrix B.dd_i)
                );
    if any(keys B,i->i==symbol ss) then (
	degenmaps = hashTable for i in keys (B.ss.map) list i => (
	    map(modules#(i_0+1), modules#(i_0), matrix B.ss_i)
                );
	    result =  simplicialModule(modules,facemaps,degenmaps,B.topDegree);
	    result.cache.coimage = f;
	    return result;
	    );
    result = simplicialModule(modules,facemaps,B.topDegree);
    result.cache.coimage = f;
    result
    )

tensorwithComponents = method();
tensorwithComponents(Module,Module) := (M,N) -> (
    T := if M.cache.?indexComponents then flatten table(keys (M.cache.indexComponents),toList(0..length components N-1),(u,v) -> {u,v})
    else  flatten table(toList(0..length components M-1),toList(0..length components N-1),(u,v) -> {u,v});
    result := M**N;
    if M.cache.?indexComponents then result.cache.components = apply(T,i->(M.cache.components#(M.cache.indexComponents#(i_0)))**(N.cache.components#(i_1)))
    else result.cache.components =  apply(T,i->(M.cache.components#(i_0))**(N.cache.components#(i_1)));
    result.cache.indexComponents = hashTable for i to length(T)-1 list flatten (T_i) => i;
    result
    )

tensorwithComponents(Matrix,Matrix) := (A,B) -> (
    srcA := source A;
    srcB := source B;
    trgA := target A;
    trgB := target B;
    map(tensorwithComponents(trgA,trgB),tensorwithComponents(srcA,srcB),A**B)
    )

tensorwithComponents List := L -> (
    if length L == 2 then return tensorwithComponents(L_0,L_1);
    Ln := for i from 2 to length L-1 list L_i;
    tensorwithComponents({tensorwithComponents(L_0,L_1)}|Ln)
    )



minimalPresentation SimplicialModule := 
prune SimplicialModule := SimplicialModule => opts -> (cacheValue symbol minimalPresentation)(C -> (
	R := ring C;
    -- opts is ignored here
    -- to be cached: in the input C: cache the result D
    --               in the result: cache pruningMap: D --> C
    faceKeys = keys C.dd.map;
    if any(keys C,i->i==symbol ss) then degenKeys := keys C.ss.map; 
    prunedMods = new MutableHashTable from for i to C.topDegree list i => R^0;
    prunedFaceMaps = new MutableHashTable from for i in faceKeys list i => map(R^0,R^0,0);
    if any(keys C,i->i==symbol ss) then prunedDegenMaps := new MutableHashTable from for i in degenKeys list i=>map(R^0,R^0,0); 
    nonzeros := select(0..C.topDegree, i -> minimalPresentation C_i != 0);
    D = if #nonzeros === 0 
         then (
             simplicialModule((ring C)^0,C.topDegree)
             )
         else (
             lo = min nonzeros;
             hi = max nonzeros;
	     for i from lo to hi do prunedMods#i = minimalPresentation(C_i);
             for i in select(faceKeys,i->(i_0>lo and i_0<=hi)) do prunedFaceMaps#i = minimalPresentation C.dd_i;
	     if any(keys C,i->i==symbol ss) then (
		 for i in select(degenKeys,i->(i_0>=lo and i_0<hi)) do prunedDegenMaps#i = minimalPresentation C.ss_i;
		 );
             nmMods := new HashTable from for i in keys prunedMods list i => prunedMods#i;
	     nmFaces := new HashTable from for i in keys prunedFaceMaps list i=> prunedFaceMaps#i;
	     if any(keys C,i->i==symbol ss) then (
		 nmDegens := new HashTable from for i in keys prunedDegenMaps list i=> prunedDegenMaps#i;
		 );
	     if any(keys C,i->i==symbol ss) then simplicialModule(nmMods,nmFaces,nmDegens,C.topDegree) else simplicialModule(nmMods,nmFaces,C.topDegree)
                 );
    -- create the isomorphism D --> C
    pruning := hashTable for i from 0 to C.topDegree list i => (minimalPresentation C_i).cache.pruningMap;
    D.cache.pruningMap = map(C,D,pruning);
    D.cache.pruningMap.cache.isCommutative = true;
    D
    ))


minimalPresentation SimplicialModuleMap := 
prune SimplicialModuleMap := SimplicialModuleMap => opts -> f -> (
    C := source f;
    if not C.cache.?pruningMap then f = f * (minimalPresentation C).cache.pruningMap;
    D := target f;
    if not D.cache.?pruningMap then f = (minimalPresentation D).cache.pruningMap^-1 * f;
    f
    )

inducedMap(SimplicialModule, SimplicialModule) := SimplicialModuleMap => opts -> (D,C) -> (
    -- compute f : C --> D the map induced by the identity matrix.
    deg := if opts.Degree === null then 0 else opts.Degree;
    (loC,hiC) := (0,C.topDegree);
    (loD,hiD) := (0,D.topDegree);
    maps := hashTable for i from max(loC,loD-deg) to min(hiC,hiD-deg) list i => inducedMap(D_(i+deg),C_i, Verify => opts.Verify);
    map(D,C,maps,Degree=>deg)
    )


