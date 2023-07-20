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
    	maxK = min (d, max(length source phi,length target phi));
	for k to maxK do (
	mapList#(d,k) = directSum toList(
	    binomial(d,k):(phi_k));
	    );  
    for i in (sort keys mapList) list (i,mapList#i)
    )

combineSFactors = method();
combineSFactors(SimplicialModule,ZZ) := (S,d) -> (directSum for i to min(d,S.complexLength) list (S.module)#(d,i))
    
--H1 is the face maps, H2 is the degeneracy maps
simplicialModule = method(Options => {Base=>0,Degeneracy => false})
simplicialModule(Complex,HashTable,HashTable,ZZ) := SimplicialModule => opts -> (C,H1,H2,d) -> (
    spots := sort keys H1;
    if #spots === 0 then
      error "expected at least one map";
    R := ring C;
    moduleList := new MutableHashTable;
    for b to d do (
    	maxK = min (b, length C);
	for k to maxK do (
	moduleList#(b,k) = directSum toList(binomial(b,k):(C_k));
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
	moduleList#(b,k) = directSum toList(binomial(b,k):(C_k));
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
    --print("made it to face");
    S.dd = map(S,S,H1,Degree=>-1);
    --print("we are almost done");
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
    --print("made it to face");
    S.dd = map(S,S,H1,Degree=>-1);
    --print("we are almost done");
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
 
 simplicialModule(Complex) := SimplicialModule => opts -> C -> (simplicialModule(C,length C))
 
 
SimplicialModule _ Sequence := Module => (S,p) -> (
    if #p =!= 2 then
    	error ("Expected a pair of integer indices");
    if S.module#?(p#0,p#1) then S.module#(p#0,p#1) else (ring S)^0
    )
SimplicialModule _ ZZ := Module => (S,n) -> directSum for i in keys (S.module) list (if first(splice(toSequence{i}))==n then S.module#(i) else continue)

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
    src := simplicialModule(source phi,d);
    trg := simplicialModule(target phi,d);
    map(trg,src,new HashTable from for i to d list i => directSum apply(mapMaker(phi,i),j->j_1),Degree => degree phi)
    )
    

-*map(SimplicialModule, SimplicialModule, HashTable) := SimplicialModuleMap => opts -> (tar, src, maps) -> (
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
        if not instance(k_0, ZZ) then error "expected integer keys";
        f := maps#k;
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
        );
    new SimplicialModuleMap from {
        symbol source => src,
        symbol target => tar,
        symbol degree => deg,
        symbol map => maps',
        symbol cache => new CacheTable
        }
    )*-

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
        if instance(k, Sequence) then (
        f := maps#k;
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
	    f := maps#k;
        -- note: we use != instead of =!= in the next 2 tests,
        -- since we want to ignore any term order differences
	--print(k);
	--print(source f);
	--print(src_(first k));
        if rank source f != rank src_(k) then (
            error ("map with index "|toString(k)|" has inconsistent source");
	);
        if rank target f !=  rank tar_(k+deg) then
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

map(SimplicialModule, SimplicialModule, List) := SimplicialModuleMap => opts -> (tar, src, maps) -> (
    -- case 1: maps is a (single) list of matrices (maps between components of the complex)
    -- case 2: maps is a double list of ComplexMap's
    --    in this case, if the maps all commute with differentials, and are diagonal, then
    --    we could declare the result to be commutative as well. Should we do this?
    --  Can tell, depending on the class of maps#0.
    (lo,hi) := (0,topDegree tar);
    if instance(maps#0, List) then (
        mapHash := hashTable for i from lo to hi list i => (
            h := maps#(i-lo);
            if h_0 == 0 then continue else h
            );
	--print("makde it here");
        return map(tar,src,mapHash,opts)
        );
    -- At this point, the first entry of 'maps' is a List.
    -- Check: it is a table of ComplexMap
    -*R := ring tar;
    if R =!= ring src then error "expected complexes over the same ring";
    if not isTable maps then error "expected a table of SimplicialModuleMaps";
    -- check: all entries which are SimplicialModuleMaps have the same homological degree
    deg := if opts.Degree === null 
           then null
           else if instance(opts.Degree, ZZ) then 
             opts.Degree
           else
             error "expected integer degree";
    degs := unique for f in flatten maps list 
        if instance(f,ZZdFactorizationMap) 
            then degree f 
            else continue;
    if #degs > 1 then error "expected all ZZdFactorizationMaps to have the same degree";
    if deg =!= null and #degs == 1 and degs#0 =!= deg then error "Supplied degree is incompatible with the ComplexMaps";
    if deg === null then deg = (if #degs == 1 then degs#0 else 0);
    -- At this point, we need to create (block) matrices for each component of the complex.
    mapHash = hashTable for i from lo to hi list i => (
        newmaps := applyTable(maps, f -> if instance(f,SimplicialModuleMap) then f_i else f);
        h := map(tar_(i+deg), src_i, matrix newmaps);
        if h == 0 then continue else h
        );
    map(tar,src,mapHash,opts, Degree=>deg)*-
    )

-*flatten(SimplicialModuleMap) := SimplicialModuleMap => phi -> (if instance((keys phi.map)#0,ZZ) then return phi
    else (*-
	

SimplicialModuleMap _ ZZ := Matrix => (f,i) -> (
    if f.map#?i then f.map#i else map((target f)_(i + degree f), (source f)_i, 0))

SimplicialModuleMap _ Sequence := Matrix => (f,s) -> (
    if f.map#?s then f.map#s else map((target f)_(s#0 + degree f), (source f)_(s#0), 0))


    
    
-*expression SimplicialModuleMap := Expression => f -> (
    d := degree f;
    s := sort keys f.map;
    if #s === 0 then 
        new ZeroExpression from {0}
    else new VerticalList from for i in s list
        RowExpression {(i#0+d,i#1), ":", MapExpression { target f_i, source f_i, f_i }, ":", i}
    )*-

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

-*net SimplicialModuleMap := Net => f -> (
     v := between("",
            for i in sort keys f.map list (
                horizontalJoin(
		            net ((i#0+f.degree,i#1)), " : ", net target f_i, " <--",
		            lineOnTop net f_i,
		            "-- ", net source f_i, " : ", net i
                    )
                ));
     if # v === 0 then net "0"
     else stack v
     )*-
 
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
 
 
 
