ZZdFactorization = new Type of MutableHashTable --
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
  
 ZZdFactorizationMap = new Type of HashTable
  -- keys:
  --   degree: ZZ
  --   source: ZZ/d-graded factorization over a ring R
  --   target: ZZ/d-graded factorization over the same ring R
  --   maps themselves (HashTable of Matrices), keys lying in the concentration period of the source.
  --    not all of the keys maps#i, need be present.
  --    missing ones are presumed to be zero maps.
  --   cache: a CacheTable
  --    cache.isCommutative: whether this map commutes with the differentials
  --      not set until needed.  unset means we have not checked yet, 
  --          and the user hasn't declared it to be true/false yet.

ZZdFactorization.synonym = "ZZ/d graded factorization"
ZZdFactorizationMap.synonym = "map of ZZ/d-graded factorizations"
 
  
ring ZZdFactorization := Ring => C -> C.ring




ZZdfactorization = method(Options => {Base=>0})
ZZdfactorization HashTable := ZZdFactorization => opts -> maps -> (
    spots := sort keys maps;
    if #spots === 0 then
      error "expected at least one matrix";
    if not all(spots, k -> instance(k,ZZ)) then
      error "expected matrices to be labelled by integers";
    if not all(spots, k -> instance(maps#k,Matrix)) then
      error "expected hash table or list of matrices";
    R := ring maps#(spots#0);
    if not all(values maps, f -> ring f === R) then
      error "expected all matrices to be over the same ring";
    moduleList := new MutableHashTable;
    for k in spots do (
        if not moduleList#?(k-1) 
          then moduleList#(k-1) = target maps#k;
        moduleList#k = source maps#k;
        );
    C := new ZZdFactorization from {
           symbol ring => R,
           symbol module => new HashTable from moduleList,
           symbol period => #spots,
           symbol cache => new CacheTable
           };
    C.dd = map(C,C,maps,Degree=>-1);
    C
    )
ZZdfactorization List := ZZdFactorization => opts -> L -> (
    -- L is a list of matrices or a list of modules
    if not instance(opts.Base, ZZ) then
      error "expected Base to be an integer"; 
    if all(L, ell -> instance(ell,Matrix)) then (
        trg := target L#0;
	mapHash := hashTable for i from 0 to #L-1 list opts.Base+i+1 => map(trg, trg, L#i);
        return ZZdfactorization(mapHash, opts)
        );
    if all(L, ell -> instance(ell,Module)) then (
        R := ring L#0;
        if any(L, ell -> ring ell =!= R) then
            error "expected modules all over the same ring";
        moduleHash := hashTable for i from 0 to #L-1 list opts.Base + i => L#i;
        C := new ZZdFactorization from {
            symbol ring => R,
            symbol period => #L,
            symbol module => moduleHash,
            symbol cache => new CacheTable
            };
        C.dd = map(C,C,0,Degree=>-1);
        return C;
        );
    error "expected a list of matrices or a list of modules";
    )
ZZdfactorization Module := ZZdFactorization => opts -> (M) -> (
    if not instance(opts.Base, ZZ) then
      error "complex: expected base to be an integer";
    if M.cache.?ZZdFactorization and opts.Base === 0 then return M.cache.ZZdFactorization;
    C := new ZZdFactorization from {
           symbol ring => ring M,
           symbol period => opts.Base,
           symbol module => hashTable {opts.Base => M},
           symbol cache => new CacheTable
           };
    if opts.Base === 0 then M.cache.ZZdFactorization = C;
    C.dd = map(C,C,0,Degree=>-1);
    C
    ) 
ZZdfactorization Ring := ZZdFactorization => opts -> R -> ZZdfactorization(R^1, opts)
ZZdfactorization Ideal := ZZdFactorization => opts -> I -> ZZdfactorization(module I, opts)

-*ZZdfactorization ZZdFactorizationMap := Complex => opts -> f -> (
    if degree f === -1 then (
        if source f =!= target f then error "expected a differential";
        (lo,hi) := concentration source f;
        if lo === hi then return complex((source f)_lo, Base=>lo);
        newmaps := hashTable for i from lo+1 to hi list i => f_i;
        return complex newmaps
        );
    -- TODO: keep this??  implement it??  what is it??
    -- f : C --> C, degree -1, then return (C,f) as a complex.
    -- f : C --> C[-1], return (C,f) as a complex
    -- complex (0 * dd^C)
    )*-

ZZdFactorization _ ZZ := Module => (C,i) -> C.module#(i%C.period)
--C.module#?i then C.module#i else C.module#(i%C.period)
--(ring C)^0
ZZdFactorization ^ ZZ := Module => (C,i) -> C_(-i)


Symbol ^ ZZdFactorization := ZZdFactorizationMap => (sym, C) -> (
    if sym === dd then C.dd
    else error "expected symbol to be 'dd'"
    )

net ZZdFactorization := C -> (
     (lo,hi) := (0,C.period);
     if lo > hi then 
         error "In a complex, lo <= hi should always hold in the concentration"
         --"0"
     else if lo == hi and C_lo === 0 then 
         "0"
     else
         horizontalJoin between(" <-- ", 
             for i from lo to hi list
                 stack (net C_i, " ", net i))
     )



lineOnTop := (s) -> concatenate(width s : "-") || s

source ZZdFactorizationMap := ZZdFactorization => f -> f.source
target ZZdFactorizationMap := ZZdFactorization => f -> f.target
ring ZZdFactorizationMap := ZZdFactorization => f -> ring source f
degree ZZdFactorizationMap := ZZ => f -> f.degree

isHomogeneous ZZdFactorizationMap := (f) -> all(values f.map, isHomogeneous)

map(ZZdFactorization, ZZdFactorization, HashTable) := ZZdFactorizationMap => opts -> (tar, src, maps) -> (
    R := ring tar;
    if ring src =!= R or any(values maps, f -> ring f =!= R) then
        error "expected source, target and maps to be over the same ring";
    deg := if opts.Degree === null 
           then 0 
           else if instance(opts.Degree, ZZ) then 
             opts.Degree
           else
             error "expected integer degree";
    (lo,hi) := (0,src.period);
    maps' := hashTable for k in keys maps list (
        if not instance(k, ZZ) then error "expected integer keys";
        f := maps#k;
        -- note: we use != instead of =!= in the next 2 tests,
        -- since we want to ignore any term order differences
        if source f != src_k then
            error ("map with index "|k|" has inconsistent source");
        if target f != tar_(k+deg) then
            error ("map with index "|k|" has inconsistent target");
        if k < lo or k > hi then continue else (k,f)
        );
    new ZZdFactorizationMap from {
        symbol source => src,
        symbol target => tar,
        symbol degree => deg,
        symbol map => maps',
        symbol cache => new CacheTable
        }
    )
map(ZZdFactorization, ZZdFactorization, List) := ZZdFactorizationMap => opts -> (tar, src, maps) -> (
    -- case 1: maps is a (single) list of matrices (maps between components of the complex)
    -- case 2: maps is a double list of ComplexMap's
    --    in this case, if the maps all commute with differentials, and are diagonal, then
    --    we could declare the result to be commutative as well. Should we do this?
    --  Can tell, depending on the class of maps#0.
    (lo,hi) := (0,src.period);
    if not instance(maps#0, List) then (
        mapHash := hashTable for i from lo to hi list i => (
            h := maps#(i-lo);
            if h == 0 then continue else h
            );
        return map(tar,src,mapHash,opts)
        );
    -- At this point, the first entry of 'maps' is a List.
    -- Check: it is a table of ComplexMap
    R := ring tar;
    if R =!= ring src then error "expected complexes over the same ring";
    if not isTable maps then error "expected a table of ComplexMaps";
    -- check: all entries which are ComplexMaps have the same homological degree
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
        newmaps := applyTable(maps, f -> if instance(f,ZZdFactorizationMap) then f_i else f);
        h := map(tar_(i+deg), src_i, matrix newmaps);
        if h == 0 then continue else h
        );
    map(tar,src,mapHash,opts, Degree=>deg)
    )

map(ZZdFactorization, ZZdFactorization, Function) := ZZdFactorizationMap => opts -> (D,C,f) -> (
    deg := if opts.Degree === null then 0 else opts.Degree;
    (loC,hiC) := (0,period C);
    (loD,hiD) := (0,period D);
    maps := hashTable for i from max(loC,loD-deg) to min(hiC,hiD-deg) list (
        if C_i == 0 or D_(i+deg) == 0 then continue;
        g := f(i);
        if g === null or g == 0 then continue;
        i => g
        );
    map(D,C,maps,Degree=>deg)
    )

map(ZZdFactorization, ZZdFactorization, ZZ) := ZZdFactorizationMap => opts -> (D, C, j) -> (
    if j === 0 then (
        result := map(D,C,hashTable{},opts);
        result.cache.isCommutative = true;
        return result
        );
    if j === 1 then (
        if C == D and (opts.Degree === null or opts.Degree === 0) then
            return id_C;
        error "expected source and target to be the same";
        );
    error "expected integer to be zero or one";
    )

map(ZZdFactorization, ZZdFactorization, ZZdFactorizationMap) := ZZdFactorizationMap => opts -> (tar, src, f) -> (
    deg := if opts.Degree === null then degree f else opts.Degree;
    H := hashTable for k in keys f.map list k => map(tar_(deg+k), src_k, f.map#k);
    map(tar,src,H, Degree=>deg)
    )


ZZdFactorizationMap _ ZZ := Matrix => (f,i) -> (
    --f.map#(i%(source f).period))
    if i%((source f).period)==0 then return f.map#((source f).period);
    if f.map#?i then f.map#i else f.map#(i%(source f).period))
	--map((target f)_(i + degree f), (source f)_i, 0))
ZZdFactorizationMap ^ ZZ := ZZdFactorizationMap => (f,n) -> (
    (lo,hi) := (0,(source f).period);
    df := degree f;
    if n === -1 then (
        maps := hashTable for i from lo to hi list (i+df) => (
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
      maps = hashTable for i from lo to hi list i => (
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


expression ZZdFactorizationMap := Expression => f -> (
    d := degree f;
    s := sort keys f.map;
    if #s === 0 then 
        new ZeroExpression from {0}
    else new VerticalList from for i in s list
        RowExpression {i+d, ":", MapExpression { target f_i, source f_i, f_i }, ":", i}
    )

net ZZdFactorizationMap := Net => f -> (
     v := between("",
            for i in sort keys f.map list (
                horizontalJoin(
		            net (i+f.degree), " : ", net target f_i, " <--",
		            lineOnTop net f_i,
		            "-- ", net source f_i, " : ", net i
                    )
                ));
     if # v === 0 then net "0"
     else stack v
     )
 
 
--want a FOLD command. This will make the tensor work correctly
congClasses = method()
congClasses(ZZ,ZZ) := (d,e) -> (L=apply(0..e-1,i->(i,i%d));
    Ln := new MutableList;
    for i from 0 to d-1 do Ln#(#Ln) = toList apply(select(L,j->j_1==i),j->j_0);
    toList Ln
    )

--Fold = method()
--Fold(ZZdFactorization,ZZ) := (C,d) ->     


 
component = method()
component(Module,Thing) := (M,k) -> (
    if not M.cache.?indexComponents then error "expected Module to be a direct sum with indexed components";
    if not M.cache.indexComponents#?k then error("expected "|toString k|" to be the index of a component");
    (components M)#(M.cache.indexComponents#k)
    )


 
 tensor(ZZdFactorization, ZZdFactorization) := ZZdFactorization => {} >> opts -> (C, D) -> (
    Y := youngest(C,D);
    if Y.cache#?(tensor,C,D) then return Y.cache#(tensor,C,D);
    R := ring C;
    if ring D =!= R then error "expected factorizations over the same ring";
    (loC,hiC) := (0,C.period);
    (loD,hiD) := (0,D.period);
    modules := hashTable for i from 0 to C.period list i => (
        directSum for j from loC to 2*hiC-1 list (
            if (i-j) >= loD and (i-j)%d <= hiD then
                {j,i-j} => C_(j%d) ** D_((i-j)%d)
            else
                continue
            )
        );
    if loC === hiC and loD === hiD then (
        result := complex(modules#(loC+loD), Base => loC+loD);
        result.cache.tensor = (C,D);
        Y.cache#(tensor,C,D) = result;
        return result;
        );
    maps := hashTable for i from 1 to C.period list i => (
        map(modules#(i-1),
            modules#i,
            matrix table(
                indices modules#(i-1),
                indices modules#i,
                (j,k) -> (
                    tar := component(modules#(i-1), j);
                    src := component(modules#i, k);
                    m := map(tar, src, 
                        if k-j === {0,1} then (-1)^(k#0) * (C_(k#0) ** dd^D_(k#1))
                        else if k-j === {1,0} then (dd^C_(k#0) ** D_(k#1))
                        else 0);
                    m
                    ))));
    result = ZZdfactorization maps;
    result.cache.tensor = (C,D);
    Y.cache#(tensor,C,D) = result;
    result
    )
ZZdFactorization ** ZZdFactorization := ZZdFactorization => (C,D) -> tensor(C,D)
Module ** ZZdFactorization := ZZdFactorization => (M,D) -> (ZZdfactorization M) ** D
ZZdFactorization ** Module := ZZdFactorization => (C,N) -> C ** (ZZdfactorization N)

ZZdFactorization ** Matrix := ZZdFactorizationMap => (C, f) -> (
    if ring C =!= ring f then error "expected Complex and Matrix over the same ring";
    src := C ** source f;
    tar := C ** target f;
    map(tar, src, i -> map(tar_i, src_i, C_i ** f))
    )
Matrix ** ZZdFactorization := ZZdFactorizationMap => (f, C) -> (
    if ring C =!= ring f then error "expected Complex and Matrix over the same ring";
    src := source f ** C;
    tar := target f ** C;
    map(tar, src, i -> map(tar_i, src_i, f ** C_i))
    )

ZZdFactorization ** Ring := ZZdFactorization => (C,R) -> (
    (lo,hi) := (0,C.period);
    moduleHash := hashTable for i from lo to hi list i => C_i ** R;
    if lo === hi then 
        return complex(moduleHash#lo, Base=>lo);
    mapHash := hashTable for i from lo+1 to hi list i => 
        map(moduleHash#(i-1), moduleHash#i, (cover dd^C_i) ** R);
    ZZdfactorization mapHash
    )
Ring ** ZZdFactorization := ZZdFactorization => (R,C) -> C ** R

RingMap ZZdFactorization := ZZdFactorization => (phi,C) -> (
    (lo,hi) := (0,C.period);
    moduleHash := hashTable for i from lo to hi list i => phi C_i;
    if lo === hi then 
        return ZZdfactorization(moduleHash#lo, Base=>lo);
    mapHash := hashTable for i from lo+1 to hi list i => 
        map(moduleHash#(i-1), moduleHash#i, phi dd^C_i);
    ZZdfactorization mapHash
    )

tensor(RingMap, ZZdFactorization) := ZZdFactorization => {} >> opts -> (phi, C) -> (
    if source phi =!= ring C then error "expected the source of the ring map to be the ring of the complex";
    (lo,hi) := (0,C.period);
    modules := hashTable for i from lo to hi list i => tensor(phi, C_i);
    if lo === hi then 
        return ZZdfactorization(modules#lo, Base=>lo);
    maps := hashTable for i from lo+1 to hi list i => 
        map(modules#(i-1), modules#i, tensor(phi, matrix dd^C_i));
    ZZdfactorization maps
    )
tensor(ZZdFactorization, RingMap) := ZZdFactorization => {} >> opts -> (C, phi) -> tensor(phi, C)

RingMap ** ZZdFactorization := ZZdFactorization => (phi, C) -> tensor(phi, C)
ZZdFactorization ** RingMap := ZZdFactorization => (C, phi) -> tensor(phi, C)
-- Hello from juan
