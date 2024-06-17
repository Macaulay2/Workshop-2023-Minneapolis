doc ///
    Key
        (ZZdfactorization, HashTable)
    Headline
        make a ZZ/d-graded factorization
    Usage
        ZZdfactorization H
    Inputs
        H:HashTable
            each key is an integer indexing a differential, and the 
            value at that key is the map
        Base => ZZ
            ignored when the input is a hash table
    Outputs
        :ZZdFactorization
    Description
        Text
            A complex is a sequence of objects (e.g. modules),
            connected by maps called differentials.  The composition
            of any two consecutive maps is zero.
      
            The same data type is used for both chain and cochain
            complexes.  If {\tt C} is a complex, then we have 
            {\tt C^i = C_{-i}}.

            We construct the Koszul complex on the generators for the
            ideal of the twisted cubic curve.
        Example
            S = ZZ/101[a..d]
            I = ideal(b^2-a*c, b*c-a*d, c^2-b*d)
            F1 = map(S^1,,matrix{{I_0, I_1, I_2}})
            F2 = map(source F1,,matrix{
                    {0, I_2, -I_1},
                    {-I_2, 0, I_0},
                    {I_1, -I_0, 0}
                    })
            F3 = map(source F2,,matrix{{I_0}, {I_1}, {I_2}})
            C = complex hashTable{1 => F1, 2 => F2, 3 => F3}
            isWellDefined C
        Text
            This is the primary constructor used by all of the more
            user friendly methods for constructing a chain complex.
    Caveat
        This constructor minimizes computation
        and does very little error checking. To verify that a complex
        is well constructed, use @TO (isWellDefined, Complex)@.
    SeeAlso
        "Making chain complexes"
        (isWellDefined, Complex)
        (complex, List)
        (complex, Module)
///

doc ///
    Key
        ZZdfactorization
        (ZZdfactorization, List)
    Headline
        make a ZZ/d-graded factorization
    Usage
        ZZdfactorization L
    Inputs
        L:List
            of maps
        Base => ZZ
            the index of the target of the first map 
            in the differential.
    Outputs
        :ZZdFactorization
    Description
        Text
            A complex is a sequence of objects (e.g. modules),
            connected by maps called differentials.  The composition
            of any two consecutive maps is zero.
      
            The same data type is used for both chain and cochain
            complexes.  If {\tt C} is a complex, then we have
            {\tt C^i = C_{-i}}.

            Often, a complex is most easily described by giving a list
            of consecutive maps which form the differential.

            We construct the Koszul complex on the generators for the
            ideal of the twisted cubic curve.
        Example
            S = ZZ/101[a..d]
            I = ideal(b^2-a*c, b*c-a*d, c^2-b*d)
            F1 = map(S^1,,matrix{{I_0, I_1, I_2}})
            F2 = map(source F1,,matrix{
                    {0, I_2, -I_1},
                    {-I_2, 0, I_0},
                    {I_1, -I_0, 0}
                    })
            F3 = map(source F2,,matrix{{I_0}, {I_1}, {I_2}})
            C = complex {F1, F2, F3}
            isWellDefined C
        Text
            To start a complex at a base different from zero, use the
            optional argument {\tt Base}.
        Example
            C1 = complex({F1, F2, F3}, Base => 1)
            isWellDefined C1
        Text
            Notice that this changes the homological degrees of the
            maps, but is not the same as the shift of the complex
            (which for odd shifts negates the maps).
        Example
            dd^C1
            dd^(C[-1])
        Text
            Having constructed this complex, we can access individual
            terms and maps.
        Example
            C_2
            C^(-1)
            C^(-1) == C_1
            C_7
            dd^C
            dd^C_2
            length C
        Text
            By computing the homology of this complex, we see that
            these generators do not form a regular sequence, because
            $H_1(C)$ is non-zero.
        Example
            HH C
            prune HH C
            prune HH_1 C
    Caveat
        This constructor minimizes computation
        and does very little error checking. To verify that a complex
        is well constructed, use @TO (isWellDefined, Complex)@.
    SeeAlso
        "Making chain complexes"
        (isWellDefined, ZZdFactorization)
        (ZZdfactorization, HashTable)
        (ZZdfactorization, Module)
        (symbol SPACE, ZZdFactorization, Array)
///



doc ///
    Key
        (ring, ZZdFactorization)
        (ring, ZZdFactorizationMap)
    Headline
        access the ring of a ZZ/d-graded factorization or a factorization map
    Usage
        ring C
    Inputs
        C:ZZdFactorization
            or a @TO "ZZdFactorizationMap"@
    Outputs
        :Ring
    Description
        Text
            Every complex or complex map has a base ring.  This
            function access that information.
        Example
            S = ZZ/101[a,b,c,d];
            C = freeResolution coker vars S
            ring C
            assert(ring C === S)
            ring id_C
            assert(ring id_C === S)
    SeeAlso
        "Basic invariants and properties"
        ring
	--KELLER: do example where the ring changes after adjoining root
	--do an example in char p where you declare 1 to be the root of unity
///


doc ///
    Key
        period
        (period, ZZdFactorization)
    Headline
        The period of a ZZ/d-graded factorization
    Usage
        p = period C
    Inputs
        C:ZZdFactorization
    Outputs
        :ZZ
            the integer d for which C is a ZZ/d-graded factorization
    Description
        Text
            In this package, each factorization has a period d. 
        Example
            S = ZZ/101[a..c];
            C = freeResolution coker vars S
            concentration C
            D = C ++ C[5]
            concentration D
            min D
            max D
            assert((min D, max D) === concentration D)
        Text
            Indices that are outside of the concentration automatically
            return the zero object.
        Example
            C_-1
            D_4
        Text
            The function {\tt concentration} does no computation.
            To eliminate extraneous zeros, use @TO (prune, Complex)@.
        Example
            f1 = a*id_C  
            E = ker f1
            concentration E
            concentration prune E
        Text
            The concentration of a zero complex can be arbitrary, however,
            after pruning, its concentration will be {\tt (0,0)}.
        Example      
            C0 = (complex S^0)[4]
            concentration C0
            prune C0
            concentration oo
    SeeAlso
        "Basic invariants and properties"
        (symbol _, Complex, ZZ)
        (concentration, ComplexMap)
///


///
    Key
        (adjoinRoot,ZZ,Ring,Symbol)
	(adjoinRoot,ZZ,Ring,RingElement)
    Headline
        Adjoin a distinguished dth root of unity to a ring, with a specified name
    Usage
        adjoinRoot(d,R,t)
    Inputs
        d:ZZ
	    the order of the adjoined root of unity
	R:Ring
	    the ring to which the root of unity is adjoined
	t:Symbol or RingElement
	    the desired name of the distinguished root of unity
    Outputs
        S:Ring
	    the ring R but with the distinguished root of unity t adjoined
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (adjoinRoot,ZZdFactorization,Symbol)
	(adjoinRoot,ZZdFactorization,RingElement)
    Headline
        Adjoin a distinguished dth root of unity to a ZZ/d-graded factorization, with a specified name
    Usage
        adjoinRoot(C,t)
    Inputs
	C:ZZdFactorization
	    the ZZ/d-graded factorization to which the root of unity is adjoined
	t:Symbol or RingElement
	    the desired name of the distinguished root of unity
    Outputs
        D:Ring
	    the factorization C but with the distinguished root of unity t adjoined
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (adjoinRoot,ZZdFactorizationMap,Symbol)
	(adjoinRoot,ZZdFactorizationMap,RingElement)
    Headline
        Adjoin a distinguished dth root of unity to a ZZ/d-graded factorization map, with a specified name
    Usage
        adjoinRoot(phi,t)
    Inputs
	phi:ZZdFactorizationMap
	    the ZZ/d-graded factorization to which the root of unity is adjoined
	t:Symbol or RingElement
	    the desired name of the distinguished root of unity
    Outputs
        D:Ring
	    the factorization C but with the distinguished root of unity t adjoined
    Description
        Text
        Example
    Caveat
    SeeAlso
///


doc ///
   Key
     (isWellDefined, ZZdFactorization)
   Headline
     whether a ZZ/d-graded factorization is well-defined
   Usage
     isWellDefined C
   Inputs
     C:ZZdFactorization
   Outputs
     :Boolean
       that is true when {\tt C} determines a well defined complex
   Description
    Text
      This routine checks that the differential of {\tt C} composes to zero.
      Additionally, it checks that the underlying data in {\tt C} is a properly formed
      Complex object in Macaulay2. If the variable {\tt debugLevel} is set to a value greater than zero,
      then information about the nature of any failure is displayed.
    Text

      As a first example, we construct by hand the free resolution of the twisted
      cubic.  One must work with maps rather than matrices, because the source and the target
      of adjacent maps must be the same (including degree information).
    Example
      R = QQ[a..d];
      f0 = matrix {{-b^2+a*c, b*c-a*d, -c^2+b*d}}
      f1 = map(source f0,, {{d, c}, {c, b}, {b, a}})
      C = complex {f0, f1}
      isWellDefined C
      dd^C
      (dd^C)^2
    Text
      The zero complex is well-defined.
    Example
      C = complex R^0
      isWellDefined C
    Text
    
      The next example demonstrates the case when the sequence maps do not compose to 0.
    Example
      g1 = map(source f0,, {{-d, c}, {c, b}, {b, a}})
      C = complex {f0, g1}
      isWellDefined C
      debugLevel = 1
      isWellDefined C
      (dd^C)^2
   SeeAlso
     (isWellDefined, ComplexMap)
     map
///

doc ///
   Key
     (symbol _, ZZdFactorization, ZZ)
     (symbol ^, ZZdFactorization, ZZ)
   Headline
     access individual object in a ZZ/d-graded factorization
   Usage
     C_i
     C^i
   Inputs
     C:Complex
     i:ZZ
       either the homological or cohomological index
   Outputs
     :Module
       the {\tt i}-th object
   Description
    Text
       ZZ/d-graded factorizations can be either chain complexes or cochain complexes.  Subscripts
       refer to homological indices, and superscripts refer to
       cohomological indices.
     
       In this package homological indices are used by default.  For
       example, the @TO "concentration"@ references homological indices.
       Nevertheless, we always have the equation $C^i = C_{-i}$.
    Example
      S = ZZ/101[a..c]
      C = freeResolution coker vars S
      C_2
      C^(-2)
      C_2 == C^(-2)
    Text
      Indices that are outside of the concentration automatically
      return the zero object.
    Example
      C_-7
   SeeAlso
///

doc ///
   Key
     (symbol ==, ZZdFactorization, ZZdFactorization)
     (symbol ==, ZZdFactorization, ZZ)
     (symbol ==, ZZ, ZZdFactorization)
   Headline
     whether two ZZ/d-graded factorizations are equal
   Usage
     C == D
     C == 0
   Inputs
     C:ZZdFactorization
     D:ZZdFactorization
   Outputs
     :Boolean
       that is true when {\tt C} and {\tt D} are equal
   Description
    Text
      Two complexes are equal if the corresponding 
      objects and corresponding maps at each index are equal.
    Example
      S = ZZ/101[a..c]
      C = freeResolution coker vars S
      D = C[3][-3]
      C === D
      C == D
    Text
      Both the maps and the objects must be equal.
    Example
      (lo,hi) = concentration C
      E = complex for i from lo+1 to hi list 0*dd^C_i
      dd^E
      C == E
      E == 0
    Text
      A complex is equal to zero if all the objects and maps are zero.
      This could require computation to determine if something that
      is superficially not zero is in fact zero.
    Example
      f = id_C
      D = coker f
      D == 0
    Example
      C0 = complex S^0
      C1 = C0[4]
      concentration C0 == concentration C1
      C0 == C1
      C0 == 0
      C1 == 0
    Text
      Testing for equality is not the same testing for isomorphism.
      In particular, different presentations of a complex need not be equal.
    Example
      R = QQ[a..d];
      f0 = matrix {{-b^2+a*c, b*c-a*d, -c^2+b*d}}
      f1 = map(source f0,, {{d, c}, {c, b}, {b, a}})
      C = complex {f0, f1}
      HH C != complex coker f0
      prune HH C == complex coker f0
   Caveat
   SeeAlso
///

doc ///
    Key
        "differential of a ZZ/d-graded factorization"
        (symbol^, Symbol, ZZdFactorization)
    Headline
        get the maps between the terms in a ZZ/d-graded factorization
    Usage
        dd^C
        dd_C
    Inputs
        C:ZZdFactorization
    Outputs
        :ZZdFactorizationMap
            a map of degree -1
    Description
        Text
            A chain complex is a sequence of modules connected
            by homomorphisms, called differentials, such that
            the composition of any two consecutive maps is zero.
        Text
            One can access the differential of a complex as follows.
        Example
            R = QQ[a..d];
            I = ideal(a*d-b*c, b^2-a*c, c^2-b*d);
            C = freeResolution(R^1/I)
            dd^C
            C.dd
            assert(dd^C === C.dd)
            assert(source dd^C === C)
            assert(target dd^C === C)
            assert(degree dd^C === -1)
        Text
            The composition of the differential with itself is zero.
        Example
            (dd^C)^2 == 0
        Text
            The individual maps between terms are indexed by their
            source.
        Example
            dd^C_2
            assert(source dd^C_2 === C_2)
            assert(target dd^C_2 === C_1)
    SeeAlso
        "Making maps between chain complexes"
        (symbol_, ComplexMap, ZZ)
        (symbol_, Complex, ZZ)
        (source, ComplexMap)
        (target, ComplexMap)
        (degree, ComplexMap)
///

doc ///
    Key
        (gradedModule, ZZdFactorization)
    Headline
        a new ZZ/d-graded factorization in which the differential is zero
    Usage
        gradedModule C
    Inputs
        C:ZZdFactorization
    Outputs
        :ZZdFactorization
            whose differential is the zero map
    Description
        Text
            This routine isolates the terms in the complex
            and forgets the differentials
        Example
            R = ZZ/101[a,b,c,d,e];
            I = intersect(ideal(a,b),ideal(c,d,e))
            C = (dual freeResolution I)[-4]
            dd^C
            G = gradedModule C
            dd^G
            assert(isWellDefined G)
            assert(G != C)
        Text
            The homology of a complex already has zero differential.
        Example
            H = HH C
            prune H
            dd^H == 0
            assert(H == gradedModule H)
    SeeAlso
        (homology, Complex)
///

doc ///
   Key
     (homology, ZZdFactorization)
   Headline
     homology of a ZZ/d-graded factorization
   Usage
     H = HH C
   Inputs
     C:ZZdFactorization
   Outputs
     H:ZZdFactorization
   Description
    Text
      The homology factorization $H$ is defined by {\tt ker dd^C}/{\tt image dd^C}.
      The differential of the homology complex is the zero map.
      
      The first example is the complex associated to
      a triangulation of the real projective plane, having
      6 vertices, 15 edges, and 10 triangles.
    Example
      d1 = matrix {
          {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
          {-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}, 
          {0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 1, 0, 0, 0}, 
          {0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 1, 1, 0}, 
          {0, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, -1, 0, 1}, 
          {0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, -1, -1}}
      d2 = matrix {
          {-1, -1, 0, 0, 0, 0, 0, 0, 0, 0}, 
          {0, 0, -1, -1, 0, 0, 0, 0, 0, 0}, 
          {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
          {0, 1, 0, 0, -1, 0, 0, 0, 0, 0}, 
          {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, 
          {0, 0, 0, 0, 0, -1, -1, 0, 0, 0}, 
          {-1, 0, 0, 0, 0, 0, 0, -1, 0, 0}, 
          {0, -1, 0, 0, 0, 1, 0, 0, 0, 0}, 
          {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, 
          {0, 0, -1, 0, 0, 0, 0, 0, -1, 0}, 
          {0, 0, 0, 0, 0, -1, 0, 0, 1, 0}, 
          {0, 0, 0, -1, 0, 0, -1, 0, 0, 0}, 
          {0, 0, 0, 0, 0, 0, 0, 0, -1, -1}, 
          {0, 0, 0, 0, 0, 0, 0, -1, 0, 1}, 
          {0, 0, 0, 0, -1, 0, 0, 0, 0, -1}}
      C = complex {d1,d2}
      dd^C
      H = HH C
      dd^H == 0
    Text
      To see that the first homology group has torsion,
      we compute a minimal presentation of the homology.
    Example
      Hpruned = prune HH C
      dd^Hpruned == 0
    Text
      By dualizing the minimal free resolution of a monomial ideal,
      we get a free complex with non-trivial homology.  This particular
      complex is related to the local cohomology supported at the
      monomial ideal.
    Example
      S = ZZ/101[a..d, DegreeRank=>4];
      I = intersect(ideal(a,b),ideal(c,d))
      C = freeResolution (S^1/I)
      prune HH C
      Cdual = dual C
      prune HH Cdual
      prune HH_(-2) Cdual
   SeeAlso
     (dual, ZZdFactorization)
     (prune, ZZdFactorization)
///

doc ///
   Key
     (homology,ZZ,ZZdFactorization)
   Headline
     homology module of a ZZ/d-graded factorization
   Usage
     HH_i C
     HH^i C
   Inputs
     i:ZZ
     C:ZZdFactorization
   Outputs
     :Module
       the $i$-th homology or cohomology of the ZZ/d-graded factorization
   Description
    Text
      The $i$-th homology of a complex $C$ is the quotient
      ({\tt ker dd^C_i/image dd^C_(i+1)}).

      The first example is the complex associated to
      a triangulation of the real projective plane, having
      6 vertices, 15 edges, and 10 triangles.
    Example
      d1 = matrix {
          {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
          {-1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}, 
          {0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 1, 0, 0, 0}, 
          {0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 1, 1, 0}, 
          {0, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, -1, 0, 1}, 
          {0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, -1, 0, -1, -1}};
      d2 = matrix {
          {-1, -1, 0, 0, 0, 0, 0, 0, 0, 0}, 
          {0, 0, -1, -1, 0, 0, 0, 0, 0, 0}, 
          {1, 0, 1, 0, 0, 0, 0, 0, 0, 0}, 
          {0, 1, 0, 0, -1, 0, 0, 0, 0, 0}, 
          {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, 
          {0, 0, 0, 0, 0, -1, -1, 0, 0, 0}, 
          {-1, 0, 0, 0, 0, 0, 0, -1, 0, 0}, 
          {0, -1, 0, 0, 0, 1, 0, 0, 0, 0}, 
          {0, 0, 0, 0, 0, 0, 1, 1, 0, 0}, 
          {0, 0, -1, 0, 0, 0, 0, 0, -1, 0}, 
          {0, 0, 0, 0, 0, -1, 0, 0, 1, 0}, 
          {0, 0, 0, -1, 0, 0, -1, 0, 0, 0}, 
          {0, 0, 0, 0, 0, 0, 0, 0, -1, -1}, 
          {0, 0, 0, 0, 0, 0, 0, -1, 0, 1}, 
          {0, 0, 0, 0, -1, 0, 0, 0, 0, -1}};
      C = complex {d1,d2}
      dd^C
      HH C
      prune HH_0 C
      prune HH_1 C
      prune HH_2 C
    Text
      The $i$-th cohomology of a complex $C$ is the $(-i)$-th
      homology of $C$.
    Example
      S = ZZ/101[a..d, DegreeRank=>4];
      I = intersect(ideal(a,b),ideal(c,d))
      C = dual freeResolution (S^1/I)
      prune HH^1 C
      prune HH^2 C
      prune HH^3 C
   SeeAlso
     prune
     (dual, Complex)
///

doc ///
   Key
     (directSum, ZZdFactorization)
     (symbol++, ZZdFactorization, ZZdFactorization)
   Headline
     direct sum of complexes
   Usage
     D = C1 ++ C2
     D = directSum(C1,C2,...)
     D = directSum(name1 => C1, name2 => C2, ...)
   Inputs
     Ci:ZZdFactorization
   Outputs
     D:ZZdFactorization
       the direct sum of the input ZZ/d-graded factorizations
   Description
    Text
      The direct sum of two factorizations is another factorization, assuming the inputs factored the same ring element
    Example
      S = ZZ/101[a,b,c];
      C1 = freeResolution coker vars S
      C1 ++ complex(S^13)[-2]
      C2 = complex (ideal(a,b,c))
      C1 ++ C2
      assert isWellDefined(C1 ++ C2)
    Text
      The direct sum of a sequence of complexes can be computed as follows.
    Example
      C3 = directSum(C1,C2,C2[-2])
      assert isWellDefined C3
    Text
      The direct sum is an n-ary operator with projection and
      inclusion maps from each component satisfying appropriate
      identities.
    Example
      C4 = directSum(first => C1, second => C2)
      C4_[first] -- inclusion map C1 --> C4
      C4^[first] -- projection map C4 --> C1
      C4^[first] * C4_[first] == 1
      C4^[second] * C4_[second] == 1
      C4^[first] * C4_[second] == 0
      C4^[second] * C4_[first] == 0
      C4_[first] * C4^[first] + C4_[second] * C4^[second] == 1
    Text
      There are two short exact sequences associated to a direct sum.
    Example
      isShortExactSequence(C4^[first], C4_[second])
      isShortExactSequence(C4^[second], C4_[first])
    Text
      Given a complex which is a direct sum, we obtain the component
      complexes and their names (indices) as follows.
    Example
      components C3
      indices C3
      components C4
      indices C4
   SeeAlso
     (components,Complex)
     indices
     (symbol^, Complex, Array)
     (symbol_, Complex, Array)
     (isShortExactSequence, ComplexMap, ComplexMap)
     (sum, Complex)
///

doc ///
   Key
     (symbol_, ZZdFactorization, Array)
     (symbol^, ZZdFactorization, Array)
   Headline
     the canonical inclusion or projection map of a direct sum
   Usage
     i = C_[name]
     p = C^[name]
   Inputs
     C:ZZdFactorization
     name:
   Outputs
     :ZZdFactorizationMap
       {\tt i} is the canonical inclusion and {\tt p} is
       the canonical projection
   Description
    Text
      The direct sum is an n-ary operator with projection and
      inclusion maps from each component satisfying appropriate
      identities.

      One can access these maps as follows.      
    Example
      S = ZZ/101[a,b,c];
      C1 = freeResolution coker vars S
      C2 = complex (ideal(a,b,c))
      D = C1 ++ C2
      D_[0]
      D_[1]
      D^[0] * D_[0] == 1
      D^[1] * D_[1] == 1
      D^[0] * D_[1] == 0
      D^[1] * D_[0] == 0
      D_[0] * D^[0] + D_[1] * D^[1] == 1
    Text
      The default names for the components are the non-negative
      integers.  However, one can choose any name.
    Example
      E = (mike => C1) ++ (greg => C2)
      E_[mike]
      E_[greg]
      E^[mike] * E_[mike] == 1
      E^[greg] * E_[greg] == 1
      E^[mike] * E_[greg] == 0
      E^[greg] * E_[mike] == 0
      E_[mike] * E^[mike] + E_[greg] * E^[greg] == 1
    Text
      One can also access inclusion and projection maps of sub-direct sums.
    Example
      F = directSum(C1, C2, (complex S^13)[-4])
      F^[0,1]
      F_[0,2]
   SeeAlso
     (directSum, Complex)
     (components, Complex)
     indices
///

doc ///
   Key
     (components, ZZdFactorization)
   Headline
     list the components of a direct sum
   Usage
     components C
   Inputs
     C:ZZdFactorization
   Outputs
     :List
       the component factorizations of a direct sum (of ZZ/d-graded factorizations)
   Description
    Text
      A ZZ/d-graded factorization which has been constructed as a direct sum
      stores its component factorizations.
    Example
      S = ZZ/101[a,b,c];
      C1 = freeResolution coker vars S
      C2 = complex (ideal(a,b,c))
      D = C1 ++ C2
      L = components D
      L_0 === C1
      L_1 === C2
      E = (mike => C1) ++ (greg => C2)
      components E
    Text
      The names of the component complexes are called indices, 
      and are used to access the relevant inclusion and projection maps.
    Example
      indices D
      D^[0]
      indices E
      E_[greg]
   SeeAlso
     (directSum, ZZdFactorization)
     indices
     (symbol_, ZZdFactorization, Array)
     (symbol^, ZZdFactorization, Array)
///

doc ///
    Key
        (symbol SPACE, RingMap, ZZdFactorization)
    Headline
        apply a ring map
    Usage
        phi C
    Inputs
        phi:RingMap
            whose source is a ring $R$, and whose target is a ring $S$
        C:ZZdFactorization
            over the ring $R$
    Outputs
        :ZZdFactorization
            over the ring $S$
    Description
        Text
            We illustrate the image of a ZZ/d-graded factorization under a ring map.
        Example
            R = QQ[x,y,z]
            S = QQ[s,t]
            phi = map(S, R, {s, s+t, t})
            I = ideal(x^3, x^2*y, x*y^4, y*z^5)
            C = freeResolution I
            D = phi C
            isWellDefined D
            dd^D
            prune HH D
        Text
            When the ring map doesn't preserve homogeneity,
            the @TO "DegreeMap"@ option is needed to determine
            the degrees of the image free modules in the complex.
        Example
            R = ZZ/101[a..d]
            S = ZZ/101[s,t]
            phi = map(S, R, {s^4, s^3*t, s*t^3, t^4}, DegreeMap => i -> 4*i)
            C = freeResolution coker vars R
            D = phi C
            assert isWellDefined D
            assert isHomogeneous D
            prune HH D
    Caveat
        Every term in the ZZ?d-graded factorization must be free or a submodule of a free module.
        Otherwise, use @TO (tensor, RingMap, ZZdFactorization)@.
    SeeAlso
        (symbol SPACE, RingMap, ZZdFactorizationMap)
        (symbol **, RingMap, ZZdFactorization)
///

doc ///
    Key
        (symbol**, RingMap, ZZdFactorization)
        (symbol**, ZZdFactorization, RingMap)
        (tensor, RingMap, ZZdFactorization)
        (tensor, ZZdFactorization, RingMap)
        (symbol**, ZZdFactorization, Ring)
        (symbol**, Ring, ZZdFactorization)
    Headline
        tensor a ZZ/d-graded factorization along a ring map
    Usage
        phi ** C
        tensor(phi, C)
        S ** C
        C ** S
    Inputs
        phi:RingMap
            whose source is a ring $R$ and whose target is a ring $S$
        C:ZZdFactorization
            over the ring $R$
    Outputs
        :ZZdFactorization
            over the ring $S$
    Description
        Text
            These methods implement the base change of rings.  As input, one can either
            give a ring map $\phi$, or the ring $S$ (when there is a canonical map
                from $R$ to $S$).
        Text
            We illustrate the tensor product of a complex along a ring map.
        Example
            R = QQ[x,y,z];
            S = QQ[s,t];
            phi = map(S, R, {s, s+t, t})
            I = ideal(x^3, x^2*y, x*y^4, y*z^5)
            C = freeResolution I
            D = phi ** C
            assert isWellDefined D
            dd^D
            prune HH D
        Text
            If a ring is used rather than a ring map, then the implicit
            map from the underlying ring of the complex to the given ring
            is used.
        Example
            A = R/(x^2+y^2+z^2);
            C ** A
            assert(map(A,R) ** C == C ** A)
        Text
            The commutativity of tensor product is witnessed as follows.
        Example
            assert(D == C ** phi)
            assert(C ** A == A ** C)
        Text
            When the modules in the complex are not free modules,
            this is different than the image of a complex 
            under a ring map.
        Example
            use R
            I = ideal(x*y, x*z, y*z);
            J = I + ideal(x^2, y^2);
            g = inducedMap(module J, module I)
            assert isWellDefined g
            C = complex {g}
            D1 = phi C
            assert isWellDefined D1
            D2 = phi ** C
            assert isWellDefined D2
            prune D1
            prune D2
        Text
            When the ring map doesn't preserve homogeneity,
            the @TO "DegreeMap"@ option is needed to determine
            the degrees of the image free modules in the complex.
        Example
            R = ZZ/101[a..d];
            S = ZZ/101[s,t];
            f = map(S, R, {s^4, s^3*t, s*t^3, t^4}, DegreeMap => i -> 4*i)
            C = freeResolution coker vars R
            D = f ** C
            D == f C
            assert isWellDefined D
            assert isHomogeneous D
            prune HH D
            C1 = Hom(C, image vars R)
            D1 = f ** C1
            isWellDefined D1
            assert isHomogeneous D1
    SeeAlso
        (symbol **, RingMap, ZZdFactorizationMap)
        (symbol SPACE, RingMap, ZZdFactorization)
///

doc ///
    Key
        (sum, ZZdFactorization)
        (sum, ZZdFactorizationMap)
    Headline
        make the direct sum of all terms
    Usage
        sum C
        sum f
    Inputs
        C:ZZdFactorization
            or {\tt f}, @ofClass ZZdFactorizationMap@
    Outputs
        :Module
            or @ofClass Matrix@, if the input is a ZZ/d-graded factorization map
    Description
        Text
            This is the forgetful functor from the
            category of chain complexes to the category of modules.
            A chain complex $C$ is sent to the direct sum 
            $\bigoplus_i C_i$ of its terms.
            A map of chain complexes $f \colon C \to D$ is sent to the
            direct sum $\bigoplus_i f_i \colon \bigoplus_i C_i \to \bigoplus_i D_i$.
        Example
            S = ZZ/101[a,b,c];
            C = koszulComplex {a,b,c}
            sum C
            assert(rank sum C == 2^3)
        Example
            f = randomFactorizationMap(C, C, InternalDegree => 1, Cycle => true)
            g = sum f
            assert(g^2 === sum f^2)
            assert(target g === sum target f)
            assert(source g === sum source f)
            h = sum dd^C
            assert(h^2 == 0)
    SeeAlso
        "Basic invariants and properties"
        (directSum, ZZdFactorization)
        randomFactorizationMap
///





///
    Key
    Headline 
    Usage
    Inputs
    Outputs
    Description
        Text
        Example
    Caveat
    SeeAlso
///
