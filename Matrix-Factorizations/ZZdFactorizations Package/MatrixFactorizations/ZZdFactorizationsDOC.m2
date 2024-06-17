
doc ///
    Key
        ZZdFactorization
    Headline
        A package for creating and computing objects in the category of ZZ/d-graded factorizations, such as matrix factorizations
    Description
        Text
            A ZZ/d-graded factorization F  of a ring element f is a ZZ/d-graded complex of free R-modules equipped with a degree -1 (mod d) endomorphism d^F such that (d^F)^d = f * id_F. 
	    In practice, a ZZdFactorization may be visualized as a sequence of R-module maps:
	    
	    $F_0 \leftarrow F_1 \leftarrow \cdots \leftarrow F_{d-1}$
	    
	    with the caveat that $d^F_0 : F_0 \to F_{d-1}$, since one should count degree modulo d. Any 2-periodic complex may be reinterpretted as a ZZ/2-graded factorization of 0, and 
	    likewise a matrix factorization of a ring element f is equivalently a ZZ/2-graded factorization of f. Because of their similarity with complexes, much of the functionality
	    and syntax for this package closely resembles the "Complexes" package, with some key differences that we will highlight below.
        Example
            Q = QQ[x_1..x_3];
	    F1 = ZZdfactorization {x_1,x_2}
	    F1.dd
	    F2 = ZZdfactorization {x_1,x_2,x_3}
	    F2.dd
	Text
	    Notice that in the above, both the modules and the differentials are displayed modulo the period of the complex. Moreover, if the user tries to access the data of the modules
	    or differentials, this input is also taken modulo the period:
	Example
	    F1_0
	    F1_2
	    F1.dd_123
	Text
	    The package is implemented to not actually check for well-definedness of the factorization (that is, it does not check if all of the differentials actually compose to a scalar
	    multiple of the identity). The user can check this by using the isWellDefined and isFactorization commands: 
	Example
	    isdFactorization F1
	Text
	    Much of the syntax and functionality of the types ZZdFactorization and ZZdFactorizationMap are based on current functionalities for the analogous objects in the Complexes package,
	    so there should be essentially no learning curve for users already familiar with working with chain complexes.*-
    SeeAlso
        "Making ZZdFactorizations"
        "Making maps between factorizations"
        "Basic invariants and properties"
	"Preprogrammed examples and operations"
///

doc ///
    Key
        "Making ZZdFactorizations"
    Headline
        Information about the basic constructors
    Description
    	Text
    	    @SUBSECTION "Basic constructors"@
	Text
    	    @UL {
                TO (ZZdfactorization, HashTable),
                TO (ZZdfactorization, List), 
                TO (isWellDefined, ZZdFactorization),
		TO (isdFactorization, ZZdFactorization)
            }@
    	Text
    	    @SUBSECTION "Important computations creating new ZZ/d-graded factorizations"@
	Text
	    --PUT KOSZUL FACTORIZATION/RANDOM FACTORIZATION HERE
    	    @UL {
                TO (resolution, Complex),
                TO (homology, Complex)
            }@
    	Text
    	    @SUBSECTION "More advanced constructors"@
	Text
    	    @UL {
                TO (symbol++, ZZdFactorization, ZZdFactorization),
                TO (symbol**, ZZdFactorization, ZZdFactorization),
                TO (Hom, ZZdFactorization, ZZdFactorization),
                TO (dual, ZZdFactorization),
                TO (symbol SPACE, RingMap, ZZdFactorization),
                TO (symbol **, RingMap, ZZdFactorization),
                --TO (koszulComplex, Matrix),
                TO (minimalPresentation, ZZdFactorization),
                --TO (minimize, ZZdFactorization),
            }@
    	Text
    	    @SUBSECTION "Extracting ZZ/d-graded factorization from ZZ/d-graded factorization maps"@
        Text
    	    @UL {
                TO (source, ZZdFactorizationMap),
                TO (target, ZZdFactorizationMap),
                TO (kernel, ZZdFactorizationMap),
                TO (cokernel, ZZdFactorizationMap),
                TO (image, ZZdFactorizationMap),
                TO (coimage, ZZdFactorizationMap),
                TO (cone, ZZdFactorizationMap),
            }@
    SeeAlso
        "Making maps between factorizations"
        "Basic invariants and properties"
	"Preprogrammed examples and operations"  
///


doc ///
    Key
        "Basic invariants and properties"
    Headline
        information about accessing basic features
    Description
    	Text
    	    @SUBSECTION "Predicates for complexes and complex maps"@
        Text
    	    @UL {
                TO (isWellDefined, ZZdFactorization),
                --TO (isFree,ZZdFactorization),
                TO (isWellDefined, ZZdFactorizationMap),
                TO (isCommutative, ZZdFactorizationMap),
                --TO (isQuasiIsomorphism, ZZdFactorizationMap),
                --TO (isShortExactSequence, ZZdFactorizationMap, ZZdFactorizationMap),
                --TO (isNullHomotopic, ZZdFactorizationMap),
                --TO (isNullHomotopyOf, ZZdFactorizationMap, ZZdFactorizationMap)
            }@
    	Text
    	    @SUBSECTION "Other invariants for ZZdFactorizations"@
        Text
    	    @UL {
                TO (ring,ZZdFactorization),
                TO (period,ZZdFactorization),
                TO (components, ZZdFactorization)
            }@
    	Text
    	    @SUBSECTION "Other invariants for complex maps"@
        Text
    	    @UL {
                TO (source, ZZdFactorizationMap),
                TO (target, ZZdFactorizationMap),
                TO (degree, ZZdFactorizationMap),
                TO (ring, ZZdFactorizationMap),
                TO (components,ZZdFactorizationMap),
            }@
    SeeAlso
        "Making chain complexes"
        "Making maps between factorizations"
        "Preprogrammed examples and operations"
///


doc ///
    Key
        ZZdFactorization
    Headline
        the class of all ZZ/d-graded factorizations
    Description
        Text
            A ZZ/d-graded factorization is a sequence of objects $C_i$, connected by
            maps $dd^C_i : C_i \rightarrow C_{i-1}$ such that the
            composition of any d consecutive maps is equal to a fixed scalar multiple of the identity.

            TODO: more needs to be added here explaining how to used complexes.
            and links to "landing pages".
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
        indices on which a complex may be non-zero
    Usage
        p = period C
    Inputs
        C:ZZdFactorization
    Outputs
        :ZZ
            the integer d for which C is a ZZ/d-graded factorization
    Description
        Text
            In this package, each factorization has a period d.  When {\tt lo <= i <= hi}, the module
            {\tt C_i} might be zero.  The methods {\tt max} and {\tt min} 
            applied to the complex {\tt C} return {\tt lo} and {\tt hi}, respectively.
      
            This function is mainly used in programming, to loop over all
            non-zero modules or maps in the complex.  This should not be confused
            with the support of a complex.
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

-*doc ///
    Key --THIS MUST TAKE THE PERIOD INPUT AS WELL
        (ZZdfactorization, Module)
        (ZZdactorization, Ideal)
        (ZZdfactorization, Ring)
    Headline
       view a module, ring, or ideal as a ZZ/d-graded factorization with a fixed period
    Usage
        complex M
    Inputs
        M:Module
            or @TO "Ideal"@, or @TO "Ring"@.
        Base => ZZ
            index for {\tt M}
    Outputs
        :ZZdFactorization
            returns the complex whose 0-th component is {\tt M}.
    Description
        Text
            In contrast to @TO (complex,HashTable)@ and @TO
            (complex,List)@, this constructor provides a convenient
            method to construct a complex with only one non-zero
            component.
      
            We illustrate this with a free module.
        Example
            S = ZZ/101[a..d]
            C0 = complex S^2
            f = dd^C0
            source f, target f
            f == 0
            isWellDefined C0
            C0 == 0
            length C0
        Example
            C1 = complex(S^2, Base=>3)
            C1 == C0[-3]
            C1_3
            C1_0
        Text
            A ring or an ideal will be converted to a module first.
        Example
            C2 = complex S
            I = ideal(a^2-b, c^3)
            C3 = complex I
            C4 = complex (S/I)
            (ring C3, ring C4)
        Text
            The zero complex over a ring {\tt S} is most conveniently
            created by giving the zero module.
        Example
            C5 = complex S^0
            C5 == 0
            dd^C5 == 0
            C5_0
    SeeAlso
        "Making chain complexes"
        (isWellDefined, Complex)
        (complex, HashTable)
///*-


///
    Key
        (isdFactorization,ZZdFactorization)
    Headline
        Check whether all d-fold compositions of the differentials of a ZZ/d-graded factorization
	compose to a scalar multiple of the identity, and outputs this scalar multiple.
    Usage
        isdFactorization(C)
    Inputs
        C:ZZdFactorization
    Outputs
        :Boolean
	    whether or not the d-fold compositions are a scalar multiple of the identity
	RingElement:f
	    the element f that the ZZ/d-graded factorization factors
    Description
        Text
        Example
            S=QQ[x,y,z]
            X = ZZdfactorization{x,y,z}
            isdFactorization(X)

            Y = ZZdfactorization{y,y}
            Z = ZZdfactorization{z,z}
            Y ++ Z 
    Caveat
    SeeAlso
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
        (complex, ComplexMap)
    Headline
        make a complex by specifying the differential
    Usage
        C = complex d
    Inputs
        d:ComplexMap
        Base => ZZ
            unused
    Outputs
        C:Complex
            whose differential is $d$
    Description
        Text
            Given a map $d$ of complexes having degree -1 and whose source 
            and targets are equal, this method constructs the chain complex
            whose differential is $d$.  This constructor does not verify that
            $d^2 = 0$.
        Example
            S = ZZ/101[x_1..x_4];
            F = freeResolution coker vars S
            d = randomComplexMap(F, F, Cycle => true, InternalDegree => -1, Degree => -1)
            d^2
            C = complex d
            assert isWellDefined C
            assert all(0..4, i -> dd^C_i == d_i)
        Example
            e = randomComplexMap(F, F, InternalDegree => -1, Degree => -1)
            D = complex e
            debugLevel = 1
            assert not isWellDefined D
    SeeAlso
        "Making chain complexes"
        (isWellDefined, Complex)
        (complex, HashTable)
        (complex, List)
///

-- TODO: Add programming details
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
     (symbol SPACE, ZZdFactorization, Array)
     (symbol SPACE, ZZdFactorizationMap, Array)
   Headline
     shift a ZZ/d-graded factorization or map of ZZ/d-graded factorizations
   Usage
     D = C[i]
     g = f[i]
   Inputs
     C:ZZdFactorization
       or {\tt f}, a @TO ZZdFactorizationMap@
     :Array
       {\tt [i]}, where {\tt i} is an integer
   Outputs
     D:ZZdFactorization
       or {\tt g}, a @TO ZZdFactorizationMap@.
   Description
    Text
      The shifted complex $D$ is defined by $D_j = C_{i+j}$ for all $j$
      and the sign of the differential is changed if $i$ is odd.
       
      The shifted complex map $g$ is defined by $g_j = f_{i+j}$ for all $j$.
    
      The shift defines a natural automorphism on the category of complexes. 
      Topologists often call the shifted complex $C[1]$ the {\it suspension} of $C$.
    Example
      S = ZZ/101[a..d]
      C = freeResolution coker vars S
      dd^C_3
      D = C[1]
      assert isWellDefined D
      assert(dd^D_2 == -dd^C_3)
    Text
      In order to shift the complex one step, and not change the differential, one
      can do the following.
    Example
      E = complex(C, Base => -1)
      assert isWellDefined E
      assert(dd^E_2 == dd^C_3)
    Text
      The shift operator is functorial, as illustrated below.
    Example
      C2 = freeResolution (S^1/(a^2, b^2, c^2, d^2))
      C3 = freeResolution (S^1/(a^2, b^3, c^4, d^5))
      f2 = extend(C, C2, map(C_0, C2_0, 1))
      f3 = extend(C2, C3, map(C2_0, C3_0, 1))
      assert((f2*f3)[1] == (f2[1]) * (f3[1]))
      assert(source(f2[1]) == C2[1])
      assert(target(f2[1]) == C[1])
   SeeAlso
     concentration
     (complex, Complex)
     (extend, Complex, Complex, Matrix)
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
     (cohomology,ZZ,ZZdFactorization)
   Headline
     homology or cohomology module of a ZZ/d-graded factorization
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
        (isHomogeneous, ZZdFactorization)
    Headline
         whether a complex is homogeneous
    Usage
         isHomogeneous C
    Inputs
         C:Complex
    Outputs
         :Boolean
             that is true when {\tt C} is a homogeneous (i.e. graded) complex
    Description
        Text
            A complex is homogeneous (graded) if the base ring is graded,
            all of the component objects are graded, and
            all the component maps are graded of degree zero.
        Example
            S = ZZ/101[a,b,c,d];
            I = minors(2, matrix{{a,b,c},{b,c,d}})
            C = freeResolution (S^1/I)
            isHomogeneous C
            J = minors(2, matrix{{a,b,c},{b,c,d^2}})
            D = freeResolution (S^1/J)
            isHomogeneous D
    SeeAlso
        "Basic invariants and properties"
        isHomogeneous
        (isHomogeneous, ComplexMap)
///

doc ///
   Key
     (symbol**, ZZdFactorization, ZZdFactorization)
     (symbol**, Complex, ZZdFactorization)
     (symbol**, ZZdFactorization, Complex)
     (symbol**, ZZdFactorization, Module)
     (symbol**, Module, ZZdFactorization)
     (tensor, ZZdFactorization, ZZdFactorization)
   Headline
     tensor product of ZZ/d-graded factorizations
   Usage
     D = C1 ** C2
   Inputs
     C1:ZZdFactorization
       or @ofClass Module@ or @ofClass Complex@
     C2:ZZdFactorization
       or @ofClass Module@ or @ofClass Complex@
   Outputs
     D:ZZdFactorization
       tensor product of {\tt C1} and {\tt C2}
   Description
    Text
      The tensor product is a complex $D$ whose $i$th component is
      the direct sum of $C1_j \otimes C2_k$ over all $i = j+k$.
      The differential on $C1_j \otimes C2_k$ is the differential 
      $dd^{C1} \otimes id_{C2} + (-1)^j id_{C1} \otimes dd^{C2}$.
      
      As the next example illustrates, the Koszul complex can be constructed via iterated tensor products.
    Example
      S = ZZ/101[a..c]
      Ca = complex {matrix{{a}}}
      Cb = complex {matrix{{b}}}
      Cc = complex {matrix{{c}}}
      Cab = Cb ** Ca
      dd^Cab
      assert isWellDefined Cab
      Cabc = Cc ** Cab
      Cc ** Cb ** Ca
      dd^Cabc
      assert isWellDefined Cabc
    Text
      If one of the arguments is a module, it is considered as a complex concentrated in homological degree 0.
    Example
      Cabc ** (S^1/(a,b,c))
      S^2 ** Cabc
    Text
      Because the tensor product can be regarded as the total complex of a double complex,
      each term of the tensor product comes with pairs of indices, labelling the summands.
    Example
      indices Cabc_1
      components Cabc_1
      Cabc_1_[{1,0}]
      indices Cabc_2
      components Cabc_2
      Cabc_2_[{0,2}]
   SeeAlso
     indices
     components
     directSum
///

doc ///
   Key
     (Hom, ZZdFactorization, ZZdFactorization)
     (Hom, ZZdFactorization, Module)
     (Hom, Module, ZZdFactorization)     
     (Hom, ZZdFactorization, Ring)
     (Hom, Ring, ZZdFactorization)     
   Headline
     the ZZ/d-graded homomorphism factorization between two ZZ/d-graded factorizations
   Usage
     D = Hom(C1,C2)
   Inputs
     C1:ZZdFactorization
       or @ofClass Module@, or @ofClass Ring@
     C2:ZZdFactorization
       or @ofClass Module@, or @ofClass Ring@
   Outputs
     D:ZZdFactorization
       the ZZ/d-graded factorization of homomorphisms between {\tt C1} and {\tt C2}
   Description
    Text
      The complex of homomorphisms is a complex $D$ whose $i$th component is
      the direct sum of $Hom(C1_j, C2_{j+i})$ over all $j$.
      The differential on $Hom(C1_j, C2_{j+i})$ is the differential 
      $Hom(id_{C1}, dd^{C2}) + (-1)^j Hom(dd^{C1}, id_{C2})$.
      $dd^{C1} \otimes id_{C2} + (-1)^j id_{C1} \otimes dd^{C2}$.

      In particular, for this operation to be well-defined, both
      arguments must have the same underlying ring.
    Example
      S = ZZ/101[a..c]
      C = freeResolution coker vars S
      D = Hom(C,C)
      dd^D
      assert isWellDefined D
    Text
      The homology of this complex is $Hom(C, ZZ/101)$
    Example
      prune HH D == Hom(C, coker vars S)
    Text
      If one of the arguments is a module or a ring, it is considered as a complex concentrated in homological degree 0.
    Example
      E = Hom(C, S^2)
      prune HH E
    Text
      There is a simple relationship between Hom complexes and @TO2 ((symbol SPACE, Complex, Array), "shifts")@.
      Specifically, shifting the first argument is the same as the negative shift of the result.  But
      shifting the second argument is only the same as the positive shift of the result
      up to a sign.
    Example
      Hom(C[3], C) == D[-3]
      Hom(C, C[-2]) == D[-2]
      Hom(C, C[-3]) != D[-3]
      Hom(C, C[-3]) == complex(- dd^(D[-3]))
    Text
      Specific maps and morphisms between complexes can be obtained
      with @TO (homomorphism, ComplexMap)@.
    Text
      Because the Hom complex can be regarded as the total complex of a double complex,
      each term comes with pairs of indices, labelling the summands.
    Example
      indices D_-1
      components D_-1
      indices D_-2
      components D_-2
   SeeAlso
     (homomorphism, ComplexMap)
     (homomorphism', ComplexMap)
     (randomComplexMap, Complex, Complex)
     indices
     components
     (Hom, ComplexMap, ComplexMap)
///

doc ///
    Key
        (homomorphism, ZZdFactorizationMap)
    Headline
        get the homomorphism from an element of Hom
    Usage
        g = homomorphism f
    Inputs
        f:ZZdFactorizationMap
            a map of the form $f : R^1 \to Hom(C, D)$, where
            $C$ and $D$ are complexes,
            $Hom(C,D)$ has been previously computed, and $R$ is
            the underlying ring of these complexes
    Outputs
        g:ZZdFactorizationMap
            the corresponding map of chain complexes from $C$ to $D$
    Description
        Text
            As a first example, consider two Koszul complexes $C$ and $D$.
            From a random map $f : R^1 \to Hom(C, D)$, we construct 
            the corresponding map of chain complexes $g : C \to D$.
        Example
            R = ZZ/101[a,b,c]
            C = freeResolution ideal"a,b,c"
            D = freeResolution ideal"a2,b2,c2"
            H = Hom(C,D)
            f = randomFactorizationMap(H, complex R^{-2})
            isWellDefined f
            g = homomorphism f
            isWellDefined g
            assert not isCommutative g
        Text
            The map $g : C \to D$ corresponding to a random map into $Hom(C,D)$
            does not generally commute with the differentials.  However, if the
            element of $Hom(C,D)$ is a cycle, then the corresponding map does commute.
        Example
            f = randomComplexMap(H, complex R^{-2}, Cycle => true)
            isWellDefined f
            g = homomorphism f
            isWellDefined g
            assert isCommutative g
            assert(degree g === 0)
            assert(source g === C)
            assert(target g === D)
            assert(homomorphism' g == f)
        Text
            A homomorphism of non-zero degree can be encoded
            in (at least) two ways.
        Example
            f1 = randomComplexMap(H, complex R^1, Degree => -2)
            f2 = map(target f1, (source f1)[2], i -> f1_(i+2))
            assert isWellDefined f2
            g1 = homomorphism f1
            g2 = homomorphism f2
            assert(g1 == g2)
            assert isWellDefined g1
            assert isWellDefined g2
            homomorphism' g1 == f1
            homomorphism' g2 == f1
    SeeAlso
        (homomorphism, Matrix)
        (homomorphism, ZZ, Matrix, Complex)
        (homomorphism', ComplexMap)
        (Hom, Complex, Complex)
        (randomComplexMap, Complex, Complex)
///

doc ///
    Key
        (homomorphism', ZZdFactorizationMap)
    Headline
        get the element of Hom from a map of ZZ/d-graded factorizations
    Usage
        f = homomorphism g
    Inputs
        g:ZZdFactorizationMap
            from $C$ to $D$
    Outputs
        f:ZZdFactorizationMap
            a map of the form $f : R^1 \to Hom(C, D)$, where
            $R$ is the underlying ring of these ZZ/-graded factorizations
    Description
        Text
            As a first example, consider two Koszul complexes $C$ and $D$.
            From a random map $f : R^1 \to Hom(C, D)$, we construct 
            the corresponding map of chain complexes $g : C \to D$.
        Example
            R = ZZ/101[a,b,c]
            C = freeResolution ideal"a,b,c"
            D = freeResolution ideal"a2,b2,c2"
            g = randomComplexMap(D, C, InternalDegree => 2)
            isWellDefined g
            f = homomorphism' g
            isWellDefined f
        Text
            The map $g : C \to D$ corresponding to a random map into $Hom(C,D)$
            does not generally commute with the differentials.  However, if the
            element of $Hom(C,D)$ is a cycle, then the corresponding map does commute.
        Example
            g = randomComplexMap(D, C, Cycle => true, InternalDegree => 3)
            isWellDefined g
            f = homomorphism' g
            isWellDefined f
            assert isCommutative g
            assert(degree f === 0)
            assert(source f == complex(R^{-3}))
            assert(target g === D)
            assert(homomorphism f == g)
    SeeAlso
        "Working with Ext"
        (homomorphism', Matrix)
        (homomorphism, ZZdFactorizationMap)
        (Hom, ZZdFactorization, ZZdFactorization)
        (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (homomorphism, ZZ, Matrix, ZZdFactorization)
    Headline
        get the homomorphism from an element of Hom
    Usage
        g = homomorphism(i, f, E)
    Inputs
        i:ZZ
        f:Matrix
            a map of the form $f \colon R^1 \to E_i$
        E:ZZdFactorization
            having the form
            $E = \operatorname{Hom}(C, D)$ for some ZZZ/d-graded factorizations $C$ and $D$
    Outputs
        g:ZZdFactorizationMap
            the corresponding map of chain complexes from $C$ to $D$ of degree $i$
    Description
        Text
            An element of the complex $\operatorname{Hom}(C, D)$ corresponds to a map of 
            ZZ/d-graded factorizations from $C$ to $D$.  Given an element in the $i$-th term, this
            method returns the corresponding map of factorizations of degree $i$.
        Text
            As a first example, consider two Koszul complexes $C$ and $D$.
            From a random map $f \colon R^1 \to Hom(C, D)$, we construct 
            the corresponding map of chain complexes $g \colon C \to D$.
        Example
            R = ZZ/101[a,b,c];
            C = freeResolution ideal"a,b,c"
            D = freeResolution ideal"a2,b2,c2"
            E = Hom(C,D)
            f = random(E_2, R^{-5})
            g = homomorphism(2, f, E)
            assert isWellDefined g
            assert not isCommutative g
        Text
            The map $g \colon C \to D$ corresponding to a random map into $Hom(C,D)$
            does not generally commute with the differentials.  However, if the
            element of $Hom(C,D)$ is a cycle, then the corresponding map does commute.
        Example
            h = randomComplexMap(E, complex R^{-2}, Cycle => true, Degree => -1)
            f = h_0
            g = homomorphism(-1, f, E)
            assert isWellDefined g
            assert isCommutative g
            assert(degree g === -1)
            assert(source g === C)
            assert(target g === D)
            assert(homomorphism' g == h)
    SeeAlso
        "Working with Ext"
        (homomorphism, Matrix)
        (homomorphism', ComplexMap)
        (Hom, Complex, Complex)
        (randomComplexMap, Complex, Complex)
///



-- TODO: once we have Hom evaluation map,
-- let's add in the map from C to dual dual C.
doc ///
   Key
     (dual, ZZdFactorization)
   Headline
     make the dual of a ZZ/d-graded factorization
   Usage
     dual C
   Inputs
     C:ZZdFactorization
   Outputs
     :ZZdFactorization
   Description
    Text
      The dual of a ZZ/d-graded factorization $C$ is by definition $Hom(C, R)$, where $R$ is the ring of $C$.
    Example
      S = ZZ/101[a..d];
      B = intersect(ideal(a,c),ideal(b,d))
      C1 = freeResolution B
      C2 = dual C1
      assert(C2 == Hom(C1, S^1))
      C1 == dual dual C1
      prune HH C2
    Text
      The double dual is not necessarily isomorphic to the original complex.
    Example
      I = ideal(a^2, a*b, b^2)
      J = ideal(b^3, b*c, c^3)
      K = intersect(I,J)
      f = map(S^1/I ++ S^1/J, S^1/K, {{1},{1}})
      g = map(S^1/(I+J), S^1/I ++ S^1/J, {{1,-1}})
      C = complex{g,f}
      assert isWellDefined C
      assert isExact C
      assert(dual C == 0)
      assert(C != dual dual C)
   SeeAlso
     (Hom, Complex, Complex)
     (dual, Module)
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
        (minimalPresentation, ZZdFactorization)
        (prune, ZZdFactorization)
        (prune, ZZdFactorizationMap)
        (minimalPresentation, ZZdFactorizationMap)
    Headline
        minimal presentation of all terms in a ZZ/d-graded factorization
    Usage
        D = minimalPresentation C
        D = prune C
        h = minimalPresentation f
        h = prune f
    Inputs
        C:ZZdFactorization
            or $f$ @ofClass ZZdFactorizationMap@
        Exclude => 
            unused
    Outputs
        D:ZZdFactorization
            isomorphic to the input, where each term is replaced
            by a minimally presented model, or $h$ @ofClass ZZdFactorizationMap@
            where the source and target are minimally presented
    Consequences
        Item
            The isomorphism $g : D \to C$ is available as 
            @TT "g = D.cache.pruningMap"@.  The inverse isomorphism
            can be obtained as @TT "g^-1"@
    Description
        Text
            This is frequently useful to make the output of certain
            operations readable or understandable.  This operation
            is functorial, applying both to ZZ/d-graded factoriations and factorization maps.
        Text
            In particular, homology often needs to be pruned to
            be understood.  For instance, this is useful 
            for recognizing when terms given by subquotient modules 
            are actually zero.
        Example
            S = ZZ/101[a,b,c,d,e];
            I = ideal(a,b) * ideal(c,d,e)
            F = dual freeResolution I
            C = HH F
            D = prune C
            g = D.cache.pruningMap
            assert isWellDefined g
            assert isComplexMorphism g
            assert (target g == C)
            assert (source g == D)
            g^-1
            assert(g*g^-1 == 1 and g^-1*g == 1)
        Text
            The image of a map of complexes also becomes more
            understandable via pruning.
        Example
            S = ZZ/101[a,b,c];
            I = ideal(a^2,b^2,c^2);
            J = I + ideal(a*b*c);
            FI = freeResolution I
            FJ = freeResolution J
            f = randomComplexMap(FJ, FI ** S^{-1}, Cycle => true)
            C = image f
            D = prune C
            g = D.cache.pruningMap
            assert isWellDefined g
            assert isComplexMorphism g
            assert (target g == C)
            assert (source g == D)
            g^-1
            assert(g*g^-1 == 1 and g^-1*g == 1)
        Text
            One can directly prune the map of complexes $f$.
        Example
            h = prune f
            assert(source h === prune source f)
            assert(target h === prune target f)
   SeeAlso
       "Making chain complexes"
       --(minimize, ZZdFactorization)
       (minimalPresentation, Module)
       randomFactorizationMap
       isComplexMorphism
///

-*doc ///
    Key
        (minimize, ZZdFactorization)
        minimize
        minimizingMap
    Headline
        a quasi-isomorphic ZZ/d-graded factorization whose terms have minimal rank
    Usage
        D = minimize C
    Inputs
        C:ZZdFactorization
            graded, whose terms are all free modules
    Outputs
        D:ZZdFactorization
            graded, whose terms are all free modules of minimal rank
    Consequences
        Item
            The projection morphism $g : C \to D$ is available as 
            @TT "g = D.cache.minimizingMap"@.  
    Description
        Text
            This method essentially removes all scalar units 
            from the matrices in the differential of $C$.
            
            We illustrate this in a simple example.
        Example
            S = ZZ/32003[a,b];
            I = ideal(a^2-b^2, a*b)
            C = freeResolution(I, FastNonminimal=>true)
            betti C
            D = minimize C
            assert(isWellDefined D and isHomogeneous D)
            betti D
            g = D.cache.minimizingMap
            assert isWellDefined g
            assert(isComplexMorphism g and isQuasiIsomorphism g)
            assert(source g == C)
            assert(target g == D)
            assert(coker g == 0)
        Text
            The minimal complex $D$ is a direct summand of the
            original complex $C$.  The natural inclusion
            of $D$ into $C$ can be constructed as follows.
        Example
            f = liftMapAlongQuasiIsomorphism(id_D, g)
            g*f == id_D
            assert(source f == D)
            assert(target f == C)
            assert(ker f == 0)
            f*g
        Text
            The chain complex $D$ is a direct summand of $C$,
            giving rise to a split short exact sequence of
            chain complexes.
        Example
            h = prune canonicalMap(C, ker g)
            assert isShortExactSequence(g, h)
        Text
            Warning: If the input complex is not homogeneous, then
            the output is probably not what one would expect.
        Example
            S = ZZ/32003[a..d]
            J = ideal(a*b*c-b*c, a*d-c, a^3-d^2*c)
            CJ = freeResolution J
            assert not isHomogeneous CJ
            D = minimize CJ
            isWellDefined D
            prune HH D == prune HH CJ
   SeeAlso
       freeResolution
       (resolution, Complex)
       (resolutionMap, Complex)
       (minimalPresentation, Complex)
///*-




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
 
-*doc ///
  Key
    (isFree, ZZdFactorization)
    isFree
  Headline
    whether a ZZ/d-graded factorization consists of free modules
  Usage
    isFree C
  Inputs
    C:ZZdFactorization
  Outputs
    :Boolean
      that is true when each $C_i$ is a free module
  Description
    Text
      This method checks whether the given representation of each
      module $C_i$ is free. To determine whether the complex $C$ is
      isomorphic to a free complex, use @TO2((prune,ZZdFactorization), "prune")@.
    Text
      The following example demonstrates that the presentation of a module
      might not reveal the property of being free.
    Example
      S = ZZ/101[a,b];
      M = kernel vars S
      assert not isFreeModule M
      assert isFreeModule prune M
    Text
      By definition, a free resolution $C$ consists of free modules.
      In contrast, the augmented complex $C'$ might or might not
      consist of free modules.
    Example
      C = freeResolution M
      assert isFree C
      C' = cone map(complex M, C, i -> map(M, C_0, 1))[1]
      isWellDefined C'
      assert not isFree C'
      prune C'
      assert isFree prune C'
  SeeAlso
    isFreeModule
    freeResolution
    (prune, Complex)
///*-






doc ///
    Key--CHANGE THIS TO KOSZUL FACTORIZATIONS
        (koszulComplex, Matrix)
        (koszulComplex, List)
        koszulComplex
        [(koszulComplex, Matrix), Concentration]
        [(koszulComplex, List), Concentration]
    Headline
        makes the Koszul factorization of a polynomial
    Usage
        K = koszulComplex f
    Inputs
        f:Matrix
            having one row, or a @ofClass List@ of ring elements
        Concentration => Sequence
            a pair {\tt (lo, hi)} which limits the non-zero terms in the output
    Outputs
        :Complex
            the Koszul complex (or a subcomplex)
    Description
        Text
            Let $R$ be a commutative ring and let $E$ be a free $R$-module of finite rank $r$.
            Given a linear map $f \colon E \to R$, the Koszul complex associated to $f$
            is the chain complex of $R$-modules

            $\phantom{WWWW}
              0 \leftarrow R \leftarrow \bigwedge^1 E \leftarrow \bigwedge^2 E \leftarrow \dotsb \leftarrow \bigwedge^r E 
              \leftarrow 0,
            $

            where the differential is given by 

            $\phantom{WWWW}
              dd_k(e_1 \wedge e_2 \wedge \dotsb \wedge e_k) = 
              \sum_{i=1}^k (-1)^{i+1} f(e_i) \, e_1 \wedge e_2 \wedge \dotsb \wedge \widehat{e_i} \wedge \dotsb \wedge e_k,
            $
            
            and the superscript hat means the term is omitted.  For this method, the linear map $f$ is
            given as either a matrix with one row, or a list of ring elements.
        Example
            S = QQ[a..d]
            koszulComplex {a}
            C = koszulComplex {a^2+b^2,c^3}
            dd^C
            K4 = koszulComplex vars S
            dd^K4
            assert isWellDefined K4
        Text
            To obtain natural subcomplexes, use the @TT "Concentration"@ option.
        Example
            koszulComplex(vars S, Concentration => (2,3))
            koszulComplex(vars S, Concentration => (-1,5))
        Text
            The koszul complex can be constructed as an iterated tensor product.
            The maps are identical, except that the even indexed differentials
            have the opposite sign.
        Example
            C = koszulComplex {d} ** (koszulComplex {c} ** (koszulComplex {b} ** koszulComplex {a}))
            K = koszulComplex {a,b,c,d}
            netList {{dd^C, dd^K}}
    SeeAlso
        "Making chain complexes"
        (symbol**, Complex, Complex)
///

///
    Key
        (trivialFactorization,RingElement)
    Headline
        Given a degree d monomial m, outputs the trivial d-fold factorization of m
    Usage
        trivialFactorization(m)
    Inputs
        m:RingElement
	    this element should be a monomial
    Outputs
        C:ZZdFactorization
	    the trivial factorization of m which decomposes m as a product of d linear forms
    Description
        Text
        Example
    Caveat
    SeeAlso
///


///
    Key
        (linearFactorization,RingElement)
	(linearFactorization,RingElement,RingElement)
	(linearFactorization,RingElement,Symbol)
    Headline
        Construct a linear d-fold factorization of a homogeneous degree d polynomial
    Usage
        linearFactorization(f)
	linearFactorization(f,t)
    Inputs
        f:RingElement
	    a degree d polynomial
	t:Symbol or RingElement
	    the name of the distinguished root of unity, only needed if f has degree d > 2
    Outputs
        C:ZZdFactorization
	    a ZZ/d-graded factorization of f with linear differentials, obtained by taking the
	    tensor product of the trivial factorizations of each of its monomial terms
    Description
        Text
        Example
    Caveat
    SeeAlso
///


///
    Key
        (randomFactorization,Ring)
	(randomFactorization,Ring,RingElement)
	(randomFactorization,Ring,Symbol)
    Headline
        Construct a linear factorization of a random degree d homogeneous polynomial
    Usage
        randomFactorization(d,Q)
	randomFactorization(d,Q,t)
    Inputs
        d:ZZ
	    specifies the degree of the desired polynomial
	Q:Ring
        t:Symbol or RingElement
	    the name of the distinguished root of unity, only needed if d > 2
    Outputs
        C:ZZdFactorization
	    A linear factorization of a randomly chosen degree d monomial of the input ring
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (higherHomotopyFactorization, List, Complex)
	(higherHomotopyFactorization, RingElement, Complex)
    Headline
        Construct the matrix factorization induced by a system of higher homotopies
    Usage
        higherHomotopyFactorization(L,C)
	higherHomotopyFactorization(f,C)
    Inputs
        C:Complex
        L:List
	    A list of ring elements annihilating the homology of C
	f:RingElement
	    A ring element annihilating the homology of C
    Outputs
        :ZZdFactorization
	    A matrix factorization of f_1t_1 + ... + f_n t_n if L = {f_1,...,f_n} or a matrix factorization of f
	    induced by a system of higher homotopies with respect to these polynomials.
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (toBranchedCover, ZZdFactorization, Symbol)
	(toBranchedCover, ZZdFactorization, RingElement)
    Headline
        Convert a ZZ/d-graded factorization of a ring element f into a module over the
	d-fold branched cover of f
    Usage
        toBranchedCover(C,z)
    Inputs
        C:ZZdFactorization
	z:Symbol or RingElement
	    this specifies the name of the variable to use in the branched cover
    Outputs
        M:Module
	    A maximal Cohen-Macaulay module over the d-fold branched cover hypersurface
    Description
        Text
        Example
    Caveat
        In order to guarantee that this function behaves well with respect to pushing forward, one should
	give the variables of the ambient ring degree 0.
    SeeAlso
///

///
    Key
        (zeroOutDegrees,Ring)
	(zeroOutDegrees,ZZdFactorization)
    Headline
        Redefine a ring or ZZ/d-graded factorization so that all variables have degree 0
    Usage
        zeroOutDegrees(R)
	zeroOutDegrees(C)
    Inputs
        R:Ring
	C:ZZdFactorization
    Outputs
        :Ring
	    A ring that is isomorphic to the input ring, but with the variables viewed as having degree 0
	:ZZdFactorization
	    A ZZ/d-graded factorization that is isomorphic to the original, but viewed as over a ring
	    where all variables have degree 0
    Description
        This is a helper function for converting ZZ/d-graded factorizations into maximal Cohen-Macaulay modules over d-fold branched coverings,
	since functions such as PushFwd work better when the variables of the original ring are viewed as having degree 0.
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (degreeSorter,ZZ,Module)
	(degreeSorter,ZZ,ZZ,Module)
	(degreeSorter,ZZ,Matrix)
	(degreeSorter,ZZ,ZZ,Matrix)
    Headline
        Reorder the homogeneous basis elements of a free module or the rank/source of a matrix in terms of
	the congruence classes of their degrees.
    Usage
        degreeSorter(d,M)
	degreeSorter(d,offset,M)
    Inputs
        d:ZZ
	    the integer determining the congruence classes that the degrees will be ordered with respect to
	offset:ZZ
	    an offset value
	M:Module or Matrix
    Outputs
        L:List
	    if the input M is a module, then the output is the list of positions, indicating an ordering of 
	    the homogeneous basis elements grouped according to the congruence class of their degrees modulo d
	N:Matrix
	     if the input M is a matrix, then the output is the sae matrix M but with rows/columns reordered
	     to respect the ordering of the basis elements according to the congruence class of their degrees modulo d
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (branchedToMF,Module)
	(branchedToMF,Module,Ring)
	(branchedToMF,Matrix)
	(branchedToMF,Matrix,Ring)
    Headline
        Convert a maximal Cohen-Macaulay module over a d-fold branched covering into a matrix over
	the base ring which yields a well-defined ZZ/d-graded factorization
    Usage
        branchedToMF(M)
	branchedToMF(M,Q)
    Inputs
        M:Module or Matrix
	    a maximal Cohen-Macaulay module over a d-fold branched cover, or a matrix whose cokernel
	    yields a maximal Cohen-Macaulay module over a d-fold branched cover
	Q:Ring
	    optional argument to specify the base ring of the d-fold branched covering
    Outputs
        N:Matrix
	    A matrix over the base ring that yields a well-defined ZZ/d-graded factorization
    Description
        Text
        Example
    Caveat
    SeeAlso
///


///
    Key
        (mooreMF,ZZ)
    Headline
        Create a matrix factorization induced by a Moore matrix and its adjoint
    Usage
        mooreMF(p)
    Inputs
        p:ZZ
	    p is the characteristic of the underlying field
    Outputs
        :ZZdFactorization
            A matrix factorization of the determinant of the generic Moore matrix
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (rk1MCM2gen,List,ZZ)
    Headline
        Construct every possible rank 1, 2-generated maximal Cohen-Macaulay module over the hypersurface --KELLER: CHECK WHICH HYPERSURFACE
    Usage
        rk1MCM2gen(L,d)
    Inputs
        L:List
	    any permutation of the set {2,3,4}
	d:ZZ
	    any integer specifying the characteristic of the ambient field
    Outputs
        :ZZdFactorization
	    A matrix factorization corresponding to a rank 1, 2-generated maximal Cohen-Macaulay module over --KEORERK
    Description
        Text
        Example
    Caveat
    SeeAlso
///

doc ///
    Key
        "Making maps between ZZ/d-graded factorizations"
    Headline
        information about the basic constructors
    Description
    	Text
    	    @SUBSECTION "Basic constructors"@
	Text
    	    @UL {
                TO (map, ZZdFactorization, ZZdFactorization, HashTable),
                TO (map, ZZdFactorization, ZZdFactorization, ZZ),
                TO (map, ZZdFactorization, ZZdFactorization, Function),
                TO (map, ZZdFactorization, ZZdFactorization, List),
                TO (map, ZZdFactorization, ZZdFactorization, ZZdFactorizationMap),
                TO (id, ZZdFactorization),
                TO "differential of a chain complex",
                TO (symbol SPACE, ZZdFactorizationMap, Array),
                TO (isWellDefined, ZZdFactorizationMap)
            }@
    	Text
    	    @SUBSECTION "Important computations creating new ZZ/d-graded factorization maps"@
        Text 
            @UL {
                TO (homology, ZZdFactorizationMap),
                TO (symbol**, ZZdFactorization, Matrix),
                TO (extend, ZZdFactorization, ZZdFactorization, Matrix),
                TO (nullHomotopy, ZZdFactorizationMap)
           }@
    	Text
    	    @SUBSECTION "Canonical maps between ZZ/d-graded factorizations"@
        Text
            Some ZZ/d-graded factorizations come with canonical maps.
            To access the ZZ/d-graded factorization map, 
            one uses @TO (canonicalMap, ZZdFactorization, ZZdFactorization)@.
            The following operations have associated canonical maps.
	Text
    	    @UL {
                TO (kernel, ZZdFactorizationMap),
                TO (cokernel, ZZdFactorizationMap),
                TO (image, ZZdFactorizationMap),
                TO (coimage, ZZdFactorizationMap),
                TO (cone, ZZdFactorizationMap),
                TO (cylinder, ZZdFactorizationMap),
                TO (inducedMap, ZZdFactorization, ZZdFactorization)
            }@
    	Text
    	    @SUBSECTION "Random maps of ZZ/d-graded factorizations"@
        Text
            The method @TO (randomFactorizationMap, ZZdFactorization, ZZdFactorization)@
            allows one to construct random ZZ/d-graded factorization maps,
            random morphisms between factorizations, and random
            null homotopies between factorizaions.
	Text
    	    @UL {
                TO (isCommutative, ZZdFactorizationMap),
                TO (isComplexMorphism, ComplexMap), --LOOK AT THIS ONE
                TO (isNullHomotopic, ZZdFactorizationMap)
            }@
    	Text
    	    @SUBSECTION "Elementary operations on complex maps"@
        Text
    	    @UL {
                TO "arithmetic with complex maps",
                TO (symbol +, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (symbol |, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (symbol ||, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (symbol ++, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (symbol **, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (Hom, ZZdFactorizationMap, ZZdFactorizationMap),
                TO (dual, ZZdFactorizationMap),
                TO (symbol _, ZZdFactorizationMap, Array),
                TO (symbol ^, ZZdFactorizationMap, Array),
                TO (part, List, ZZdFactorizationMap),
                TO (symbol SPACE, RingMap, ZZdFactorizationMap),
                TO (symbol **, RingMap, ZZdFactorizationMap)
            }@
    SeeAlso
        "Making chain complexes"
        "Basic invariants and properties"
        "Working with Ext"
        "Working with Tor"
        "Towards computing in the derived category"
///



doc ///
  Key
    ZZdFactorizationMap
  Headline
    the class of all maps between ZZ/d-graded factorizations
  Description
    Text
      @LITERAL ////<script> macros["\\Hom"] = "\\operatorname{Hom}" </script>////@

      A map of ZZ/d-graded factorizations $f \colon C \rightarrow D$ of degree $d$ is a
      sequence of maps $f_i \colon C_i \rightarrow D_{d+i}$.  
      No relationship between the maps $f_i$ and 
      and the differentials of either $C$ or $D$ is assumed.
      
      The set of all maps from $C$ to $D$ form
      the complex $\Hom(C,D)$ where $\Hom(C,D)_d$ consists of the
      maps of degree $d$.

      The usual algebraic operations are available: addition,
      subtraction, scalar multiplication, and composition. The
      identity map from a chain complex to itself can be produced with
      @TO "id"@. An attempt to add (subtract, or compare) a ring
      element to a chain complex will result in the ring element being
      multiplied by the appropriate identity map.
  SeeAlso
    Complex
///

doc ///
    Key
        (map, ZZdFactorization, ZZdFactorization, HashTable)
    Headline
        make a map of ZZ/d-graded factorizations
    Usage
        f = map(D, C, H)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        H:HashTable
            whose keys are integers, and whose values are the maps between
            the corresponding terms
        Degree => ZZ
            the degree of the resulting map
        DegreeLift => 
            unused
        DegreeMap =>
            unused
    Outputs
        f:ZZdFactorizationMap
    Description
        Text
            A map of ZZ/d-graded factorizations $f : C \rightarrow D$ of degree $d$ is a
            sequence of maps $f_i : C_i \rightarrow D_{d+i}$.  
            No relationship between the maps $f_i$ and 
            and the differentials of either $C$ or $D$ is assumed.
            
            We construct a map of chain complexes by specifying the
            individual maps between the terms.
        Example
            R = ZZ/101[a,b,c];
            C = freeResolution coker matrix{{a^2-b^2,b^3-c^3,c^4}}
            D = freeResolution coker vars R
            H = hashTable { 0 => map(D_0, C_0, 1),
                1 => map(D_1, C_1, {{a, 0, 0}, {-b, b^2, 0}, {0, -c^2, c^3}}),
                2 => map(D_2, C_2, {{a*b^2, 0, 0}, {-a*c^2, a*c^3, 0}, {b*c^2, -b*c^3, b^2*c^3}}),
                3 => map(D_3, C_3, {{a*b^2*c^3}})
                }
            f = map(D, C, H)
            assert isWellDefined f
            assert isHomogeneous f
            assert(degree f == 0)
            assert isComplexMorphism f
        Text
            The keys in the hash table index the terms in the source of the
            map.  If a key is missing, that map is taken to be the zero map.
            We illustrate by constructing a map of chain complexes
            having nonzero degree, and omitting one key in the hash table.
        Example
            E = D[-3]
            H = hashTable { 0 => map(E_3, C_0, 1),
                1 => map(E_4, C_1, {{a, 0, 0}, {-b, b^2, 0}, {0, -c^2, c^3}}),
                3 => map(E_6, C_3, {{a*b^2*c^3}})
                }
            g = map(E, C, H, Degree => 3)
            g_2
            assert(g_1 == f_1)
            assert(g != f)
            assert isWellDefined g
            assert isHomogeneous g
            assert(degree g == 3)
            assert not isComplexMorphism g
            assert not isCommutative g
            assert(source g == C)
            assert(target g == E)
        Text
            This is the primary constructor used by all of the more user friendly
            methods for constructing a chain complex.
    Caveat
        This constructor minimizes computation
        and does very little error checking. To verify that a complex
        is well constructed, use @TO (isWellDefined, ComplexMap)@.
    SeeAlso
        ZZdFactorizationMap
        (map, ZZdFactorization, ZZdFactorization, Function)
        (isWellDefined, ZZdFactorizationMap)
        (isHomogeneous, ZZdFactorizationMap)
        (degree, ZZdFactorizationMap)
        (isComplexMorphism, ComplexMap)
        (isCommutative, ZZdFactorizationMap)
        (source, ZZdFactorizationMap)
        (target, ZZdFactorizationMap)
///

doc ///
    Key
        (map, ZZdFactorization, ZZdFactorization, List)
    Headline
        make a map of ZZ/d-graded factorizations
    Usage
        f = map(D, C, L)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        L:List
            consisting of either matrices, or lists of maps of factorizations
        Degree => ZZ
            the degree of the resulting map
        DegreeLift => 
            unused
        DegreeMap =>
            unused
    Outputs
        f:ZZdFactorizationMap
            from $C$ to $D$
    Description
        Text
            A map of complexes $f \colon C \rightarrow D$ of degree $d$ is a
            sequence of maps $f_i \colon C_i \rightarrow D_{d+i}$.  
            No relationship between the maps $f_i$ and 
            and the differentials of either $C$ or $D$ is assumed.
            
            This method has two very different usages.  The first is to 
            construct a chain complex map from a list of matrices.  The second
            constructs a chain complex map from essentially a block matrix
            whose entries are chain complex maps.
        Text
            In the first case, we construct a map of chain complexes
            by specifying the individual maps between the terms.
        Example
            R = ZZ/101[a,b,c];
            C = freeResolution coker matrix{{a^2-b^2,b^3-c^3,c^4}}
            D = freeResolution coker vars R
            L = {map(D_0, C_0, 1),
                map(D_1, C_1, {{a, 0, 0}, {-b, b^2, 0}, {0, -c^2, c^3}}),
                map(D_2, C_2, {{a*b^2, 0, 0}, {-a*c^2, a*c^3, 0}, {b*c^2, -b*c^3, b^2*c^3}}),
                map(D_3, C_3, {{a*b^2*c^3}})
                }
            f = map(D, C, L)
            assert isWellDefined f
            assert isHomogeneous f
            assert(degree f == 0)
            assert isComplexMorphism f
        Text
            In the second, we construct a map of chain complexes via a block matrix
            whose individual entries are already maps of chain complexes.
            We illustrate by constructing a mapping cone.
        Example
            f = extend(D,C,id_(R^1))
            assert(degree f == 0)
            g = map(D, C[-1], f[-1], Degree => -1) -- a variant of f having degree -1
            cf = map(E = C[-1] ++ D, E, {
                    {dd^(C[-1]),    0}, 
                    {         g, dd^D}
                    })
            assert isWellDefined cf
            assert(degree cf == -1)
        Text
            We convert this map of complexes {\tt cf} into the differential of the mapping cone.
            For the following constructor, the source and target of
            the input must be identical, in this case the chain complex $E$.
        Example
            conef = complex cf 
            assert isWellDefined conef
            assert(conef == cone f)
    SeeAlso
        "Making maps between chain complexes"
        (map, Complex, Complex, HashTable)
        (degree, ComplexMap)
        (extend, Complex, Complex, Matrix)
        (cone, ComplexMap)
///

doc ///
    Key
        (map, ZZdFactorization, ZZdFactorization, ZZ)
    Headline
        make the zero map or identity between ZZ/d-graded factorizations
    Usage
        f = map(D, C, 0)
        f = map(C, C, 1)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        0:ZZ
            or 1
        Degree => ZZ
            the degree of the resulting map
        DegreeLift => 
            unused
        DegreeMap =>
            unused
    Outputs
        f:ZZdFactorizationMap
            the zero map from $C$ to $D$ or the identity map from $C$ to $C$
    Description
        Text
            A map of ZZ/d-graded factorizations $f : C \rightarrow D$ of degree $d$ is a
            sequence of maps $f_i : C_i \rightarrow D_{d+i}$.  
            
            We construct the zero map between two
            factorizations.
        Example
            R = QQ[a,b,c]
            C = freeResolution coker vars R
            D = freeResolution coker matrix{{a^2, b^2, c^2}}
            f = map(D, C, 0)
            assert isWellDefined f
            assert isComplexMorphism f
            g = map(C, C, 0, Degree => 13)
            assert isWellDefined g
            assert(degree g == 13)
            assert not isComplexMorphism g
            assert isCommutative g
            assert isHomogeneous g
            assert(source g == C)
            assert(target g == C)
        Text
            Using this function to create the identity map
            is the same as using @TO (id, Complex)@.
        Example
            assert(map(C, C, 1) === id_C)
   SeeAlso
        ZZdFactorizationMap
        (map, ZZdFactorization, ZZdFactorization, Function)
        (isWellDefined, ZZdFactorizationMap)
        (isHomogeneous, ZZdFactorizationMap)
        (degree, ZZdFactorizationMap)
        (isCommutative, ZZdFactorizationMap)
        (source, ZZdFactorizationMap)
        (target, ZZdFactorizationMap)
        (id, ZZdFactorization)
///

doc ///
    Key
        (map, ZZdFactorization, ZZdFactorization, Function)
    Headline
        make a map of ZZ/d-graded factorizations
    Usage
        f = map(D, C, fcn)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        fcn:Function
            whose values at integers are the maps between
            the corresponding terms
        Degree => ZZ
            the degree of the resulting map
        DegreeLift => 
            unused
        DegreeMap =>
            unused
    Outputs
        f:ZZdFactorizationMap
    Description
        Text
            A map of complexes $f : C \rightarrow D$ of degree $d$ is a
            sequence of maps $f_i : C_i \rightarrow D_{d+i}$.  
            No relationship between the maps $f_i$ and 
            and the differentials of either $C$ or $D$ is assumed.
            
            We construct a map of chain complexes by specifying a
            function which determines the maps between the terms.
        Example
            R = ZZ/101[x]/x^3;
            M = coker vars R
            C = freeResolution(M, LengthLimit => 6)
            D = C[1]
            f = map(D, C, i -> 
                if odd i then 
                    map(D_i, C_i, {{x}})
                else map(D_i, C_i, {{x^2}})
                )
            assert isWellDefined f
            assert isCommutative f
            assert(source f == C)
            assert(target f == D)
    SeeAlso
        ComplexMap
        (isWellDefined, ComplexMap)
        (isCommutative, ComplexMap)
        (source, ComplexMap)
        (target, ComplexMap)
///


doc ///
    Key
        (map, ZZdFactorization, ZZdFactorization, ZZdFactorizationMap) --IS THIS IMPLEMENTED??
    Headline
        make a new map of chain complexes from an existing one
    Usage
        g = map(D, C, f)
    Inputs
        C:Complex
        D:Complex
        f:ComplexMap
            regarded as providing matrices which induce maps between the terms of $C$ and $D$
        Degree => ZZ
            the degree $d$ of the resulting map
        DegreeLift => 
            unused
        DegreeMap =>
            unused
    Outputs
        g:ComplexMap
    Description
        Text
            A map of complexes $f : C' \rightarrow D'$ is a
            sequence of maps $f_i : C'_i \rightarrow D'_{d'+i}$.  
            The new map $g : C \rightarrow D$ is the sequence of maps $g_i : C_i \rightarrow D_{d+i}$
            induced by the matrix of $f_i$.
            
            One use for this function is to get the new map of chain complexes obtained by shifting 
            the source or target of an existing chain map.  For example, one can regard the differential
            on a complex can be regarded as a map of degree zero between shifted complexes.
        Example
            R = ZZ/101[a,b,c];
            C = freeResolution coker vars R
            f = map(C[-1], C, dd^C, Degree => 0)
            assert isWellDefined f
            assert(degree f == 0)
            assert isCommutative f
            assert isComplexMorphism f
            assert not isComplexMorphism dd^C
    SeeAlso
        ComplexMap
        (map, Complex, Complex, Function)
        (isWellDefined, ComplexMap)
        (degree, ComplexMap)
        (isComplexMorphism, ComplexMap)
        (isCommutative, ComplexMap)
        (symbol SPACE, Complex, Array)
///

doc ///
    Key
        (id, ZZdFactorization)
    Headline
        the identity map of a ZZ/d-graded factorization
    Usage
        f = id_C
    Inputs
        C:ZZdFactorization
    Outputs
        f:ZZdFactorizationMap
          the identity map from $C$ to itself
    Description
        Text
            The collection of ZZ/d-graded factorizations together with ZZ/d-graded factorization morphisms
            forms a category.  In particular, every ZZ/d-graded factorization has an identity map.
        Example
            R = ZZ/101[x,y]/(x^3, y^3)
            C = freeResolution(coker vars R, LengthLimit=>6)
            f = id_C
            assert isWellDefined f
            assert isComplexMorphism f
        Text
            The identity map corresponds to an element of
            the Hom complex.
        Example
            R = ZZ/101[a,b,c]
            I = ideal(a^2, b^2, b*c, c^3)
            C = freeResolution I
            D = Hom(C, C)
            homomorphism' id_C
    SeeAlso
        (map, Complex, Complex, ZZ)
        (isWellDefined, ComplexMap)
        (isComplexMorphism, ComplexMap)
        (Hom, Complex, Complex)
        (homomorphism', ComplexMap)
///

doc /// 
    Key
        (isWellDefined, ZZdFactorizationMap)
    Headline
        whether a map of ZZ/d-graded factorizations is well-defined
    Usage
        isWellDefined f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        :Boolean
            that is true when {\tt f} determines a well defined ZZ/d-graded factorization map
    Description
        Text
            A map of chain complexes $f : C \to D$ of degree $d$ is a sequence of
            maps $f_i : C_i \to D_{d+i}$.  No relationship is required between
            these maps and the differentials in the source and target.

            This routine checks that $C$ and $D$ are well-defined
            chain complexes, and that, for each $f_i$, the source and
            target equal $C_i$ and $D_{d+i}$, respectively.  If the
            variable {\tt debugLevel} is set to a value greater than
            zero, then information about the nature of any failure is
            displayed.
        Text
            Unlike the @TO2((isWellDefined, Complex), 
                "corresponding function for Complexes")@,
            the basic constructors for complex maps are all but
            assured to be well defined. The only case that could cause
            a problem is if one constructs the source or target
            complex, and those are not well defined.
        Example
            R = ZZ/101[a,b,c];
            C = freeResolution coker matrix{{a^2-b^2,b^3-c^3,c^4}}
            D = freeResolution coker vars R
            H = hashTable { 0 => map(D_0, C_0, 1),
                1 => map(D_1, C_1, {{a, 0, 0}, {-b, b^2, 0}, {0, -c^2, c^3}}),
                2 => map(D_2, C_2, {{a*b^2, 0, 0}, {-a*c^2, a*c^3, 0}, {b*c^2, -b*c^3, b^2*c^3}}),
                3 => map(D_3, C_3, {{a*b^2*c^3}})
                }
            f = map(D, C, H)
            assert isWellDefined f
            assert isHomogeneous f
            assert(degree f == 0)
            assert isComplexMorphism f
        Text
            We construct two random maps of chain complexes,
            and check to see that, as should be the case, 
            both are well defined.
        Example
            g = randomComplexMap(D,C)
            assert isWellDefined g
            assert not isCommutative g
        Example
            h = randomComplexMap(D,C, Cycle => true)
            assert isWellDefined h
            assert isComplexMorphism h
        Text
            This method also checks the following aspects of 
            the data structure:
        Text
            @UL {
                TEX "The underlying hash table has exactly the expected keys,
                namely, {\\tt source, target, degree, map, cache}",
                "The ring of the source and target are the same",
                "The source and target are well defined complexes",
                "The degree is an integer",
                TEX "All keys in the {\\tt map} field are integers,
                in the range of the concentration of the source",
                TEX "The source and target of each $f_i$ is as expected",
                TEX "If the {\\tt isCommutative} key is present in the cache
                table, then commutativity of the map with the differentials
                is checked"
                }@
    SeeAlso
        (isWellDefined, ZZdFactorization)
        (isCommutative, ZZdFactorizationMap)
        (map, ZZdFactorization, ZZdFactorization, HashTable)
///

doc ///
    Key
        (source, ZZdFactorizationMap)
    Headline
        get the source of a map of ZZ/d-graded factorizations
    Usage
        C = source f
    Inputs
      f:ZZdFactorizationMap
    Outputs
      C:ZZdFactorization
    Description
        Text
            Given a ZZ/d-graded factorization map $f : C \to D$
            this method returns the ZZ/d-graded factorization $C$.
        Example
            R = ZZ/101[a..d]
            I = ideal(a^2, b^2, c^2)
            J = I + ideal(a*b*c)
            FI = freeResolution I
            FJ = freeResolution J
            f = randomComplexMap(FJ, FI, Cycle=>true)
            source f
            assert isWellDefined f
            assert isComplexMorphism f
            assert(source f == FI)
            assert(target f == FJ)
        Text
            The differential in a complex is a map of chain complexes.
        Example
            kk = coker vars R
            F = freeResolution kk
            source dd^F == F
            target dd^F == F
            degree dd^F == -1
   SeeAlso
       "Making chain complexes"
       (target, ZZdFactorizationMap)
       (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (target, ZZdFactorizationMap)
    Headline
        get the target of a map of ZZ/d-graded factorizations
    Usage
        C = target f
    Inputs
      f:ZZdFactorizationMap
    Outputs
      C:ZZdFactorization
    Description
        Text
            Given a ZZ/d-graded factorization map $f : C \to D$
            this method returns the ZZ/d-graded factorization $D$.
        Example
            R = ZZ/101[a..d]
            I = ideal(a^2, b^2, c^2)
            J = I + ideal(a*b*c)
            FI = freeResolution I
            FJ = freeResolution J
            f = randomComplexMap(FJ, FI, Cycle=>true)
            target f
            assert isWellDefined f
            assert isComplexMorphism f
            assert(target f == FJ)
            assert(source f == FI)
        Text
            The differential in a complex is a map of chain complexes.
        Example
            kk = coker vars R
            F = freeResolution kk
            target dd^F == F
            source dd^F == F
            degree dd^F == -1
   SeeAlso
       "Making chain complexes"
       (source, ZZdFactorizationMap)
       (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (degree, ZZdFactorizationMap)
    Headline
        get the degree of a map of ZZ/d-graded factorizations
    Usage
        degree f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        :ZZ
    Description
        Text
            A ZZ/d-graded factorization map $f : C \to D$ of degree $d$ is a sequence of
            of maps $f_i : C_i \to D_{i+d}$.
            This method returns $d$.
        Text
            The degree of the differential of a complex is always -1.
        Example
            R = ZZ/101[a..d];
            I = ideal(a^2, b^2, c^2)
            FI = freeResolution I
            assert(degree dd^FI == -1)
        Example
            FJ = freeResolution (I + ideal(a*b*c))
            f = randomComplexMap(FJ, FI, Cycle=>true, Degree => -2)
            assert(degree f == -2)
   SeeAlso
       "Basic invariants and properties"
       (source, ComplexMap)
       (target, ComplexMap)
       (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
///



doc ///
    Key
        (symbol _, ZZdFactorizationMap, ZZ)
    Headline
        access individual matrices in a ZZ/d-graded factorization map
    Usage
        f_i
    Inputs
        f:ZZdFactorizationMap
        i:ZZ
            the homological index
    Outputs
        :Matrix
            the {\tt i}-th map
    Description
        Text
            A ZZ/d-graded factorization map $f : C \to D$ of degree $d$ is a sequence of maps $f_i : C_i \to D_{i+d}$.
            This method allows one to access the individual $f_i$.
        Example
            S = ZZ/101[a..c];
            C = freeResolution coker matrix{{a^2, b^2, c^2}}
            D = freeResolution coker vars S
            f = randomComplexMap(D, C)
            f_2
            f_0
        Text
            Indices that are outside of the concentration are automatically zero.
        Example
            concentration f
            f_-1
            f_3
            f_4
    SeeAlso
        (symbol_, ZZdFactorization, ZZ)
        (period, ZZdFactorizationMap)
///



doc ///
    Key
        (components, ZZdFactorizationMap)
    Headline
        list the components of a direct sum
    Usage
        components f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        :List
            the component maps of a direct sum of maps of ZZ/d-graded factorizations
    Description
        Text
            A map of ZZ/d-graded factorizations stores its component maps.
        Example
            S = ZZ/101[a,b,c];
            C = freeResolution coker vars S
            g1 = id_C
            g2 = randomComplexMap(C[1], C[2], Boundary => true)
            f = g1 ++ g2
            assert isWellDefined f
            L = components f
            L_0 === g1
            L_1 === g2
            indices f
            f' = (greg => g1) ++ (mike => g2)
            components f'
            indices f'
        Text
            The names of the components are called indices, and are
            used to access the relevant inclusion and projection maps.
        Example
            f'_[mike]
            f'^[greg]
            f^[0]
            f_[0]
    SeeAlso
        (directSum, ZZdFactorizationMap)
        (components, ZZdFactorization)
        indices
        (symbol_, ZZdFactorizationMap, Array)
        (symbol^, ZZdFactorizationMap, Array)
///

doc ///
  Key
    (symbol*, ZZdFactorizationMap, ZZdFactorizationMap)
  Headline
    composition of homomorphisms of ZZ/d-graded factorizations
  Usage
    f = h * g
  Inputs
    h:ZZdFactorizationMap
      if a ring element or integer, then we multiply the ring element
      by the appropriate identity map
    g:ZZdFactorizationMap
  Outputs
    f:ZZdFactorizationMap
      the composition of $g$ followed by $h$
  Description
    Text
      If $g_i : C_i \rightarrow D_{d+i}$, and $h_j : D_j \rightarrow E_{e+j}$,
      then the composition corresponds to 
      $f_i := h_{d+i} * g_i : C_i \rightarrow E_{i+d+e}$.  In particular,
      the degree of the composition $f$ is the sum of the degrees of
      $g$ and $h$.
    Example
      R = ZZ/101[a..d]
      C = freeResolution coker vars R
      3 * dd^C
      0 * dd^C
      dd^C * dd^C
  SeeAlso
      "Making maps between chain complexes"
      "arithmetic with complex maps"
///

doc ///
    Key
        (symbol ^, ZZdFactorizationMap, ZZ)
    Headline
        the n-fold composition
    Usage
        f^n
    Inputs
        f:ZZdFactorizationMap
            whose source and target are the same ZZ/d-graded factorization
        n:ZZ
    Outputs
        :ZZdFactorizationMap
            the composition of $f$ with itself $n$ times.
    Description
        Text
            A ZZ/d-graded factorization map $f : C \to C$ can be composed with itself.
            This method produces these new maps of ZZ/d-graded factorizations.
        Text
            The differential on a ZZ/d-graded factorization should compose with itself d times to give a  
            scalar multiple of the identity map.
        Example
            S = ZZ/101[a..c];
            C = freeResolution coker matrix{{a^2, b^2, c^2}}
            f = dd^C
            f^2
            assert(source f == target f)
            assert(degree f == -1)
            assert(degree f^2 == -2)
        Example
            g = randomComplexMap(C, C, Degree => -1)
            g^2
            g^3
            assert(g^4 == 0)
        Text
            The zero-th power returns the identity map
        Example
            f^0 == id_C
            g^0 == id_C
        Text
            When $n$ is negative, the result is the $n$-fold power
            of the inverse complex map, if it exists.
        Example
            h = randomComplexMap(C, C)
            h^-1
            assert(h * h^-1 == id_C)
            h^-4
            assert(h^-4 * h^4 == id_C)
    SeeAlso
        (symbol^, Matrix, ZZ)
        (symbol^, ZZdFactorization, ZZ)
///

doc ///
   Key
     (symbol ==, ZZdFactorizationMap, ZZdFactorizationMap)
     (symbol ==, ZZdFactorizationMap, ZZ)
     (symbol ==, ZZ, ZZdFactorizationMap)
   Headline
     whether two ZZ/d-graded factorization maps are equal
   Usage
     f == g
     f == 0
     f == 1
   Inputs
     f:ZZdFactorizationMap
       or 0, or 1.
     g:ZZdFactorizationMap
       or 0, or 1.
   Outputs
     :Boolean
       that is true when {\tt f} and {\tt g} are equal
   Description
     Text
       Two ZZ/d-graded factorization maps are equal if they have the same source,
       the same target, and $f_i = g_i$ for all $i$.
     Example
       S = ZZ/101[a..c]
       C = freeResolution coker vars S
       f = id_C
       assert(f == 1)
       f === id_C[-1][1]
       f == id_C[-1][1]
     Text
       A complex map is equal to zero if all the maps are zero.
       This could require computation to determine if something that
       is superficially not zero is in fact zero.
     Example
       assert(0 * id_C == 0)
     Example
       g = randomComplexMap(C, C)
       h = canonicalMap(coker g, target g)
       assert(h == 0)
     Text
       Testing whether a map is equal to 1 is a shorthand for determining
       if the complex map is the identity.
       Although the matrices may appear to be the identity, the map is not the
       identity when the source and target are not equal.
     Example
       g = randomComplexMap(C, C, InternalDegree=>1, Cycle=>true)
       h = canonicalMap(coker g, target g)
       assert(h != 1)
     Text
       Testing for equality is not the same testing for isomorphism.
       In particular, different presentations of a complex need not be equal.
     Example
       D = prune image g
       p = D.cache.pruningMap
       p == 1
       assert(coker p == 0 and ker p == 0)
       assert(prune p == 1)
   SeeAlso
     (symbol ==, Complex, Complex)
     (symbol SPACE, ComplexMap, Array)
     randomComplexMap
     canonicalMap
     (prune, Complex)
///

doc ///
  Key
    (isCommutative, ZZdFactorizationMap)
  Headline
    whether a complex map commutes with the differentials
  Usage
    isCommutative f
  Inputs
    f:ZZdFactorizationMap
  Outputs
    :Boolean
      that is true when $f$ commutes with the differentials
  Description
    Text
      For a ZZ/d-graded factorization map $f : C \to D$ of degree $d$, this method
      checks whether, for all $i$, we have
      $dd^D_{i+d} * f_i = (-1)^d * (f_{i-1} * dd^C_i)$.
    Text
      We first construct a random factorization map which commutes with the differential.
    Example
      S = ZZ/101[a,b,c];
      C = freeResolution coker vars S
      D = C ** C
      f1 = randomComplexMap(D, C, Boundary => true, InternalDegree => 1)
      isCommutative f1
      assert(degree f1 == 0)
      assert isNullHomotopic f1
      assert(source f1 == C and target f1 == D)
    Text
      We next generate a complex map that is commutative and (likely) 
      induces a nontrivial map on homology.
    Example
      f2 = randomComplexMap(D, C, Cycle => true)
      isCommutative f2
      assert(degree f2 == 0)
      assert isComplexMorphism f2
    Text
      When the degree of the complex map is odd, isCommutative determines
      whether the map is anti-commutative.  We illustrate
      this for one square.
    Example
      f3 = randomComplexMap(D, C, Cycle => true, Degree=>1, InternalDegree => 1)
      isCommutative f3
      assert(degree f3 == 1)
      part1 = dd^D_3 * f3_2
      part2 = f3_1 * dd^C_2
      assert(part1 + part2 == 0)
    Text
      If the @TO "debugLevel"@ is greater than zero, then
      the location of the first failure of commutativity is displayed.
    Example
      f4 = randomComplexMap(D, C)
      isCommutative f4
      debugLevel = 1
      isCommutative f4
  SeeAlso
    isComplexMorphism
    randomComplexMap
    freeResolution
///

doc ///
  Key
    (isComplexMorphism, ComplexMap) --IS THIS IMPLEMENTED??
    isComplexMorphism
  Headline
    whether a complex map is a morphism of complexes
  Usage
    isComplexMorphism f
  Inputs
    f:ComplexMap
  Outputs
    :Boolean
      that is true when $f$ commutes with the differentials and has degree $0$
  Description
    Text
      For a complex map $f : C \to D$ of degree $d$, this method
      checks whether $d = 0$ and, for all $i$, we have
      $dd^D_{i+d} * f_i = (-1)^d * (f_{i-1} * dd^C_i)$.
    Text
      We first construct a random complex morphism.
    Example
      S = ZZ/101[a,b,c];
      C = freeResolution coker vars S
      D = C ** C
      f1 = randomComplexMap(D, C, Boundary => true, InternalDegree => 1)
      isComplexMorphism f1
      assert(degree f1 == 0)
      assert isNullHomotopic f1
      assert(source f1 == C and target f1 == D)
    Text
      We next generate a complex morphism that (likely) 
      induces a nontrivial map on homology.
    Example
      f2 = randomComplexMap(D, C, Cycle => true)
      isComplexMorphism f2
      assert(degree f2 == 0)
      assert isComplexMorphism f2
    Text
      When the degree is non-zero, the map is not a complex morphism.
      If the @TO "debugLevel"@ is greater than zero, then
      information about the failure is displayed.
    Example
      f3 = randomComplexMap(D, C, Cycle => true, Degree=>1, InternalDegree => 1)
      assert(degree f3 == 1)
      isComplexMorphism f3
      debugLevel = 1
      isComplexMorphism f3
      assert isCommutative f3
    Example
      f4 = randomComplexMap(D, C)
      assert(degree f4 == 0)
      debugLevel = 0
      isComplexMorphism f4
      debugLevel = 1
      isComplexMorphism f4
  SeeAlso
    (isCommutative, ComplexMap)
    randomComplexMap
    freeResolution
///

doc ///
    Key
        (Hom, ZZdFactorizationMap, ZZdFactorizationMap)
        (Hom, ZZdFactorization, ZZdFactorizationMap)
        (Hom, ZZdFactorization, Matrix)
        (Hom, ZZdFactorizationMap, Module)
        (Hom, ZZdFactorizationMap, Matrix)
        (Hom, ZZdFactorizationMap, ZZdFactorization)
        (Hom, ZZdFactorizationMap, Ring)
        (Hom, Matrix, ZZdFactorization)
        (Hom, Matrix, ZZdFactorizationMap)
        (Hom, Module, ZZdFactorizationMap)
        (Hom, Ring, ZZdFactorizationMap)
	(Hom, ZZdFactorizationMap, ZZdFactorizationMap,RingElement)
	(Hom, ZZdFactorizationMap, ZZdFactorizationMap,Symbol)
    Headline
        the map of factorizations between Hom complexes
    Usage
        h = Hom(f,g)
	h = Hom(f,g,t)
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
	t:RingElement or Symbol
	     optional input specifying the symbol used to represent the distinguished root of unity, 
	     required only if the underlying ring does not have a specified root of unity and the period is > 2
    Outputs
        h:ZZdFactorizationMap
    Description
        Text
            The maps $f : C \to D$ and $g : E \to F$ of ZZ/d-graded factorizations induces the map
            $h = Hom(f,g) : Hom(D,E) \to Hom(C,F)$ defined by $\phi \mapsto g \phi f$.
        Example
            S = ZZ/101[a..c];
            C = freeResolution coker vars S
            D = (freeResolution coker matrix{{a^2,a*b,b^3}})[-1]
            f = randomComplexMap(D,C)
            E = (dual C)[-3]
            F = (dual D)[-3]
            g = randomComplexMap(F,E)
            h = Hom(f,g)
            assert isWellDefined h
            assert(source h === Hom(D,E))
            assert(target h === Hom(C,F))
        Text
            We illustrate the defining property of the map $h$ on a random element $\phi$
            in degree zero.
        Example
            e = randomComplexMap(source h, complex(S^1))
            phi = homomorphism e
            psi = homomorphism'(g * phi * f)
            assert(h*e == psi)
        Text
            If either of the arguments is a @TO "Complex"@, that argument is
            understood to be the identity map on that complex.
        Example
            assert(Hom(f, C) == Hom(f, id_C))
            assert(Hom(C, f) == Hom(id_C, f))
        Text
            If either of the arguments is a @TO "Module"@ or a @TO "Ring"@, that argument is
            understood to be the identity map on the complex having a unique non-zero term in 
            in homological degree 0.  The ring must be the underlying ring of the map of complexes.
        Example
            assert(Hom(f, S) == Hom(f, id_(complex S)))
            assert(Hom(S, f) == Hom(id_(complex S), f))
            M = S^1/(a^2, b^2, c^2)
            assert(Hom(f, M) == Hom(f, id _ (complex M)))
            assert(Hom(M, f) == Hom(id _ (complex M), f))
        Text
            If either of the arguments is a @TO "Matrix"@, that argument is
            understood to be a map of complexes whose source and target have a unique non-zero entry
            in homological degree 0.
        Example
            m = vars S;
            h1 = Hom(f, m)
            assert(h1 == Hom(f, map(complex target m, complex source m, i -> m)))
            m = vars S;
            h2 = Hom(m, f)
            assert(h2 == Hom(map(complex target m, complex source m, i -> m), f))
        Text
            XXX write this text after writing doc for homomorphism and homomorphism'.
        Example
            e = randomComplexMap(source h, complex(S^1, Base => -1))
            phi = homomorphism e
            assert(degree phi == -1)
            psi = homomorphism'(g * phi * f)
            i = map(complex S^1, source e, id_(source e), Degree => 1)
            assert(h*e == psi*i)
            assert((degree h, degree e, degree psi, degree i) === (0, 0, -1, 1))
        Text
            This routine is functorial.
        Example
            D' = (freeResolution coker matrix{{a^2,a*b,c^3}})[-1]
            f' = randomComplexMap(D', D)
            Hom(f' * f, g) == Hom(f, id_F) * Hom(f', g)
            Hom(f' * f, g) == Hom(f, g) * Hom(f', id_E)
            F' = dual (freeResolution coker matrix{{a^2,a*b,a*c,b^3}})[-3]
            g' = randomComplexMap(F', F)
            Hom(f, g' * g) == Hom(f, g') * Hom(id_D, g)
            Hom(f, g' * g) == Hom(id_C, g') * Hom(f, g)
    SeeAlso
        (homomorphism, ComplexMap)
        (homomorphism', ComplexMap)
        (randomComplexMap, Complex, Complex)
        (Hom, Complex, Complex)
///

doc ///
    Key
        (dual, ZZdFactorizationMap)
	(dual, ZZdFactorizationMap, RingElement)
	(dual, ZZdFactorizationMap, Symbol)
    Headline
        the dual of a map of ZZ/d-graded factorizations
    Usage
        h = dual f
	h = dual(f,t)
    Inputs
        f:ZZdFactorizationMap
	t:RingElement or Symbol
	     optional input specifying the symbol used to represent the distinguished root of unity, 
	     required only if the underlying ring does not have a specified root of unity and the period is > 2 
    Outputs
        h:ZZdFactorizationMap
    Description
        Text
            The map $f : C \to D$ of ZZ/d-graded factorizations over the ring $S$ induces the map
            $h = Hom(f, S^1) : Hom(D, S^1) \to Hom(C,S^1)$ defined by $\phi \mapsto \phi f$.
        Example
            S = ZZ/101[a..c]
            C = freeResolution coker vars S
            D = (freeResolution coker matrix{{a^2,a*b,b^3}})[-1]
            f = randomFactorizationMap(D,C)
            h = dual f
            assert isWellDefined h
            assert(h == Hom(f, S^1))
            assert(source h == Hom(D,S^1))
            assert(target h == Hom(C,S^1))
        Text
            This routine is functorial.
        Example
            D' = (freeResolution coker matrix{{a^2,a*b,c^3}})[-1]
            f' = randomComplexMap(D', D)
            dual(f' * f) == dual f * dual f'
    SeeAlso
        (Hom, ZZdFactorizationMap, ZZdFactorizationMap)
        (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
        (dual, Matrix)
///

///
    Key
        (End,ZZdFactorization)
	(End,ZZdFactorization,RingElement)
	(End,ZZdFactorization,Symbol)
    Headline
        Create the endomorphism complex of a ZZ/d-graded factorization
    Usage
        End(C)
	End(C,t)
    Inputs
        C:ZZdFactorization
	t:RingElement or Symbol
	     optional input specifying the symbol used to represent the distinguished root of unity, 
	     required only if the underlying ring does not have a specified root of unity and the period is > 2 
    Outputs
        :ZZdFactorization
    Description
        Text
        Example
    Caveat
    SeeAlso
///

///
    Key
        (End,ZZdFactorizationMap)
	(End,ZZdFactorizationMap,RingElement)
	(End,ZZdFactorizationMap,Symbol)
    Headline
        Create the endomorphism map of a ZZ/d-graded factorization map
    Usage
        End(phi)
	End(phi,t)
    Inputs
        phi:ZZdFactorizationMap
	t:RingElement or Symbol
	     optional input specifying the symbol used to represent the distinguished root of unity, 
	     required only if the underlying ring does not have a specified root of unity and the period is > 2 
    Outputs
        :ZZdFactorizationMap
    Description
        Text
        Example
    Caveat
    SeeAlso
///



doc ///
    Key
        (symbol**, ZZdFactorization, Matrix)
        (symbol**, Matrix, ZZdFactorization)
    Headline
        create the tensor product of a ZZ/d-graded factorization and a map of modules
    Usage
        h = C ** f
        h = f ** C
    Inputs
        C:ZZdFactorization
            over a ring $R$
        f:Matrix
            defining a homomorphism from the $R$-module $M$ to the $R$-module $N$
    Outputs
        h:ZZdFactorizationMap
            from $C \otimes M$ to $C \otimes N$
    Description
        Text
            For any ZZ/d-graded factorization map $C$, a map $f \colon M \to N$ of $R$-modules induces a
            morphism $C \otimes f$ of ZZ/d-graded factorizations
            from $C \otimes M$ to $C \otimes N$.  This method returns this map of factorizations.
        Example
            R = ZZ/101[a..d];
            I = ideal(c^2-b*d, b*c-a*d, b^2-a*c)
            J = ideal(I_0, I_1)
            C = koszulComplex vars R
            f = map(R^1/I, R^1/J, 1)
            C ** f
            f ** C
            f' = random(R^2, R^{-1, -1, -1})
            C ** f'
            f' ** C
            assert isWellDefined(C ** f')
            assert isWellDefined(f' ** C)
        Text
            Tensoring with a complex defines a functor from the category
            of $R$-modules to the category of complexes over $R$.
        Example
            f'' = random(source f', R^{-2,-2})
            assert((C ** f') * (C ** f'') == C ** (f' * f''))
            assert(C ** id_(R^{-1,-2,-3}) == id_(C ** R^{-1,-2,-3}))
    SeeAlso
        "Making maps between chain complexes"
        (symbol**, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (symbol**, ZZdFactorizationMap, ZZdFactorizationMap)
        (tensor, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol**, ZZdFactorization, ZZdFactorizationMap)
        (symbol**, ZZdFactorizationMap, ZZdFactorization)
        (symbol**, ZZdFactorizationMap, Module)
        (symbol**, Module, ZZdFactorizationMap)
    Headline
        the map of ZZ/d-graded factorizations between tensor factorizations
    Usage
        h = f ** g
        h = tensor(f, g)
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
    Outputs
        h:ZZdFactorizationMap
    Description
        Text
            The maps $f : C \to D$ and $g : E \to F$ of ZZ/d-graded factorizations induce the map
            $h = f \otimes g : C \otimes E \to D \otimes F$ defined by $c \otimes e \mapsto f(c) \otimes g(e)$.
        Example
            S = ZZ/101[a..c]
            C = freeResolution coker vars S
            D = (freeResolution coker matrix{{a^2,a*b,b^3}})[-1]
            f = randomComplexMap(D,C)
            E = (dual C)[-3]
            F = (dual D)[-3]
            g = randomComplexMap(F,E)
            h = f ** g
            assert isWellDefined h
            assert(source h === C ** E)
            assert(target h === D ** F)
        Text
            If one argument is a Complex or Module,
            then the identity map of the corresponding complex is used.
        Example
            fE = f ** E
            assert(fE === f ** id_E)
            k = coker vars S
            gk = g ** k
            assert(gk == g ** id_(complex k))
        Text
            This routine is functorial.
        Example
            D' = (freeResolution coker matrix{{a^2,a*b,c^3}})[-1]
            f' = randomComplexMap(D', D)
            (f' * f) ** g == (f' ** g) * (f ** id_E)
            (f' * f) ** g == (f' ** id_F) * (f ** g)
            F' = dual (freeResolution coker matrix{{a^2,a*b,a*c,b^3}})[-3]
            g' = randomComplexMap(F', F)
            f ** (g' * g) == (f ** g') * (id_C ** g)
            f ** (g' * g) == (id_D ** g') * (f ** g)
    SeeAlso
        (symbol**, ZZdFactorization, ZZdFactorization)
        (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
        (Hom, ZZdFactorization, ZZdFactorization)
///




doc ///
    Key
        (symbol SPACE, RingMap, ZZdFactorizationMap)
    Headline
        apply a ring map to a map of ZZ/d-graded factorizations
    Usage
        phi f
    Inputs
        phi:RingMap
            whose source is a ring $R$, and whose target is a ring $S$
        f:ZZdFactorizationMap
            over the ring $R$
    Outputs
        :ZZdFactorizationMap
            over the ring $S$
    Description
        Text
            We illustrate the image of a ZZ/d-graded factorization map along a ring map.
        Example
            R = QQ[a,b,c,d];
            S = QQ[s,t];
            phi = map(S, R, {s, s+t, t, s-t})
            I = ideal(a*b, b*c, c*d)
            J = I + ideal(a^2, b^2, c^2, d^2)
            CI = freeResolution I
            CJ = freeResolution J
            f = extend(CJ, CI, map(CJ_0, CI_0, 1))
            assert isWellDefined f
            g = phi f
            assert isWellDefined g
            dd^(source g)
            dd^(target g)
            prune HH g
    SeeAlso
        (symbol SPACE, RingMap, ZZdFactorization)
        (symbol **, RingMap, ZZdFactorizationMap)
///

doc ///
    Key
        (symbol**, RingMap, ZZdFactorizationMap)
        (symbol**, Ring, ZZdFactorizationMap)
        (symbol**, ZZdFactorizationMap, RingMap)
        (symbol**, ZZdFactorizationMap, Ring)
        (tensor, RingMap, ZZdFactorizationMap)
        (tensor, ZZdFactorizationMap, RingMap)
    Headline
        tensor a map of ZZ/d-graded factorizations along a ring map
    Usage
        phi ** f
        tensor(phi, f)
        S ** f
        f ** S
    Inputs
        phi:RingMap
            whose source is a ring $R$, and whose target is a ring $S$
        f:ZZdFactorizationMap
            over the ring $R$
    Outputs
        :ZZdFactorizationMap
            over the ring $S$
    Description
        Text
            These methods implement the base change of rings.  As input, one can either
            give a ring map $\phi$, or the ring $S$ (when there is a canonical map
                from $R$ to $S$).
        Text
            We illustrate the tensor product of a map of ZZ/d-graded factorizations along a ring map.
        Example
            R = QQ[a,b,c,d];
            S = QQ[s,t];
            phi = map(S, R, {s, s+t, t, s-t})
            I = ideal(a*b, b*c, c*d)
            J = I + ideal(a^2, b^2, c^2, d^2)
            CI = freeResolution I
            CJ = freeResolution J
            f = extend(CJ, CI, map(CJ_0, CI_0, 1))
            assert isWellDefined f
            g = phi ** f
            assert isWellDefined g
            dd^(source g)
            dd^(target g)
            prune HH g
    SeeAlso
        (symbol **, RingMap, ZZdFactorization)
        (symbol SPACE, RingMap, ZZdFactorizationMap)
///
 




doc ///
    Key
        canonicalMap --ISTHIS IMPLEMENTED???
        (canonicalMap, ZZdFactorization, ZZdFactorization)
        [canonicalMap, UseTarget]
        UseTarget
    Headline
        gets the natural map arising from various constructions
    Usage
        g = canonicalMap(D, C)
    Inputs
        C:Complex
        D:Complex
        UseTarget => Boolean
            determines the choice of canonical map
            when $D$ is a cylinder of a map $f$
            and the source and target of $f$ are the same
    Outputs
        g:ZZdFactorizationMap
    Description
        Text
            A canonical map, also called a natural map, is a 
            map that arises naturally from the definition or
            the construction of the object.
            
            The following six constructions are supported: kernel, 
            cokernel, image, coimage, cone, and cylinder.
        Text
            The @TO2((kernel, ComplexMap), "kernel of a complex map")@ 
            comes with a natural injection into the source complex.
            This natural map is always a complex morphism.
        Example
            R = ZZ/101[a,b,c,d];
            D = freeResolution coker vars R
            C = (freeResolution coker matrix"a,b,c")[1]
            f = randomComplexMap(D, C, Cycle=>true)
            assert isComplexMorphism f
            K1 = kernel f
            g = canonicalMap(source f, K1)
            degree g
            assert(isWellDefined g and isComplexMorphism g)
        Example
            f2 = randomComplexMap(D, C)
            assert not isComplexMorphism f2
            K2 = kernel f2
            g2 = canonicalMap(source f2, K2)
            assert(isWellDefined g2 and isComplexMorphism g2)
        Text
            The @TO2((cokernel, ComplexMap), "cokernel of a complex map")@
            comes with a natural surjection from the target complex.
        Example
            Q = cokernel f
            g3 = canonicalMap(Q, target f)
            assert(isWellDefined g3 and isComplexMorphism g3)
        Text
            The @TO2((image, ComplexMap), "image of a complex map")@ 
            comes with a natural injection into the target complex.
        Example
            I = image f
            g4 = canonicalMap(target f, I)
            assert(isWellDefined g4 and isComplexMorphism g4)
        Text
            The @TO2((coimage, ComplexMap), "coimage of a complex map")@
            comes with a natural surjection from the source complex.
            This natural map is always a complex morphism.
        Example
            J = coimage f
            g5 = canonicalMap(J, source f)
            assert(isWellDefined g5 and isComplexMorphism g5)
        Example
            J2 = coimage f2
            g6 = canonicalMap(J2, source f2)
            assert(isWellDefined g6 and isComplexMorphism g6)
        Text
            The @TO2((cone, ComplexMap), "cone of a complex morphism")@
            comes with two natural maps.  Given a map $f : C \to D$,
            let $E$ denote the cone of $f$.  The first is a natural
            injection from the target $D$ of $f$ into $E$.  The
            second is a natural surjection from $E$ to $C[-1]$.
            Together, these maps form a short exact sequence of
            complexes.
        Example
            E = cone f
            g = canonicalMap(E, target f)
            h = canonicalMap((source f)[-1], E)
            assert(isWellDefined g and isWellDefined h)
            assert(isComplexMorphism g and isComplexMorphism h)
            assert isShortExactSequence(h,g)
        Text
            The @TO2((cylinder, ComplexMap), "cylinder of a complex
            map")@ comes with four natural maps.  Given a map $f : C
            \to D$, let $F$ denote the cylinder of $f$.  The first is
            the natural injection from the source $C$ of $f$ into the
            cylinder $F$.  Together these two maps form a short exact
            sequence of complexes.
        Example
            F = cylinder f
            g = canonicalMap(F, source f)
            h = canonicalMap(E, F)
            assert(isWellDefined g and isWellDefined h)
            assert(isComplexMorphism g and isComplexMorphism h)
            assert isShortExactSequence(h,g)
        Text
            The third is the natural injection from the target $D$ of $F$
            into the cylinder $F$.
            The fourth is the natural surjection from the cylinder $F$ to the 
            target $D$ of $f$.
            However, these two maps do not form a short exact
            sequence of complexes.
        Example
            g' = canonicalMap(F, target f)
            h' = canonicalMap(target f, F)
            assert(isWellDefined g' and isWellDefined h')
            assert(isComplexMorphism g' and isComplexMorphism h')
            assert not isShortExactSequence(h',g')
        Text
            When $D == C$, the optional argument {\tt UseTarget} 
            selects the appropriate natural map.
        Example
            f' = id_C
            F' = cylinder f'
            g = canonicalMap(F', C, UseTarget=>true)
            h = canonicalMap(F', C, UseTarget=>false)
            assert(isWellDefined g and isWellDefined h)
            assert(g != h)
            assert(isComplexMorphism g and isComplexMorphism h)
    SeeAlso
        (inducedMap, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (inducedMap, ZZdFactorization, ZZdFactorization)
    Headline
        make the map of ZZ/d-graded factorizations induced at each term by the identity map
    Usage
        f = inducedMap(D, C)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        Degree => ZZ
            specify the degree of the map of ZZ/d-graded factorizations, if not 0
        Verify => Boolean
            if true, check that the resulting maps are well-defined
    Outputs
        f:ZZdFactorizationMap
    Description
        Text
            Let $d$ be the value of the optional argument {\tt
            Degree}, or zero, if not given.  For each $i$, the terms
            $D_{i+d}$ and $C_i$ must be subquotients of the same
            ambient free module.  This method returns the complex map
            induced by the identity on each of these free modules.
            
            If {\tt Verify => true} is given, then this method
            also checks that these identity maps induced well-defined 
            maps.  This can be a relatively expensive computation.
        Text
            We illustrate this method by truncating a free resolution
            at two distinct internal degrees.  We check that 
            the various induced maps compose to give another
            induced map.
        Example
            needsPackage "Truncations"
            kk = ZZ/32003
            R = kk[a,b,c]
            F = freeResolution (ideal gens R)^2
            C1 = truncate(3, F)
            C2 = truncate(4, F)
            assert isWellDefined C1
            assert isWellDefined C2
            f = inducedMap(C1, C2)
            assert isWellDefined f
            f1 = inducedMap(F, C1)
            f2 = inducedMap(F, C2)
            assert isWellDefined f1
            assert isWellDefined f2
            assert(f2 == f1 * f)
    SeeAlso
        (inducedMap, Module, Module)
        "Truncations :: truncate(ZZ,Complex)"
///

doc ///
    Key
        "arithmetic with ZZ/d-graded factorization maps"
        (symbol+, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol+, RingElement, ZZdFactorizationMap)
        (symbol+, Number, ZZdFactorizationMap)
        (symbol+, ZZdFactorizationMap, RingElement)
        (symbol+, ZZdFactorizationMap, Number)
        (symbol-, ZZdFactorizationMap)
        (symbol-, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol-, RingElement, ZZdFactorizationMap)
        (symbol-, Number, ZZdFactorizationMap)
        (symbol-, ZZdFactorizationMap, RingElement)
        (symbol-, ZZdFactorizationMap, Number)
        (symbol*, RingElement, ZZdFactorizationMap)
        (symbol*, Number, ZZdFactorizationMap)
        (symbol*, ZZdFactorizationMap, RingElement)
        (symbol*, ZZdFactorizationMap, Number)
    Headline
        perform arithmetic operations on ZZ/d-graded factorization maps
    Usage
        f + g
        a + f
        f + a
        -f
        f - g
        a - f
        f - a
        a * f
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
        a:RingElement
          that is, an element in the underlying ring or a number
    Outputs
        :ZZdFactorizationMap
    Description
        Text
            The set of ZZ/d-graded factorizations maps forms a module over the underlying @TO2((ring, ZZdFactorizationMap), "ring")@.
            These methods implement the basic operations of addition, subtraction, and scalar multiplication.
        Example
            R = ZZ/101[a..d];
            C = freeResolution coker matrix{{a*b, a*c^2, b*c*d^3, a^3}}
            D = freeResolution coker matrix{{a*b, a*c^2, b*c*d^3, a^3, a*c*d}}
            f = randomFactorizationMap(D, C, Cycle => true)
            g = randomFactorizationMap(D, C, Boundary => true)
        Example
            f+g
            f-g
            -f
            3*f
            0*f
            a*f
            assert(0*f == 0)
            assert(1*f == f)
            assert((-1)*f == -f)
            assert(-(f-g) == g-f)
            assert((a+b)*f == a*f + b*f)
            assert(a*(f+g) == a*f + a*g)
            assert isFactorizationMorphism (f+g)
        Text
            Adding or subtracting a scalar is the same as adding or subtracting the
            scalar multiple of the identity.  In particular, the source and target must be equal.
        Example
            h = randomFactorizationMap(C, C)
            h+1
            assert(h+1 == h + id_C)
            assert(h+a == h + a*id_C)
            assert(1-h == id_C - h)
            assert(b-c*h == -c*h + b*id_C)
            assert(b-h*c == -h*c + id_C*b)
        Text
            Arithmetic on differentials can be a useful method
            for constructing new chain complexes.
        Example
            E = complex(-dd^C)
            isWellDefined E
            assert(dd^E == map(E, E, -dd^C))
    SeeAlso
        "Making maps between chain complexes"
        randomFactorizationMap
        (ZZdfactorization, ZZdFactorizationMap)
        (map, ZZdFactorization, ZZdFactorization, ZZdFactorizationMap)
///

doc ///
    Key
        (symbol|, ZZdFactorizationMap, ZZdFactorizationMap)
    Headline
        join or concatenate maps horizontally
    Usage
        f | g
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
    Outputs
        :ZZdFactorizationMap
    Description
        Text
            Given ZZ/d-graded factorization maps with the same target,
            this method constructs the associated map
            from the direct sum of the sources to the target.

            First, we define some non-trivial maps of chain complexes.
        Example
            R = ZZ/101[a..d];
            C1 = (freeResolution coker matrix{{a,b,c}})[1]
            C2 = freeResolution coker matrix{{a*b,a*c,b*c}}
            D = freeResolution coker matrix{{a^2,b^2,c*d}}
            f = randomComplexMap(D, C1)
            g = randomComplexMap(D, C2)
        Example
            h = f|g
            assert isWellDefined h
            assert(source h === source f ++ source g)
            assert(target h === target f)
        Text
            This is really a shorthand for constructing complex maps via block matrices.
        Example
            assert(h === map(D, C1 ++ C2, {{f,g}}))
    SeeAlso
        (symbol++, ZZdFactorization, ZZdFactorization)
        (symbol++, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol||, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol|, Matrix, Matrix)
        (map, ZZdFactorization, ZZdFactorization, List)
///

doc ///
    Key
        (symbol||, ZZdFactorizationMap, ZZdFactorizationMap)
    Headline
        join or concatenate maps vertically
    Usage
        f || g
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
    Outputs
        :ZZdFactorizationMap
    Description
        Text
            Given ZZ/d-graded factorization maps with the same source,
            this method constructs the associated map
            from the source to the direct sum of the targets.

            First, we define some non-trivial maps of chain complexes.
        Example
            R = ZZ/101[a..d];
            D1 = (freeResolution coker matrix{{a,b,c}})[1]
            D2 = freeResolution coker matrix{{a*b,a*c,b*c}}
            C = freeResolution coker matrix{{a^2,b^2,c*d}}
            f = randomComplexMap(D1, C)
            g = randomComplexMap(D2, C)
        Example
            h = f||g
            assert isWellDefined h
            assert(target h === target f ++ target g)
            assert(source h === source f)
        Text
            This is really a shorthand for constructing complex maps via block matrices.
        Example
            assert(h === map(D1 ++ D2, C, {{f},{g}}))
    SeeAlso
        (symbol++, ZZdFactorization, ZZdFactorization)
        (symbol++, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol|, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol||, Matrix, Matrix)
        (map, ZZdFactorization, ZZdFactorization, List)
///

doc ///
    Key
        (symbol++, ZZdFactorizationMap, ZZdFactorizationMap)
        (directSum, ZZdFactorizationMap)
    Headline
        direct sum of ZZ/d-graded factorization maps
    Usage
        h = f ++ g
        h = directSum(f,g,...)
        h = directSum(name1 => f, name2 => g, ...)
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
    Outputs
        h:ZZdFactorizationMap
          that is the direct sum of the input ZZ/d-graded factorization maps
    Description
        Text
            The direct sum of two ZZ/d-graded factorization maps is a a ZZ/d-graded factorization map
            from the direct sum of the sources to the direct sum of
            the targets.

            First, we define some non-trivial maps of ZZ/d-graded factorizations.
        Example
            R = ZZ/101[a..d];
            C1 = (freeResolution coker matrix{{a,b,c}})[1]
            C2 = freeResolution coker matrix{{a*b,a*c,b*c}}
            D1 = (freeResolution coker matrix{{a,b,c}})
            D2 = freeResolution coker matrix{{a^2, b^2, c^2}}[-1]
            f = randomComplexMap(D1, C1, Cycle => true)
            g = randomComplexMap(D2, C2, Cycle => true)
        Example
            h = f ++ g
            assert isWellDefined h
            assert(h == map(D1 ++ D2, C1 ++ C2, {{f,0},{0,g}}))
        Text
            The direct sum of any sequence of complex maps can be 
            computed as follows.
        Example
            directSum(f, g, f[2])
            h2 = directSum(mike => f, greg => g, dan => f[2])
            h2_[greg,dan]
            assert(source oo == C2 ++ C1[2])
        Text
            One can easily obtain the compositions with canonical
            injections and surjections.
        Example
            h_[0]^[0] == f
            h_[1]^[1] == g
            h_[0]^[1] == 0
            h_[1]^[0] == 0
        Example
            h_[0] == h * (C1 ++ C2)_[0]
            h_[1] == h * (C1 ++ C2)_[1]
            h^[0] == (D1 ++ D2)^[0] * h
            h^[1] == (D1 ++ D2)^[1] * h
    SeeAlso
        (symbol++, ZZdFactorization, ZZdFactorization)
        (symbol**, ZZdFactorizationMap, ZZdFactorizationMap)
        (Hom, ZZdFactorizationMap, ZZdFactorizationMap)
        (symbol_, ZZdFactorizationMap, Array)
///

-- TODO for the next 4 nodes:
-- make better examples
-- add text (when are these actually defined?  Maybe change code)
-- also add SeeAlso canonicalMap's.
doc ///
  Key
    (image, ZZdFactorizationMap)
  Headline
    make the image of a map of ZZ/d-graded factorizations
  Usage
    E = image f
  Inputs
    f : ZZdFactorizationMap
  Outputs
    E : ZZdFactorization
  Description
    Text
      If $f : C \to D$ is a map of ZZ/d-graded factorizations of degree $d$,
      then the image is the ZZ/d-graded factorization $E$ whose $i-th$ is $image(f_{i-d})$,
      and whose differential is induced from the differential 
      on the target.
    Text
      In the following example, we first construct a random
      complex morphism $f : C \to D$.  We consider 
      the exact sequence $0 \to D \to cone(f) \to C[-1] \to 0$.
      For the maps $g : D \to cone(f)$ and $h : cone(f) \to C[-1]$,
      we compute the image.
    Example
      S = ZZ/101[a,b,c,d];
      C = freeResolution ideal(b^2-a*c, b*c-a*d, c^2-b*d)
      D = freeResolution ideal(a,b,c)
      f = randomComplexMap(D, C, Cycle => true, InternalDegree => 0)
      Cf = cone f
      g = canonicalMap(Cf, D)
      h = canonicalMap(C[-1], Cf)
      prune image g == D
      prune image h == C[-1]
    Text
      There is a canonical map of complexes from the image to the target.
    Example
      g1 = canonicalMap(target g, image g)
      ker g1 == 0
      image g1 == image g
      h1 = canonicalMap(target h, image h)
      ker h1 == 0
      image h1 == image h
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    image
    (coimage, ComplexMap)
    (kernel, ComplexMap)
    (cokernel, ComplexMap)
    canonicalMap
///

doc ///
  Key
    (coimage, ZZdFactorizationMap)
  Headline
    make the coimage of a map of ZZ/d-graded factorizations
  Usage
    coimage f
  Inputs
    f : ZZdFactorizationMap
  Outputs
    : ZZdFactorization
  Description
    Text
      The coimage of a ZZ/d-graded factorization map $f : C \to D$
      is the ZZ/d-graded factorization $E$ whose $i-th$ term is $coimage(f_i)$,
      and whose differential is induced from the differential 
      on the source.
    Text
      In the following example, we first construct a random
      complex morphism $f : C \to D$.  We consider 
      the exact sequence $0 \to D \to cone(f) \to C[-1] \to 0$.
      For the maps $g : D \to cone(f)$ and $h : cone(f) \to C[-1]$,
      we compute the coimage.
    Example
      S = ZZ/101[a,b,c,d];
      C = freeResolution ideal(b^2-a*c, b*c-a*d, c^2-b*d)
      D = freeResolution ideal(a,b,c)
      f = randomComplexMap(D, C, Cycle => true, InternalDegree => 0)
      Cf = cone f
      g = canonicalMap(Cf, D)
      h = canonicalMap(C[-1], Cf)
      coimage g == D
      prune coimage h == C[-1]
    Text
      There is a canonical map of complexes from the source to the coimage.
    Example
      g1 = canonicalMap(coimage g, source g)
      coimage g1 == coimage g
      coker g1 == 0
      h1 = canonicalMap(coimage h, source h)
      coimage h1 == coimage h
      coker h1 == 0
  Caveat
    The coimage is more computationally intensive than @TO (image, ComplexMap)@
    because, unlike {\tt image}, it computes kernels of maps of modules.
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    coimage
    (image, ZZdFactorizationMap)
    (kernel, ZZdFactorizationMap)
    (cokernel, ZZdFactorizationMap)
    canonicalMap
///

doc ///
  Key
    (kernel, ZZdFactorizationMap)
  Headline
    make the kernel of a map of ZZ/d-graded factorizations
  Usage
    kernel f
    ker f
  Inputs
    f : ZZdFactorizationMap
  Outputs
    : ZZdFactorization
  Description
    Text
      The kernel of a chain complex map $f : C \to D$
      is the complex $E$ whose $i-th$ term is $kernel(f_i)$,
      and whose differential is induced from the differential 
      on the source.
    Text
      In the following example, we first construct a random
      complex morphism $f : C \to D$.  We consider 
      the exact sequence $0 \to D \to cone(f) \to C[-1] \to 0$.
      For the maps $g : D \to cone(f)$ and $h : cone(f) \to C[-1]$,
      we compute the kernel.
    Example
      S = ZZ/101[a,b,c,d];
      C = freeResolution ideal(b^2-a*c, b*c-a*d, c^2-b*d)
      D = freeResolution ideal(a,b,c)
      f = randomComplexMap(D, C, Cycle => true, InternalDegree => 0)
      Cf = cone f
      g = canonicalMap(Cf, D)
      h = canonicalMap(C[-1], Cf)
      ker g == 0
      prune ker h == D
    Text
      There is a canonical map of complexes from the kernel to the source.
    Example
      h1 = canonicalMap(source h, ker h)
      ker h == image h1
      ker h1 == 0
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    ker
    (image, ZZdFactorizationMap)
    (coimage, ZZdFactorizationMap)
    (cokernel, ZZdFactorizationMap)
    canonicalMap
///

doc ///
  Key
    (cokernel, ZZdFactorizationMap)
  Headline
    make the cokernel of a map of ZZ/d-graded factorizations
  Usage
    cokernel f
    coker f
  Inputs
    f : ZZdFactorizationMap
  Outputs
    : ZZdFactorization
  Description
    Text
      If $f : C \to D$ is a map of chain complexes of degree $d$,
      then the cokernel is the complex $E$ whose $i-th$ is $cokernel(f_{i-d})$,
      and whose differential is induced from the differential 
      on the target.
    Text
      In the following example, we first construct a random
      complex morphism $f : C \to D$.  We consider 
      the exact sequence $0 \to D \to cone(f) \to C[-1] \to 0$.
      For the maps $g : D \to cone(f)$ and $h : cone(f) \to C[-1]$,
      we compute the kernel.
    Example
      S = ZZ/101[a,b,c,d];
      C = freeResolution ideal(b^2-a*c, b*c-a*d, c^2-b*d)
      D = freeResolution ideal(a,b,c)
      f = randomComplexMap(D, C, Cycle => true, InternalDegree => 0)
      Cf = cone f
      g = canonicalMap(Cf, D)
      h = canonicalMap(C[-1], Cf)
      prune coker g == C[-1]
      coker h == 0
    Text
      There is a canonical map of complexes from the target to the cokernel.
    Example
      g1 = canonicalMap(coker g, target g)
      coker g == image g1
      coker g1 == 0
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    cokernel
    (image, ZZdFactorizationMap)
    (coimage, ZZdFactorizationMap)
    (kernel, ZZdFactorizationMap)
    canonicalMap
///


doc ///
  Key
    (cone, ZZdFactorizationMap)
  Headline
    make the mapping cone of a morphism of ZZ/d-graded factorizations
  Usage
    cone f
  Inputs
    f:ZZdFactorizationMap
      which is a morphism of ZZ/d-graded factorizations
  Outputs
    :ZZdFactorization
  Description
    Text
      Given a morphism $f \colon B \to C$, the mapping cone is the complex
      whose $i$-th term is $B_{i-1} \oplus C_i$, and whose $i$-th 
      differential is given by
      \[ \begin{bmatrix} -\operatorname{dd}^{B[-1]} & 0 \\ f[-1] & \operatorname{dd}^C \end{bmatrix}. \]
    Text
      A map between modules induces a map between their free resolutions,
      and we compute the associated mapping cone.
    Example
      S = ZZ/32003[x,y,z];
      M = ideal vars S
      B = freeResolution(S^1/M^2)
      C = freeResolution(S^1/M)
      f = extend(C,B,id_(S^1))
      Cf = cone f
      dd^Cf
      prune HH Cf
      assert(prune HH_1 Cf == prune(M/M^2))
    Text
      The mapping cone fits into a canonical short exact
      sequence of chain complexes:
      $$0 \to C \to \operatorname{cone}(f) \to B[-1] \to 0.$$
    Example
      g = canonicalMap(Cf,C)
      h = canonicalMap(B[-1],Cf)
      assert(isWellDefined g and isWellDefined h)
      assert(isShortExactSequence(h,g))
    Text
      The most important application of mapping cones is to 
      identify quasi-isomorphisms: $f$ is a quasi-isomorphism 
      if and only if the mapping cone is acyclic.
    Example
      aug = augmentationMap C
      assert isWellDefined aug
      cone aug
      assert(0 == prune HH cone aug)
      assert isQuasiIsomorphism aug
    Text
      Mapping cones can also be used to construct free resolutions
      of subschemes linked via a complete intersection to a
      arithmetically Cohen-Macaulay subscheme;
      see Peskine-Szpiro, Liaison des variétés algébriques I, 
          {\it Invent. math.} {\bf 26} (1974) 271-302.
    Text
      Here, we consider a random complete intersection of 2 cubics
      contained in the ideal of the twisted cubic curve, and we
      compute a free resolution of the linked curve of degree 6.
    Example
      S = ZZ/32003[a..d]
      I = monomialCurveIdeal(S, {1,2,3})
      K = ideal((gens I) * random(source gens I, S^{-3,-3}))
      C = freeResolution(S^1/I)
      B = freeResolution(S^1/K)
      f = dual extend(C,B,id_(S^1))
      Cf = (cone f)[-2]
      prune HH Cf
      Cf' = minimize Cf
      J = ideal dd^Cf'_1
      freeResolution J
      assert(degree J == 6)
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    (cylinder, ZZdFactorizationMap)
    (extend, ZZdFactorization, ZZdFactorization, Matrix)
    canonicalMap
    isQuasiIsomorphism
    isShortExactSequence
///


doc ///
  Key
    cylinder
    (cylinder, ZZdFactorizationMap)
  Headline
    make the mapping cylinder of a morphism of ZZ/d-graded factorizations
  Usage
    cylinder f
  Inputs
    f:ZZdFactorizationMap
      which is a morphism of complexes
  Outputs
    :ZZdFactorization
  Description
    Text
      Given a morphism $f : B \to C$, the mapping cylinder 
      is the complex whose the $i$-th term is $B_{i-1} \oplus B_i \oplus C_i$
      and whose $i$-th differential is given in block form by
              {\tt matrix \{\{ - dd^B_{i-1}, 0, 0 \}, 
                \{ -id_{B_{i-1}}, dd^B_i, 0 \},
                \{ f_{i-1}, 0, dd^C_i\}\}}.
      Alternatively, the cylinder is the
      mapping cone of the morphism $g : B \to B \oplus C$ given in block form
      by
        {\tt matrix\{\{-id_B\}, \{f\}\}}.
    Text
      A map between modules induces a map between their free resolutions,
      and we compute the associated mapping cylinder.
    Example
      S = ZZ/32003[x,y,z];
      M = ideal vars S
      B = freeResolution(S^1/M^2)
      C = freeResolution(S^1/M)
      f = extend(C,B,id_(S^1))
      cylf = cylinder f
      dd^cylf
      assert isWellDefined cylf
    Text
      The mapping cylinder fits into a canonical short exact
      sequence of chain complexes,
      $$0 \to B \to cyl(f) \to cone(f) \to 0.$$
    Example
      Cf = cone f
      g = canonicalMap(cylf, B)
      h = canonicalMap(Cf, cylf)
      assert(isWellDefined g and isWellDefined h)
      assert(isShortExactSequence(h,g))
    Text
      The alternative interpretation of the cylinder, defined above,
      can be demonstrated as follows.
    Example
      g = map(B ++ C, B, {{-id_B},{f}})
      cone g == cylf
  SeeAlso
    "Making chain complexes"
    "Making maps between chain complexes"
    (cone, ComplexMap)
    (extend, Complex, Complex, Matrix)
    (freeResolution, Module)
    canonicalMap
    isShortExactSequence
///

doc ///
    Key
        (symbol^, ZZdFactorizationMap, Array)
        (symbol_, ZZdFactorizationMap, Array)
    Headline
        the composition with the canonical inclusion or projection map
    Usage
        i = f_[name]
        p = f^[name]
    Inputs
        f:ZZdFactorizationMap
        name:
    Outputs
        :ZZdFactorizationMap
            {\tt i} is the composition of {\tt f} with the canonical inclusion and {\tt p} is
            the composition of the canonical projection with {\tt f}
    Description
        Text
            The direct sum is an n-ary operator with projection and
            inclusion maps from each component satisfying appropriate
            identities.

            One can access these maps as follows.  First, we define
            some non-trivial maps of chain complexes.
        Example
            R = ZZ/101[a..d];
            C1 = (freeResolution coker matrix{{a,b,c}})[1]
            C2 = freeResolution coker matrix{{a*b,a*c,b*c}}
            D1 = (freeResolution coker matrix{{a,b,c}})
            D2 = freeResolution coker matrix{{a^2, b^2, c^2}}[-1]
            f = randomComplexMap(D1, C1, Cycle => true)
            g = randomComplexMap(D2, C2, Cycle => true)
        Example
            h = f ++ g
        Text
            The four basic maps are the inclusion from each summand of the source
            and the projection to each summand of the target.
        Example
            h_[0] == h * (C1 ++ C2)_[0]
            h_[1] == h * (C1 ++ C2)_[1]
            h^[0] == (D1 ++ D2)^[0] * h
            h^[1] == (D1 ++ D2)^[1] * h
        Text
            These can be combined to obtain the blocks of the map of chain complexes.
        Example
            h_[0]^[0] == f
            h_[1]^[1] == g
            h_[0]^[1] == 0
            h_[1]^[0] == 0
            assert(h == map(D1 ++ D2, C1 ++ C2, {{f,0},{0,g}}))
        Text
            The default names for the components are the non-negative
            integers.  However, one can choose any name.
        Example
            h = (mike => f) ++ (greg => g)
            h_[mike]^[mike] == f
            h_[greg]^[greg] == g
    SeeAlso
        (symbol++, ZZdFactorization, ZZdFactorization)
        (symbol^, ZZdFactorization, Array)
        (symbol_, ZZdFactorization, Array)
        (directSum, ZZdFactorization)
        (components, ZZdFactorization)
        indices
///

doc ///
    Key
        (randomFactorizationMap, ZZdFactorization, ZZdFactorization)
        randomFactorizationMap
        [randomFactorizationMap, Boundary]
        [randomFactorizationMap, Cycle]
        [randomFactorizationMap, Degree]
        [randomFactorizationMap, InternalDegree]
        Cycle
        Boundary
        InternalDegree
    Headline
        a random map of ZZ/d-graded factorizations
    Usage
        f = randomFactorizationMap(C,D)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        Boundary => Boolean
            whether the constructed {\tt f} is a null homotopy
        Cycle => Boolean
            whether the constructed {\tt f} commutes with the differentials
        Degree => ZZ
            the degree of the constructed map of chain complexes
        InternalDegree => List
            or @ofClass ZZ@
    Outputs
        f:ZZdFactorizationMap
    Description
        Text
            A random ZZ/d-graded factorization map $f : C \to D$ is obtained from a random element
            in the ZZ/d-graded factorization @TO2 ((Hom,ZZdFactorization,ZZdFactorization), "$Hom(C,D)$")@.
        Example
            S = ZZ/101[a..c]
            C = freeResolution coker matrix{{a*b, a*c, b*c}}
            D = freeResolution coker vars S
            f = randomComplexMap(D,C)
            assert isWellDefined f
            assert not isCommutative f
            assert not isNullHomotopic f
        Text
            When the random element in the complex $Hom(C,D)$ lies in the kernel
            of the differential, the associated map of complexes commutes
            with the differential.
        Example
            g = randomComplexMap(D,C, Cycle => true)
            assert isWellDefined g
            assert isCommutative g
            assert isComplexMorphism g
            assert not isNullHomotopic g
        Text
            When the random element in the complex $Hom(C,D)$ lies in the image
            of the differential, the associated map of complexes is a null
            homotopy.
        Example
            h = randomComplexMap(D,C, Boundary => true)
            assert isWellDefined h
            assert isCommutative h
            assert isComplexMorphism h
            assert isNullHomotopic h
            nullHomotopy h
        Text
            When the degree of the random element in the complex $Hom(C,D)$ is non-zero,
            the associated map of complexes has the same degree.
        Example
            p = randomComplexMap(D, C, Cycle => true, Degree => -1)
            assert isWellDefined p
            assert isCommutative p
            assert not isComplexMorphism p
            assert(degree p === -1)
        Text
            By default, the random element is constructed as a random linear combination of
            the basis elements in the appropriate degree of $Hom(C,D)$.  Given an internal
            degree, the random element is constructed as maps of modules with this degree.
        Example
            q = randomComplexMap(D, C, Boundary => true, InternalDegree => 2)
            assert all({0,1,2}, i -> degree q_i === {2})
            assert isHomogeneous q
            assert isWellDefined q
            assert isCommutative q
            assert isComplexMorphism q
            source q === C
            target q === D
            assert isNullHomotopic q
    SeeAlso
        (homomorphism, ZZdFactorizationMap)
        (homomorphism', ZZdFactorizationMap)
        (Hom, ZZdFactorization, ZZdFactorization)
///

doc ///
    Key
        (homology, ZZdFactorizationMap)
        (homology, ZZ, ZZdFactorizationMap)
        (cohomology, ZZ, ComplexMap)
    Headline
        induced map on homology or cohomology --DO WE WANT COHOMOLOGY?
    Usage
        h = HH f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        h:ZZdFactorizationMap
    Description
        Text
            Homology defines a functor from the category of chain complexes
            to itself.  Given a map of chain complexes $f : C \to D$,
            this method returns the induced map $HH f : HH C \to HH D$.
        Text
            To directly obtain the $n$-th map in $h$, use {\tt HH_n f} or
            {\tt HH^n f}.  By definition {\tt HH^n f === HH_(-n) f}.
            This can be more efficient, as it will compute only the desired
            induced map.
        Text
            If $f$ commutes with the differentials, then these induced
            maps are well defined.
        Example
            S = ZZ/101[a..d]
            I = ideal(a*b, a*d, c*b, c*d)
            C = (dual freeResolution I)[1]
            D = dual complex for i from 0 to 4 list koszul(i,gens I)
            assert isWellDefined D
            f = randomComplexMap(D, C, Cycle => true)
            assert isCommutative f
            h = HH f
            assert isWellDefined h
            prune h
            assert(source h == HH C)
            assert(target h == HH D)
        Example
            f2 = randomComplexMap(D, C, Cycle => true, Degree => -1)
            h2 = HH f2
            assert isWellDefined h2
            prune h2
        Text
            A boundary will always induce the zero map.
        Example
            f3 = randomComplexMap(D, C, Boundary => true)
            h3 = HH f3
            assert isWellDefined h3
            assert(h3 == 0)
    SeeAlso
        (homology, ZZdFactorization)
        (homology, ZZ, ZZdFactorization)
        (cohomology, ZZ, ZZdFactorization)
        (prune, ZZdFactorizationMap)
///


doc ///
    Key
        (tensorCommutativity, ZZdFactorization, ZZdFactorization)
    Headline
        make the canonical isomorphism arising from commutativity of the tensor product operation
    Usage
        tensorCommutativity(C, D)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
            both over the same ring $R$
    Outputs
        :ZZdFactorizationMap
            that is an isomorphism from $C \otimes_R D$ to 
            $D \otimes_R C$
    Description
        Text
            The commutativity of tensor products of modules induces
            the commutativity of tensor products of ZZ/d-graded factorizations. The main difference
	    for factorizations is that there are coefficients depending on the degrees of the terms
	    being commuted.
            This method implements this isomorphism for ZZ/d-graded factorizations.
        Text
            Using two term complexes of small rank,
            we see that this isomorphism need not be the identity map.
        Example
            S = ZZ/101[x_0..x_8];
            C = complex{genericMatrix(S,x_0,2,1)}
            D = complex{genericMatrix(S,x_2,1,2)}
            F = C ** D
            G = D ** C
            f = tensorCommutativity(C,D)
            assert isWellDefined f
            assert isComplexMorphism f
            assert(source f === F)
            assert(target f === G)
            assert(f_1 != id_(source f_1))
            assert(prune ker f == 0)
            assert(prune coker f == 0)
            g = f^-1
            assert isWellDefined g
            assert(g * f == 1)
            assert(f * g == 1)
        Text
            We illustrate this isomorphism on complexes, none
            of whose terms are free modules.
        Example
            ses = (I,J) -> (
                complex{
                    map(S^1/(I+J), S^1/I ++ S^1/J, {{1,1}}),
                    map(S^1/I ++ S^1/J, S^1/(intersect(I,J)), {{1},{-1}})
                    }
                )
            C = ses(ideal(x_0,x_1), ideal(x_1,x_2))
            D = ses(ideal(x_3,x_4,x_5), ideal(x_6,x_7,x_8))
            h = tensorCommutativity(C, D);
            assert isWellDefined h
            assert isComplexMorphism h
            assert(ker h == 0)
            assert(coker h == 0)
            k = h^-1;
            assert(h*k == 1)
            assert(k*h == 1)
            h_2
            assert(source h_2 != target h_2)
        Text
            Interchanging the arguments gives the inverse map.
        Example
            h1 = tensorCommutativity(D, C)
            assert isComplexMorphism h1
            assert(h1*h == id_(C**D))
            assert(h*h1 == id_(D**C))
        Text
            Interchanging the factors in a tensor product of
            two complex maps can be accomplished as follows.
        Example
            C = freeResolution ideal(x_0^2, x_1^2)
            D = freeResolution ideal(x_0, x_1)
            f = extend(D, C, map(D_0, C_0, 1))
            E = freeResolution ideal(x_2^3, x_3^3, x_4^3)
            F = freeResolution ideal(x_2, x_3, x_4)
            g = extend(F, E, map(F_0, E_0, 1))
            assert(tensorCommutativity(D,F) * (f**g) == (g**f) * tensorCommutativity(C,E))
            assert isComplexMorphism tensorCommutativity(D,F)
            assert isComplexMorphism tensorCommutativity(C,E)
    SeeAlso
        "Working with Tor"
        (tensorCommutativity, Module, Module)
        (tensorAssociativity, Complex, Complex, Complex)
        (isComplexMorphism, ComplexMap)
///

doc ///
    Key
        (tensorAssociativity, ZZdFactorization, ZZdFactorization, ZZdFactorization)
    Headline
        make the canonical isomorphism arising from associativity
    Usage
        tensorAssociativity(C, D, E)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        E:ZZdFactorization
    Outputs
        :ZZdFactorizationMap
            which is an isomorphism from {\tt C ** (D ** E)} to 
            {\tt (C ** D) ** E} 
    Description
        Text
            The associativity of tensor products of modules induces
            the associativity of tensor products of ZZ/d-graded factorizations.
            This method implements this isomorphism for ZZ/d-graded factorizations.
        Text
            Using two term complexes of small rank,
            we see that this isomorphism need not be the identity map.
        Example
            S = ZZ/101[x_0..x_11]
            C = complex{genericMatrix(S,x_0,2,1)}
            D = complex{genericMatrix(S,x_4,1,2)}
            E = complex{genericMatrix(S,x_8,2,2)}
            F = (C ** D) ** E
            G = C ** (D ** E)
            f = tensorAssociativity(C,D,E)
            assert isWellDefined f
            assert(source f === G)
            assert(target f === F)
            f_1
            assert(f_1 != id_(source f_1))
            assert(prune ker f == 0)
            assert(prune coker f == 0)
            g = f^-1
            assert isWellDefined g
            assert(g * f == 1)
            assert(f * g == 1)
        Text
            We illustrate this isomorphism on complexes, none
            of whose terms are free modules.
        Example
            ses = (I,J) -> (
                complex{
                    map(S^1/(I+J), S^1/I ++ S^1/J, {{1,1}}),
                    map(S^1/I ++ S^1/J, S^1/(intersect(I,J)), {{1},{-1}})
                    }
                )
            C = ses(ideal(x_0,x_1), ideal(x_1,x_2))
            D = ses(ideal(x_3,x_4,x_5), ideal(x_6,x_7,x_8))
            E = ses(ideal(x_1^2, x_1*x_2), ideal(x_1*x_3,x_9^2))
            h = tensorAssociativity(C, D, E);
            assert isWellDefined h
            assert(ker h == 0)
            assert(coker h == 0)
            k = h^-1;
            assert(h*k == 1)
            assert(k*h == 1)
            h_2
            assert(source h_2 != target h_2)
    SeeAlso
        "Working with Tor"
        (tensorCommutativity, Complex, Complex)
        (tensorAssociativity, Module, Module, Module)
///

doc ///
    Key
        (isShortExactSequence, ZZdFactorizationMap, ZZdFactorizationMap)
    Headline
        whether a pair of ZZ/d-graded factorization maps forms a short exact sequence
    Usage
        isShortExactSequence(g, f)
    Inputs
        f:ZZdFactorizationMap
        g:ZZdFactorizationMap
    Outputs
        :Boolean
            that is @TO true@ if these form a short exact sequence
    Description
        Text
            A short exact sequence of complexes 
            \[ 0 \to B \xrightarrow{f} C \xrightarrow{g} D \to 0\]
            consists of two morphisms of complexes
            $f \colon B \to C$ and $g \colon C \to D$ such that
            $g f = 0$, $\operatorname{image} f = \operatorname{ker} g$, 
            $\operatorname{ker} f = 0$, and $\operatorname{coker} g = 0$.
        Text
            From a complex morphism $h \colon B \to C$, one obtains a
            short exact sequence
            \[ 0 \to \operatorname{image} h \to C \to \operatorname{coker} h \to 0. \]
        Example
            R = ZZ/101[a,b,c];
            B = freeResolution coker matrix{{a^2*b, a*b*c, c^3}}
            C = freeResolution coker vars R
            h = randomComplexMap(C, B, Cycle => true)
            f = canonicalMap(C, image h)
            g = canonicalMap(coker h, C)
            assert isShortExactSequence(g,f)
        Text
            A short exact sequence of modules gives rise to a short
            exact sequence of complexes.  These complexes arise
            as free resolutions of the modules.
        Example
            I = ideal(a^3, b^3, c^3)
            J = I + ideal(a*b*c)
            K = I : ideal(a*b*c)
            SES = complex{
                map(comodule J, comodule I, 1),
                map(comodule I, (comodule K) ** R^{-3}, {{a*b*c}})
                }
            assert isWellDefined SES
            assert isShortExactSequence(dd^SES_1, dd^SES_2)
            (g,f) = horseshoeResolution SES
            assert isShortExactSequence(g,f)
    SeeAlso
        "Basic invariants and properties"
        canonicalMap
        (cone, ComplexMap)
        (longExactSequence, ComplexMap, ComplexMap) --YOU SHOULD IMPLEMENT THE TRIANGLE
///




doc ///
    Key
        (isQuasiIsomorphism, ZZdFactorizationMap)
	[isQuasiIsomorphism, Concentration]
	isQuasiIsomorphism
    Headline
         whether a map of ZZ/d-graded factorizations is a quasi-isomorphism
    Usage
         isQuasiIsomorphism f
    Inputs
         f:ZZdFactorizationMap
         Concentration => Sequence
             restricts attention to the induced maps indexed
             by elements in the given interval 
    Outputs
         :Boolean
             that is true when $f$ is a morphism of ZZ/d-graded factorizations
             such that the induced maps on homology are all
             isomorphisms
    Description
        Text
            The @TO2((cone, ZZdFactorizationMap), "cone")@ of a 
            map $f \colon C \to D$ is acyclic
            exactly when $f$ is a quasi-isomorphism.
        Example
            S = ZZ/32003[x,y,z];
            C = freeResolution coker vars S
            f = augmentationMap C
            assert isQuasiIsomorphism f
            assert(0 == prune HH cone f)
            assert isIsomorphism HH_0 f
            assert isIsomorphism HH_1 f
        Text
            XXX TODO. Free resolutions of complexes produce quasi 
            isomorphisms. (use example to doc of (resolution, Complex)).
        Example
            D = complex{random(S^2, S^{-3,-3,-4})}
            prune HH D
    SeeAlso
        "Basic invariants and properties"
        (cone, ZZdFactorizationMap)
        liftMapAlongQuasiIsomorphism
///

doc ///
    Key
        (isNullHomotopyOf, ZZdFactorizationMap, ZZdFactorizationMap)
        isNullHomotopyOf
    Headline
        whether the first map of chain complexes is a null homotopy for the second
    Usage
        isNullHomotopyOf(h, f)
    Inputs
        h:ZZdFactorizationMap
        f:ZZdFactorizationMap
    Outputs
        :Boolean
            that is true when $h$ is a null homotopy of $f$
    Description
        Text
            A map of ZZ/d-graded factorizations $f \colon C \to D$ is
            null-homotopic if there exists a map of chain
            complexes $h : C \to D$ of degree $\deg(f)+1$,
            such that we have the equality 
            \[ f = \operatorname{dd}^D h 
              + (-1)^{\deg(f)} h \operatorname{dd}^C.
            \]
        Text
            As a first example, we construct a map of chain complexes
            in which the null homotopy is given by the identity.
        Example
            R = ZZ/101[x,y,z];
            M = cokernel matrix{{x,y,z^2}, {y^2,z,x^2}}
            C = complex {id_M}
            h = map(C, C, i -> if i == 0 then id_M, Degree => 1)
            isWellDefined h
            assert isNullHomotopyOf(h, id_C)
            assert isNullHomotopic id_C
        Text
            A random map of chain complexes, arising as a boundary
            in the associated Hom complex, is automatically
            null homotopic.  We use the method @TO nullHomotopy@
            to construct a witness and verify it is a null homotopy.
        Example
            C = (freeResolution M) ** R^1/ideal(x^3, z^3-x)
            f = randomComplexMap(C, C[1], Boundary => true)
            assert isNullHomotopic f
            h = nullHomotopy f
            assert isNullHomotopyOf(h, f)
        Text
            By assigning @TO "debugLevel"@ a positive value,
            this method provides some information about the nature
            of the failure to be a null homotopy.
        Example
            g1 = randomComplexMap(C, C[1], Degree => 1)
            g2 = randomComplexMap(C, C[1], Degree => -1)
            debugLevel = 1
            assert not isNullHomotopyOf(g1, f)
            assert not isNullHomotopyOf(g2, f)
    SeeAlso
        "Basic invariants and properties"
        (isNullHomotopic, ComplexMap)
        (nullHomotopy, ComplexMap)
        randomComplexMap
        (Hom, Complex, Complex)
///

doc ///
    Key
        (isNullHomotopic, ZZdFactorizationMap)
        isNullHomotopic
    Headline
        whether a map of ZZ/d-graded factorizations is null-homotopic
    Usage
        isNullHomotopic f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        :Boolean
            that is true when $f$ is null-homotopic
    Description
        Text
            A map of ZZ/d-graded factorizations $f \colon C \to D$ is
            null-homotopic if there exists a map of chain
            complexes $h : C \to D$ of degree $\deg(f)+1$,
            such that we have the equality 
            \[ f = \operatorname{dd}^D h 
              + (-1)^{\deg(f)} h \operatorname{dd}^C.
            \]
        Text
            As a first example, we construct a map of chain complexes
            in which the null homotopy is given by the identity.
        Example
            R = ZZ/101[x,y,z];
            M = cokernel matrix{{x,y,z^2}, {y^2,z,x^2}}
            C = complex {id_M}
            assert isNullHomotopic id_C
            h = nullHomotopy id_C
            assert(h_0 == id_M)
            assert isNullHomotopyOf(h, id_C)
        Text
            A random map of chain complexes, arising as a boundary
            in the associated Hom complex, is automatically
            null homotopic.
        Example
            C = (freeResolution M) ** R^1/ideal(x^3, z^3-x)
            f = randomComplexMap(C, C[1], Boundary => true)
            assert isNullHomotopic f
            h = nullHomotopy f
            assert isNullHomotopyOf(h, f)
            g = randomComplexMap(C, C[1])
            assert not isNullHomotopic g
        Text
            This procedure also works for complex maps
            whose degree is non-zero.
        Example
            f = randomComplexMap(C, C[2], Boundary => true, Degree => 1)
            assert isNullHomotopic f
            h = nullHomotopy f
            assert isNullHomotopyOf(h, f)
    SeeAlso
        "Basic invariants and properties"
        (isNullHomotopyOf, ComplexMap, ComplexMap)
        (nullHomotopy, ComplexMap)
        randomComplexMap
        (Hom, Complex, Complex)
///

doc ///
    Key
        (nullHomotopy, ZZdFactorizationMap)
        nullHomotopy
    Headline
        a map which is a candidate for being a null homotopy
    Usage
        h = nullHomotopy f
    Inputs
        f:ZZdFactorizationMap
    Outputs
        h:ZZdFactorizationMap
    Description
        Text
            A map of chain complexes $f \colon C \to D$ is
            null-homotopic if there exists a map of chain
            complexes $h : C \to D$ of degree $\deg(f)+1$,
            such that we have the equality 
            \[ f = \operatorname{dd}^D h 
              + (-1)^{\deg(f)} h \operatorname{dd}^C.
            \]
            Given $f$, this method returns a map $h$ of chain complexes
            that will be a null-homotopy if one exists.
        Text
            As a first example, we construct a map of chain complexes
            in which the null homotopy is given by the identity.
        Example
            R = ZZ/101[x,y,z];
            M = cokernel matrix{{x,y,z^2}, {y^2,z,x^2}}
            C = complex {id_M}
            assert isNullHomotopic id_C
            h = nullHomotopy id_C
            assert(h_0 == id_M)
            assert isNullHomotopyOf(h, id_C)
        Text
            A random map of chain complexes, arising as a boundary
            in the associated Hom complex, is automatically
            null homotopic.
        Example
            C = (freeResolution M) ** R^1/ideal(x^3, z^3-x)
            f = randomComplexMap(C, C[1], Boundary => true)
            assert isNullHomotopic f
            h = nullHomotopy f
            assert isNullHomotopyOf(h, f)
        Text
            When a map of chain complexes is not null-homotopic,
            this method nevertheless returns a map $h$ of
            chain complexes, having the correct source, target
            and degree, but cannot be a null homotopy.
        Example
            g = randomComplexMap(C, C[1])
            assert not isNullHomotopic g
            h' = nullHomotopy g
            assert isWellDefined h'
            assert(degree h' === degree g + 1)
            assert not isNullHomotopyOf(h', g)
        Text
            For developers: when the source of $f$ is a free complex,
            a procedure, that is often faster, is attempted.  In the
            general case this method uses the Hom complex.
    Caveat
        The output is only a null homotopy when one exists.
    SeeAlso
        "Making maps between chain complexes"
        (isNullHomotopic, ZZdFactorizationMap)
        (isNullHomotopyOf, ZZdFactorizationMap, ZZdFactorizationMap)
        randomFactorizationMap
        (Hom, ZZdFactorization, ZZdFactorization)
///



doc ///
    Key
        (extend, ZZdFactorization, ZZdFactorization, Matrix, Sequence)
        (extend, ZZdFactorization, ZZdFactorization, Matrix)
    Headline
        extend a map of modules to a map of ZZ/d-graded factorizations, if possible
    Usage
        g = extend(D, C, f, p)
        g = extend(D, C, f)
    Inputs
        C:ZZdFactorization
        D:ZZdFactorization
        f:Matrix
        p:Sequence
            consisting of a pair of integers $(j,i)$, such that the
            matrix $f$ defines a map from $C_i$ to $D_j$; the default 
            value is $i = j = 0$
        Verify => Boolean
            currently, this option is ignored
    Outputs
        :ZZdFactorizationMap
    Description
        Text
            Let $C$ be a chain complex such that each term is a free
            module.  Let $D$ be a chain
            complex which is exact at the $k$-th term for all $k > j$.
            Given a map of modules $f \colon C_i \to D_j$ such that
            the image of $f \circ \operatorname{dd}^C_{i+1}$ is
            contained in the image of $\operatorname{dd}^D_{j+1}$,
            this method constructs a morphism of chain complexes $g
            \colon C \to D$ of degree $j-i$ such that $g_i = f$.

            $\phantom{WWWW}
            \begin{array}{cccccc}
            0 & \!\!\leftarrow\!\! & C_{i} & \!\!\leftarrow\!\! & C_{i+1} & \!\!\leftarrow\!\! & C_{i+2} & \dotsb \\
              &            & \downarrow \, {\scriptstyle f} & & \downarrow \, {\scriptstyle g_{i+1}} && \downarrow \, {\scriptstyle g_{i+2}} \\
            0 & \!\!\leftarrow\!\! & D_{j} & \!\!\leftarrow\!\! & D_{j+1} &  \!\!\leftarrow\!\! & D_{j+2} & \dotsb \\
            \end{array}
            $
        Text
            A map between modules extends to a map between their free resolutions.
        Example
            S = ZZ/101[a..d];
            I = ideal(a*b*c, b*c*d, a*d^2)
            C = S^{{-3}} ** freeResolution (I:a*c*d)
            D = freeResolution I
            f = map(D_0, C_0, matrix{{a*c*d}})
            g = extend(D, C, f)
            assert isWellDefined g
            assert isComplexMorphism g
            assert(g_0 == f)
            E = cone g
            dd^E
        Text
            Extension of maps to complexes is also useful in 
            constructing a free resolution of a linked ideal.
        Example
            I = monomialCurveIdeal(S, {1,2,3})
            K = ideal(I_1^2, I_2^2)
            FI = freeResolution I
            FK = freeResolution K
            f = map(FI_0, FK_0, 1)
            g = extend(FI, FK, f)
            assert isWellDefined g
            assert isComplexMorphism g
            assert(g_0 == f)
            C = cone (dual g)[- codim K]
            dd^C
            dd^(minimize C)
            assert(ideal relations HH_0 C == K:I)
        Text
            Inspired by a @TO yonedaMap@ computation, we extend a map
            of modules to a map between free resolutions having
            homological degree $-1$.
        Example
            f = map(FK_0, FI_1, matrix {{a*c^2-a*b*d, -b*c^2+a*c*d, -c^3+a*d^2}}, Degree => 1)
            assert isHomogeneous f
            assert isWellDefined f
            g = extend(FK, FI, f, (0,1))
            assert isWellDefined g
            assert isCommutative g
            assert(degree g === -1)
            assert isHomogeneous g
    SeeAlso
        "Making maps between chain complexes"
        (cone, ZZdFactorizationMap)
        (isFactorizationMorphism, ComplexMap)
        (minimize, ZZdFactorization)
///

doc ///
    Key
        (liftMapAlongQuasiIsomorphism, ZZdFactorizationMap, ZZdFactorizationMap)
        liftMapAlongQuasiIsomorphism
        (symbol//, ZZdFactorizationMap, ZZdFactorizationMap)
        (quotient, ZZdFactorizationMap, ZZdFactorizationMap)
        homotopyMap
        (homotopyMap, ZZdFactorizationMap)
    Headline
        lift a map of ZZ/d-graded factorizations along a quasi-isomorphism
    Usage
        f' = liftMapAlongQuasiIsomorphism(f, g)
        f' = f // g
    Inputs
        f:ZZdFactorizationMap
            where each term in the source of $f$ is a free module
        g:ZZdFactorizationMap
            a quasi-isomorphism having the same target as $f$
    Outputs
        f':ZZdFactorizationMap
            a map from the source of $f$ to the source of $g$
    Consequences
        Item
            the homotopy relating $f$ and $g \circ f'$ is
            available as {\tt homotopyMap f'}.
    Description
        Text
            Let $f \colon P \to C$ be a morphism of chain 
            complexes, where each term in $P$ is a free module.
            Given a quasi-isomorphism $g \colon B \to C$,
            this method produces a morphism $f' \colon P \to B$
            such that there exists a map $h \colon P \to C$
            of chain complexes having degree $1$
            satisfying

            $f - g \circ f' = h \circ \operatorname{dd}^P +
             \operatorname{dd}^C \circ h$.
             
        Text
            Given a morphism between complexes, we can construct 
            the corresponding map 
            between their
            free resolutions using this method.
            
            To be more precise, 
            given a morphism $\phi \colon B \to C$ of complexes,
            let $\alpha \colon P \to B$ and 
            $\beta \colon F \to C$ denote the free resolutions
            of the source and target complexes.
            Lifting the composite map $\phi \circ \alpha$ along the
            quasi-isomorphism $\beta$ gives a commutative diagram
            $\phantom{WWWW}
            \begin{array}{ccc}
            P & \!\!\rightarrow\!\! & F \\
            \downarrow \, {\scriptstyle \alpha} & & \downarrow \, {\scriptstyle \beta} \\
            B & \xrightarrow{\phi} & C
            \end{array}
            $
        Example
            S = ZZ/101[a,b,c,d];
            J = ideal(a*b, a*d, b*c);
            I = J + ideal(c^3);
            C = prune Hom(S^{2} ** freeResolution I, S^1/I)
            D = prune Hom(freeResolution J, S^1/J)
            r = randomComplexMap(D,C,Cycle=>true)
            f = r * resolutionMap C
            g = resolutionMap D
            assert isQuasiIsomorphism g
            f' = liftMapAlongQuasiIsomorphism(f, g)
            assert(f' == f//g)
            assert isWellDefined f'
            assert isComplexMorphism f'
            h = homotopyMap f'
            isNullHomotopyOf(h, g * (f//g) - f)
        Text
            TODO: XXX start here. Do triangles, invert interesting quasi-isomorphism.
            isSemiFree, and add in an example or 2.  Include
            finding an inverse for a quasi-isomorphism.
            We need some kind of better example here.
    Caveat
        The following three assumptions are not checked:
        $f$ is a morphism, the source of $f$ is semifree,
        and $g$ is a quasi-isomorphism.
    SeeAlso
        "Towards computing in the derived category"
        isQuasiIsomorphism
        isComplexMorphism
///



doc ///
    Key
        (connectingMap, ZZdFactorizationMap, ZZdFactorizationMap)
        connectingMap
        [connectingMap, Concentration]
    Headline
        construct the connecting homomorphism on homology
    Usage
        connectingMap(g, f)
    Inputs
        f:ZZdFactorizationMap
            an injective morphism $f \colon A \to B$
        g:ZZdFactorizationMap
            a surjective morphism $g \colon B \to C$ 
            whose kernel is the same as the image of $f$
        Concentration => Sequence
            not yet implemented
    Outputs
        h:ZZdFactorizationMap
            a ZZ/d-graded factorization morphism
            whose source is the homology of $C$ and whose target
            is the homology of $A$, shifted by $-1$
    Description
        Text
            Given a short exact sequence of ZZ/d-graded factorizations

            $\phantom{WWWW}
            0 \leftarrow C \xleftarrow{g} B \xleftarrow{f} A \leftarrow 0,
            $
            
            this function returns the unique morphism $h \colon H(C) \to H(A)[-1]$ of ZZ/d-graded factorization
            that naturally fits into the ZZ/d-graded sequence

            $\phantom{WWWW}
            \dotsb \leftarrow H(C)[-1] \xleftarrow{H(g)[-1]} H(B)[-1] \xleftarrow{H(f)[-1]} H(A)[-1] \xleftarrow{h} H(C) \xleftarrow{H(g)} H(B) \xleftarrow{H(f)} H(A) \leftarrow \dotsb.
            $
            
            $\phantom{WWWW}$
            
        Text
            As a first example, consider a free resolution $F$ of $S/I$.
            Applying the Hom functor $\operatorname{Hom}(F, -)$ to a short exact sequence of modules

            $\phantom{WWWW}
            0 \leftarrow S/h \leftarrow S \xleftarrow{h} S(- \deg h) \leftarrow 0
            $

            gives rise to a short exact sequence of complexes.  The corresponding long exact sequence in homology
            has the form

            $\phantom{WWWW}
            \dotsb \leftarrow \operatorname{Ext}^{d+1}(S/I, S(-\deg h))
            \xleftarrow{\delta} 
            \operatorname{Ext}^d(S/I, S/h)
            \leftarrow \operatorname{Ext}^d(S/I, S)
            \leftarrow \operatorname{Ext}^d(S/I, S(-\deg h))
            \leftarrow \dotsb.
            $
        Example
            S = ZZ/101[a..d, Degrees=>{2:{1,0},2:{0,1}}];
            h = a*c^2 + a*c*d + b*d^2;
            I = (ideal(a,b) * ideal(c,d))^[2]
            F = freeResolution comodule I;
            g = Hom(F, map(S^1/h, S^1, 1))
            f = Hom(F, map(S^1, S^{-degree h}, {{h}}))
            assert isWellDefined g
            assert isWellDefined f
            assert isShortExactSequence(g, f)
            delta = connectingMap(g, f)
            assert isWellDefined delta
            assert(degree delta == 0)            
            assert(source delta_(-1) == Ext^1(comodule I, S^1/h))
            assert(target delta_(-1) == Ext^2(comodule I, S^{{-1,-2}}))
            L = longExactSequence(g,f)
            assert isWellDefined L
            assert(HH L == 0)
            assert(dd^L_-9 === delta_-3)
            assert(dd^L_-8 === HH_-3 g)
            assert(dd^L_-7 === HH_-3 f)
            assert(dd^L_-6 === delta_-2)
            assert(dd^L_-5 === HH_-2 g)
            assert(dd^L_-4 === HH_-2 f)
            assert(dd^L_-3 === delta_-1)
        Text
            Applying the Hom functor $\operatorname{Hom}(-, S)$ to the horseshoe resolution of
            a short exact sequence of modules

            $\phantom{WWWW}
            0 \leftarrow S/(I+J) \leftarrow S/I \oplus S/J  \leftarrow S/I \cap J \leftarrow 0
            $

            gives rise to a short exact sequence of complexes.  The corresponding long exact sequence in homology
            has the form

            $\phantom{WWWW}
            \dotsb \leftarrow \operatorname{Ext}^{d+1}(S/(I+J), S)
            \xleftarrow{\delta} 
            \operatorname{Ext}^d(S/I \cap J, S)
            \leftarrow \operatorname{Ext}^d(S/I \oplus S/J, S)
            \leftarrow \operatorname{Ext}^d(S/(I+J), S)
            \leftarrow \dotsb.
            $
        Example
            S = ZZ/101[a..d];
            I = ideal(c^3-b*d^2, b*c-a*d)
            J = ideal(a*c^2-b^2*d, b^3-a^2*c)
            ses = complex{
                map(S^1/(I+J), S^1/I ++ S^1/J, {{1,1}}),
                map(S^1/I ++ S^1/J, S^1/intersect(I,J), {{1},{-1}})
                }
            assert isWellDefined ses
            assert(HH ses == 0)
            (g,f) = horseshoeResolution ses
            assert isShortExactSequence(g,f)
            (Hf, Hg) = (Hom(f, S), Hom(g, S));
            assert isShortExactSequence(Hf, Hg)
            delta = connectingMap(Hf, Hg)
            assert isWellDefined delta
            assert isComplexMorphism delta
            assert(source delta_-2 == Ext^2(comodule intersect(I,J), S))
            assert(target delta_-2 == Ext^3(comodule (I+J), S))
            L = longExactSequence(Hf, Hg)
            assert isWellDefined L
            assert(HH L == 0)
            assert(dd^L_-6 === delta_-3)
            assert(dd^L_-5 === HH_-3 Hf)
            assert(dd^L_-4 === HH_-3 Hg)
            assert(dd^L_-3 === delta_-2)
            assert(dd^L_-2 === HH_-2 Hf)
            assert(dd^L_-1 === HH_-2 Hg)
            assert(dd^L_0 === delta_-1)
    SeeAlso
        "Towards computing in the derived category"
        (longExactSequence, ZZdFactorizationMap, ZZdFactorizationMap)
///


///
    Key
       (Fold, Complex, ZZ)
       (Fold, ComplexMap, ZZ)
    Headline
        Convert any complex or complex map into a ZZ/d-graded factorization (or map) for a fixed integer d
    Usage
        Fold(C,d)
	Fold(phi,d)
    Inputs
        C:Complex
	phi:ComplexMap
	d:ZZ
	   an integer specifying the period of the resulting factorization or factorization map
    Outputs
        :ZZdFactorization
    Description
        Text
        Example
    Caveat
    SeeAlso
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



