
---------------------------------
---------------------------------
-- **DOCUMENTATION SECTION** --
---------------------------------
---------------------------------

doc ///
    Key
        MatrixSchubert
    Headline
        matrix Schubert varieties and ASM varieties
    Description
        Text
            This package provides functions for constructing and investigating matrix Schubert varieties. Many of the functions in this package can take as input either a permutation in 1-line notation, or an alternating sign matrix.
        Text
            @UL {
	    {"[CV20] Aldo Conca and Matteo Varbaro, ",
	    HREF("https://arxiv.org/abs/1805.11923", EM "Square-free Grobner degenerations"),
	    ", Inventiones mathematicae, 221(3), pp.713-730."},
	    {"[Ful92] William Fulton, ",
	    HREF("https://sites.math.washington.edu/~billey/classes/schubert.library/fulton.essential.set.pdf",
		EM "Flags, Schubert polynomials, degeneracy loci, and determinantal formulas"),
	    ", Duke Math J. 65 (1992): 381-420."},
            {"[KM05] Allen Knutson and Ezra Miller, ",
            HREF("https://arxiv.org/abs/math/0110058", EM "Grobner geometry of Schubert polynomials"),
            ", Annals of Mathematics (2005): 1245-1318."},
	    {"[KW21] Patricia Klein and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2108.08370", EM "Bumpless pipe dreams encode Grobner geometry of Schubert polynomials"),
	    ", arxiv preprint 2108.08370."},
            {"[PSW21] Oliver Pechenik, David Speyer, and Anna Weigandt, ",
            HREF("https://arxiv.org/abs/2111.10681", EM "Castelnuovo-Mumford regularity of matrix Schubert varieties"),
            ", arxiv preprint 2111.10681."},
            {"[Wei17] Anna Weigandt, ",
            HREF("https://arxiv.org/abs/1708.07236", EM "Prism tableaux for alternating sign matrix varieties"),
            ", arXiv preprint 1708.07236."}
            }@  
        Text
            @SUBSECTION "Contributors"@
        Text
            The following people have generously contributed code, improved existing code, or enhanced the documentation:
            @HREF("https://sites.google.com/illinois.edu/shiliang-gao", "Shiliang Gao")@,
            @HREF("https://www.math.tamu.edu/directory/graduate.html", "Pooja Joshi")@, and
            @HREF("https://www.clemson.edu/science/academics/departments/mathstat/about/profiles/arakoto", "Antsa Tantely Fandresena Rakotondrafara")@.
    SeeAlso 
        "Investigating matrix Schubert varieties"
        "Investigating ASM varieties"
        "Initial ideals of ASM ideals"
        "Functions for investigating permutations"
///
------------------------------
------------------------------
-- * See Also Doc Nodes * --
------------------------------
------------------------------

doc ///
    Key 
        "Investigating matrix Schubert varieties"
    Headline
        basic functions for Schubert determinantal ideals
    Description
    	Text
	    Matrix Schubert varieties were introduced by Fulton [Ful92] in the study of Schubert
	    varieties in the complete flag variety. Their defining ideals are called Schubert
	    determinantal ideals.
        Text
            @UL {
	    {"[Ful92] William Fulton, ",
	    HREF("https://sites.math.washington.edu/~billey/classes/schubert.library/fulton.essential.set.pdf",
		EM "Flags, Schubert polynomials, degeneracy loci, and determinantal formulas"),
	    " , Duke Math J. 65 (1992): 381-420."},
	    {"[HPW22] Zachary Hamaker, Oliver Pechenik, and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2003.13719", EM "Grobner geometry of Schubert polynomials through ice"),
	    " , Advances in Mathematics 398 (2022): 108228."},
            {"[KM05] Allen Knutson and Ezra Miller, ",
            HREF("https://arxiv.org/abs/math/0110058", EM "Grobner geometry of Schubert polynomials"),
            " , Annals of Mathematics (2005): 1245-1318."},
            {"[PSW21] Oliver Pechenik, David Speyer, and Anna Weigandt, ",
            HREF("https://arxiv.org/abs/2111.10681", EM "Castelnuovo-Mumford regularity of matrix Schubert varieties"),
            " , arxiv preprint 2111.10681."}
            }@
	Text
	    The general method for creating a
	    Schubert determinantal Ideal is @TO schubDetIdeal@.
	    The input is a permutation in the form of a list.
	    This package contains functions for investigating the rank matrix,
	    the Rothe diagram, and the essential cells of the Rothe diagram
	    as defined by Fulton in [Ful92].
	Example
	    p = {2,1,6,3,5,4};
	    rotheDiagram p
	    essentialSet p
	    rankTable p
	    netList fultonGens p	    
	Text
	    The default presentation given by @TO schubDetIdeal@ is given by the Fulton generators of the ideal.
	    In order to access a minimal generating set, use @TO trim@.
	Example
	    I = schubDetIdeal p;
	    # (I_*)
	    # ((trim I)_*)
	Text
	    After creating a Schubert determinantal ideal, the permutation associated to it is stored in its cache table.
	Example
	    peek I.cache
	Text
	    This package also contains methods for investigating antidiagonal initial ideals
	    of Schubert determinantal ideals and their associated Stanley-Reisner complexes,
	    which are a kind of subword complex. Subword complexes were introduced in [KM05] 
	    in the study of antidiagonal initial ideal of Schubert determinantal ideals
	    and their relation to Schubert polynomials.
	Example
	    antiDiagInit p
	    subwordComplex p
	Text
	    Given a list of permutations, this package also contains functions for intersecting and adding
	    the Schubert determinantal ideals associated to the list of permutations.
	Example
	    L = {{3,1,5,4,2},{2,5,3,4,1}} -- a list of 2 permutations
	    schubAdd L
	    schubIntersect L
	Text
	    Finally, this package contains functions for investigating homological invariants of matrix Schubert
	    varieties efficiently through combinatorial algorithms produced in [PSW21].
	Example
	    time schubReg p
	    time regularity comodule I 
	Text
	    @SUBSECTION "Functions for investigating matrix Schubert varieties"@
	Text
	    @UL {
		TO (antiDiagInit, List),
		TO (rankTable, List),
		TO (rotheDiagram, List),
		TO (augmentedRotheDiagram, List),
		TO (essentialSet, List),
		TO (augmentedEssentialSet, List),
		TO (schubDetIdeal, List),
		TO (fultonGens, List),
		TO (subwordComplex, List),
		TO (schubIntersect, List),
		TO (schubAdd, List),
		TO (schubReg, List),
		TO (schubCodim, List)
		}@
///

doc ///
    Key 
        "Investigating ASM varieties"
    Headline
        basic functions for alternating sign matrix ideals
    Description
    	Text
	    Alternating sign matrix (ASM) varieties were introduced by Weigandt [Wei17]. ASM varieties  generalize
	    matrix Schubert varieties. This package contains functions for investigating ASM varieties and 
	    their defining ideals.
        Text
            @UL {
	    {"[CV20] Aldo Conca and Matteo Varbaro, ",
	    HREF("https://arxiv.org/abs/1805.11923", EM "Square-free Grobner degenerations"),
	    " , Inventiones mathematicae, 221(3), pp.713-730."},
	    {"[Wei17] Anna Weigandt, ",
            HREF("https://arxiv.org/abs/1708.07236", EM "Prism tableaux for alternating sign matrix varieties"),
            " , arXiv preprint 1708.07236."}
            }@
	Text
	    The general method for defining the ideal of an ASM variety is @TO schubDetIdeal@.
	    The input can be an alternating sign matrix or a partial alternating sign matrix,
	    which always define the same ideal as some other ASM by results in [Wei17].
	    This package contains functions for checking if a matrix is a partial ASM,
	    extending a partial ASM to an ASM, and computing the rank matrix for an ASM.
	Example
	    A = matrix{{0,0,0},{0,1,0},{1,-1,0}} --Example 3.15 in [Wei17]
	    isPartialASM A
	    A' = partialASMToASM A
	Text
	    This package contains functions for investigating the rank matrix,
	    the Rothe diagram, and the essential cells of the Rothe diagram
	    of a partial ASM as defined by Fulton in [Ful92] and Weigandt in [Wei17].
	Example
	    rotheDiagram A
	    essentialSet A
	    rankTable A
	    netList fultonGens A	    
	Text
	    The default presentation given by @TO schubDetIdeal@ is given by the Fulton generators of the ideal.
	    In order to access a minimal generating set, use @TO trim@.
	Example
	    I = schubDetIdeal A;
	    # (I_*)
	    # ((trim I)_*)
	Text
	    After creating an ASM ideal, the corresponding ASM is stored in the cache table.
	Example
	    peek I.cache
	Text
	    This package also contains methods for investigating initial ideals of ASM ideals.
	Example
	    antiDiagInit A
	Text
	    Every ASM ideal can be written as the intersection of Schubert determinantal ideals.
	    The function @TO schubDecomposition@ outputs the list of permutations whose corresponding
	    Schubert determinantal ideals form the prime components of the ASM ideal associated to A.
	Example
	    schubDecomposition I
	Text
	    Given a list of partial ASMs, this package also contains functions for intersecting and adding
	    the ASM ideals associated to the list of partial ASMs.
	    The function @TO schubAdd@ automatically stores the ASM associated to this new ideal
	    in its cache table.
	Example
	    B = matrix{{0,0,1,0,0},{0,1,-1,1,0},{1,-1,1,0,0},{0,1,0,-1,1},{0,0,0,1,0}}
	    L = {A, B} -- a list of 2 partial ASMs
	    J = schubAdd L
	    peek J.cache
	    K = schubIntersect L
	Text
	    Given an ideal, we can obtain the ASM corresponding to that ideal using @TO getASM@.
	    Once we do this, that ASM is stored in the ideal's cache table for easy access.
	Example
    	    getASM K
	    peek K.cache
	Text
	    Finally, this package contains functions for investigating homological invariants of ASM ideals
	    efficiently by computing the associated invariants for their antidiagonal initial ideals,
	    which are known to be squarefree by [Wei17]. Therefore the extremal Betti numbers
	    (which encode regularity, depth, and projective dimension) of ASM ideals coincide 
	    with those of their antidiagonal initical ideals by [CV20].
	Example
	    time schubReg B
	    time regularity comodule schubDetIdeal B
	Text
	    @SUBSECTION "Functions for investigating ASM varieties"@
	Text
	    @UL {
		TO (isPartialASM, Matrix),
		TO (partialASMToASM, Matrix),
		TO (rankTable, Matrix),
		TO (rankTableFromMatrix, Matrix),
		TO (rankTableToASM, Matrix),
		TO (isASMIdeal, Ideal),
		TO (isASMUnion, List),
		TO (rotheDiagram, Matrix),
		TO (augmentedRotheDiagram, Matrix),
		TO (essentialSet, Matrix),
		TO (augmentedEssentialSet, Matrix),
		TO (schubDetIdeal, Matrix),
		TO (fultonGens, Matrix),
		TO (schubIntersect, List),
		TO (schubAdd, List),
		TO (schubReg, Matrix),
		TO (schubCodim, Matrix),
		TO (permOverASM, Matrix),
		TO (toOneLineNotation, Matrix)
		}@
///
doc ///
    Key 
        "Initial ideals of ASM ideals"
    Headline
        basic functions for investigating initial ideals of ASM varieties
    Description
    	Text
	    By work of Knutson and Miller [KM05] and Weigandt [Wei17],
	    the Fulton generators form a Grobner basis for an ASM ideal with respect to any antidiagonal term order.
	    However, the Grobner bases for ASM ideals with respect to other term orders, including diagonal ones, 
	    remain largely mysterious.
        Text
	    @UL {
	    {"[CV20] Aldo Conca and Matteo Varbaro, ",
	    HREF("https://arxiv.org/abs/1805.11923", EM "Square-free Grobner degenerations"),
	    " , Inventiones mathematicae, 221(3), pp.713-730."},
	    {"[HPW22] Zachary Hamaker, Oliver Pechenik, and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2003.13719", EM "Grobner geometry of Schubert polynomials through ice"),
	    " , Advances in Mathematics 398 (2022): 108228."},
	    {"[Kle23] Patricia Klein, ",
	    HREF("https://arxiv.org/abs/2008.01717", EM "Diagonal degenerations of matrix Schubert varieties"),
	    " , Algebraic Combinatorics 6 (2023) no. 4, 1073-1094."},	
            {"[KM05] Allen Knutson and Ezra Miller, ",
            HREF("https://arxiv.org/abs/math/0110058", EM "Grobner geometry of Schubert polynomials"),
            " , Annals of Mathematics (2005): 1245-1318."},
	    {"[Knu09] Allen Knutson, ",
            HREF("https://arxiv.org/abs/0911.4941", EM "Frobenius splitting, point-counting, and degeneration"),
            " , arxiv preprint 0911.4941."},
	    {"[KW21] Patricia Klein and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2108.08370", EM "Bumpless pipe dreams encode Grobner geometry of Schubert polynomials"),
	    " , arxiv preprint 2108.08370."},
            {"[Wei17] Anna Weigandt, ",
            HREF("https://arxiv.org/abs/1708.07236", EM "Prism tableaux for alternating sign matrix varieties"),
            " , arXiv preprint 1708.07236."}
            }@
	Text
	    Given a permutation or a partial ASM, one may compute its antidiagonal initial ideal.
	    By [KM05] and [Wei17] or [Knu09], the Fulton generators form a Grobner basis for any ASM ideal with respect to
	    any antidiagonal term order.
        Example
	    w = {2,4,5,1,3};
	    I = schubDetIdeal w;
	    inI = antiDiagInit w;
	    (netList sort inI_*, netList sort (trim I)_*)
	Text
	    In particular, by work of Conca and Varbaro [CV21], we know that the extremal Betti numbers (such as regularity
	    and projective dimension) of an ASM ideal and its antidiagonal initial ideal must coincide because the antidiagonal initial ideal is squarefree.
	    Curiously, for this running example, all of the Betti numbers coincide (not just the extremal ones).
	Example
	    (betti res I, betti res inI)
	Text
	    By work of Knutson and Miller [KM05], the prime components of the antidiagonal initial ideal
	    of a Schubert determinantal ideal for a permutation w are in bijection with the pipe dreams associated
	    to the permutation w.
	    If w is a permutation of n, then a pipe dream is a filling of an nxn grid with
	    elbow tiles (denoted /) and cross tiles (denoted +) such that no two pipes
	    cross twice.
	Example
	    # pipeDreams w == # (decompose inI)
	Text
	    In particular, to read off an associated prime of the antidiagonal initial ideal from a pipe dream,
	    one simply needs to read off the + tiles from the grid.
	Example
	    netList ((pipeDreams w)_0)
	    (decompose inI)_0
	Text
	    Diagonal term orders are less well understood for Schubert determinantal ideals and ASM ideals. 
	    Their study was initiated by [Kle23] and has also been investigated in [KW21] and [HPW22].
	    This package provides functionality for investigating three diagonal term orders:
	    One which uses lex and orders the variables diagonally reading from the southeast corner @TO diagLexInitSE@,
	    one which uses lex and orders the variables diagonaly reading from the northwest corner @TO diagLexInitNW@,
	    and one which uses revlex on the variables diagonally reading from the northwest corner @TO diagRevLexInit@.
	Text
	    Two diagonal term orders may give two distinct initial ideals which are not squarefree.
	Example
	    v = {2,1,4,3,6,5}
	    diagLexInitSE v
	    netList (decompose oo)
	    diagLexInitNW v
	    netList (decompose oo)
	Text
	    In this case, @TO diagRevLexInit@ and @TO diagLexInitSE@ give the same initial ideal.
	    It is unknown if this is the case in general.
	Example
	    diagRevLexInit v
	Text
	    @SUBSECTION "Functions for studying initial ideals of ASM ideals"@
	Text
	    @UL {
		TO (antiDiagInit, Matrix),
		TO (diagLexInitSE, Matrix),
		TO (diagLexInitNW, Matrix),
		TO (diagRevLexInit, Matrix),
		TO (pipeDreams, List)
		}@
	
///

doc ///
    Key 
        "Functions for investigating permutations"
    Headline
        basic functions for permutations
    Description
        Text 
            This package provides significantly expanded functionality for studying permutations
	    in Macaulay2.
	Text
	    @UL {
	    {"[HPW22] Zachary Hamaker, Oliver Pechenik, and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2003.13719", EM "Grobner geometry of Schubert polynomials through ice"),
	    " , Advances in Mathematics 398 (2022): 108228."},
	    {"[Kle23] Patricia Klein, ",
	    HREF("https://arxiv.org/abs/2008.01717", EM "Diagonal degenerations of matrix Schubert varieties"),
	    " , to appear in Algebraic Combinatorics."},	
            {"[KM05] Allen Knutson and Ezra Miller, ",
            HREF("https://arxiv.org/abs/math/0110058", EM "Grobner geometry of Schubert polynomials"),
            " , Annals of Mathematics (2005): 1245-1318."},
	    {"[KW21] Patricia Klein and Anna Weigandt, ",
	    HREF("https://arxiv.org/abs/2108.08370", EM "Bumpless pipe dreams encode Grobner geometry of Schubert polynomials"),
	    " , arxiv preprint 2108.08370."}
	    }@
	Text
	    Given a permutation as a list of integers, one can check if it is indeed a permutation,
	    what its descent set is, what its inverse is, and the Coxeter length of the permutation.
	Example
	    v = {2,1,6,3,5,4};
	    isPerm v --checks if v is indeed a permutation
	    lastDescent v
	    firstDescent v
	    inverseOf v
	Text
	    This package also allows one to quickly compute certain combinatorial polynomials
	    associated to a permutation, such as the (double) Schubert polynomial and the
	    Grothendieck polynomial.
	Example
	    u = {3,1,4,2}
	    schubertPoly u
	    doubleSchubertPoly u
	    grothendieckPoly u
	Text
	    Moreover, this package contains functionality for checking whether
	    a permutation avoids a set of patterns.
	    For instance, @TO isCDG@ checks whether a permutation is
	    CDG; @TO isVexillary@ checks whether a permutation
	    is 2143-avoiding; and @TO isCartwrightSturmfels@ checks whether
	    a permutation is Cartwright-Sturmfels.
	Example
	    w = {1,2,3,9,8,4,5,6,7};
	    isPatternAvoiding(w,{4,1,2,3})
	    isVexillary w
	    isCartwrightSturmfels w
	    isCDG w
	Text
	    Finally, this package contains functionality for studying
	    both reduced and nonreduced pipe dreams of a permutation.
	Example
	    decompose antiDiagInit u
	    apply(pipeDreams u, i -> netList i)	    
	    apply(pipeDreamsNonReduced u, i-> netList i)
	Text
	    @SUBSECTION "Functions for studying permutations"@
	Text
	    @UL {
		TO (isPerm, List),
		TO (permToMatrix, List),
		TO (lastDescent, List),
		TO (firstDescent, List),
		TO (permLength, List),
		TO (inverseOf, List),
		TO (longestPerm, ZZ),
		TO (toOneLineNotation, List, ZZ),
		TO (composePerms, List, List),
		TO (isPatternAvoiding, List, List),
		TO (isVexillary, List),
		TO (avoidsAllPatterns, List, List),
		TO (isCartwrightSturmfels, List),
		TO (isCDG, List),
		TO (rajCode, List),
		TO (rajIndex, List),
		TO (grothendieckPoly,List),
		TO (schubertPoly, List),
		TO (doubleSchubertPoly, List),
		TO (pipeDreams, List),
		TO (pipeDreamsNonReduced, List),
		TO (netPD, List)
		}@
///

doc ///
    Key
        isPartialASM
        (isPartialASM, Matrix)
    Headline
        whether a matrix is a partial alternating sign matrix
    Usage
        isPartialASM(M)
    Inputs
        M:Matrix
    Outputs
    	:Boolean
    Description
        Text
            Given an integer matrix, checks that the matrix is a partial alternating sign matrix. A partial alternating sign matrix is a matrix with entries in $\{-1,0,1\}$ such that:

            - The nonzero entries in each row and column alternate in sign,

            - The first nonzero entry of any row or column (if there is one) is $1$.
        Example
            M = matrix{{0,0,1,0,0,0,0,0},{1,0,1,0,1,0,0,0},{0,0,0,1,-1,0,0,1},{0,0,1,-1,1,0,0,0},{0,0,0,0,0,0,1,0},{0,0,0,0,0,1,0,0},{0,1,-1,1,0,0,0,0},{0,0,1,0,0,0,0,0}}
            isPartialASM M
            N = matrix{{0,-1,0,1,1},{1,-1,1,-1,1},{0,1,1,0,-1},{1,1,-1,1,-1},{-1,1,0,0,1}}
            isPartialASM N
///

doc ///
    Key
        partialASMToASM
        (partialASMToASM, Matrix)
    Headline
        extend a partial alternating sign matrix to an alternating sign matrix
    Usage
        partialASMToASM(A)
    Inputs
        A:Matrix
    Outputs
	:Matrix
    Description
        Text
            Given a partial alternating sign matrix, returns the unique smallest alternating sign matrix with the same essential set and same rank conditions at each element of the essential set.
        Example
            A = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            partialASMToASM(A)

	    
///

doc ///
    Key
        antiDiagInit
	(antiDiagInit, List)
        (antiDiagInit, Matrix)
	[antiDiagInit, CoefficientRing]
	[antiDiagInit, Variable]
    Headline
        compute the (unique) antidiagonal initial ideal of an ASM ideal
    Usage
        antiDiagInit(w)
        antiDiagInit(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
        CoefficientRing => Ring
	Variable => Symbol
    Outputs
    	:MonomialIdeal
    Description
        Text
            Let $Z = (z_{i,j})$ be a generic matrix and $R=k[Z]$ a polynomial ring in the entries of $Z$ over the field $k$.  We call a term order on $R$ antidiagonal if the lead term of the determinant of each submatrix $Z'$ of $Z$ is the product of terms along the antidiagonal of $Z'$. 

            This method computes the antidiagonal initial ideal of an ASM ideal by directly forming the ideal of the lead terms of the Fulton generators.  

            @UL {{"[KM05]: Knutson and Miller, Gröbner geometry of Schubert polynomials (see ", arXiv "0110058", ")."},}@ 
            
            tells us that the Fulton generators of each Schubert determinantal ideal form a Gröbner basis.  For an extension to ASM ideals, see  
	    
	    @UL {{"[KW]: Klein and Weigant, Bumpless pipe dreams encode Gröbner geometry of Schubert polynomials (see ", arXiv "2108.08370", ")."},}@
            
            @UL {{"[Wei]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
            
            @UL {{"[Knu]: Knutson, Frobenius splitting, point-counting, and degeneration (see ", arXiv "0911.4941", ")."},}@
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
           
        Example
            antiDiagInit({1,3,2},CoefficientRing=>ZZ/3001)
            antiDiagInit(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}, Variable => t)
///

doc ///
    Key
        rankTable
        (rankTable, List)
        (rankTable, Matrix)
    Headline
        compute a table of rank conditions that determines a Schubert determinantal ideal or, more generally, an alternating sign matrix ideal.
    Usage
        rankTable(w)
        rankTable(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
	:Matrix
    Description
        Text
            Given an alternating sign matrix or a permutation in 1-line notation, outputs the matrix of rank conditions associated to that alternating sign matrix or permutation.
        Example
            rankTable({1,3,2})
            rankTable(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        rotheDiagram
        (rotheDiagram, List)
        (rotheDiagram, Matrix)
    Headline
        find the Rothe diagram of a partial alternating sign matrix
    Usage
        rotheDiagram(w)
        rotheDiagram(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given a permutation in 1-line notation or a partial alternating sign matrix returns the Rothe diagram.
        Example
            w = {2,5,4,1,3}
            rotheDiagram(w)
            A = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            rotheDiagram(A)

	    
///

doc ///
    Key
        augmentedRotheDiagram
        (augmentedRotheDiagram, List)
        (augmentedRotheDiagram, Matrix)
    Headline
        find the Rothe diagram of a partial alternating sign matrix together with the rank table determining the alternating sign matrix variety
    Usage
        augmentedRotheDiagram(w)
        augmentedRotheDiagram(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given a permutation in 1-line notation or a partial alternating sign matrix returns the Rothe diagram (as a list) with the rank conditions at each diagram box.
        Example
            w = {2,5,4,1,3}
            augmentedRotheDiagram(w)
            A = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            augmentedRotheDiagram(A)
///


doc ///
    Key
        essentialSet
        (essentialSet, List)
        (essentialSet, Matrix)
    Headline
        compute the essential set in the Rothe Diagram for a partial alternating sign matrix or a permutation.
    Usage
        essentialSet(w)
        essentialSet(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given an alternating sign matrix or a permutation in 1-line notation, outputs Fulton's essential set, i.e., the maximally southeast elements of the Rothe diagram of that alternating sign matrix or permutation. 
        Example
            essentialSet({1,3,2})
            essentialSet(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        augmentedEssentialSet
        (augmentedEssentialSet, List)
        (augmentedEssentialSet, Matrix)
    Headline
        find the essential set of a partial alternating sign matrix or a permutation together with the rank conditions determining the alternating sign matrix variety
    Usage
        augmentedEssentialSet w
        augmentedEssentialSet A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given a permutation in 1-line notation or a partial alternating sign matrix, returns the essential set together with the rank condition at each element of the set.
        Example
            w = {2,5,4,1,3}
            augmentedEssentialSet w
            A = matrix{{0,1,0},{1,-1,0},{0,0,0}}
           augmentedEssentialSet A
///


doc ///
    Key
        schubDetIdeal
        (schubDetIdeal, List)
        (schubDetIdeal, Matrix)
	[schubDetIdeal, CoefficientRing]
	[schubDetIdeal, Variable]
    Headline
        compute an alternating sign matrix ideal (for example, a Schubert determinantal ideal)
    Usage
        schubDetIdeal w
        schubDetIdeal A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:Ideal
    Description
        Text
            Given a permutation in 1-line notation or, more generally, a partial alternating sign matrix, outputs the associated alternating sign matrix ideal (which is called a Schubert determinantal ideal in the case of a permutation).  (The convention throughout this package is that the permutation matrix of a pemutation w has 1's in positions (i,w(i)).)
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
        Example
            schubDetIdeal({1,3,2},CoefficientRing=>ZZ/3001)
            schubDetIdeal(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}, Variable => y)

///


doc ///
    Key 
        fultonGens
        (fultonGens, List)
        (fultonGens, Matrix)
	[fultonGens, CoefficientRing]
	[fultonGens, Variable]
    Headline
        compute the Fulton generators of an ASM ideal (for example, a Schubert determinantal ideal)
    Usage
        fultonGens w
        fultongens A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
        CoefficientRing => Ring
	Variable => Symbol
    Outputs
    	:List
    Description
        Text
            Given a partial alternating sign matrix or permutation in 1-line notation, returns the list of Fulton generators for the corresponding Schubert determinantal ideal or, more generally, ASM ideal.
        Example 
            netList fultonGens({2,5,4,1,3}, CoefficientRing => ZZ/101)
            netList fultonGens(matrix{{0,1,0},{1,-1,1},{0,1,0}}, Variable => a)
///

doc ///
    Key 
        diagLexInitSE
    	(diagLexInitSE, List)
	(diagLexInitSE, Matrix)
        [diagLexInitSE, CoefficientRing]
	[diagLexInitSE, Variable]
    Headline
        Diagonal initial ideal of an ASM ideal with respect to lex, starting from SE corner
    Usage
    	diagLexInitSE w
	diagLexInitSE A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
	CoefficientRing => Ring
	Variable => Symbol
    Outputs
    	:MonomialIdeal
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, return the diagonal initial ideal of the corresponding ASM ideal or Schubert determinantal ideal with respect to lexicographic order, 
	    where the variables are ordered reading from right to left and bottom-to-top (starting in the southeast corner).
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
	Example
	    diagLexInitSE({1,3,2},CoefficientRing=>ZZ/3001)
            diagLexInitSE(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}, Variable => c)
///


doc ///
    Key 
	diagLexInitNW
        (diagLexInitNW, List)
	(diagLexInitNW, Matrix)
	[diagLexInitNW, CoefficientRing]
	[diagLexInitNW, Variable]
    Headline
        Diagonal initial ideal of an ASM ideal with respect to lex, starting from NW corner
    Usage
    	diagLexInitNW w
	diagLexInitNW A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
        CoefficientRing => Ring
	Variable => Symbol  
    Outputs
    	:MonomialIdeal  
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, returns the diagonal initial ideal of the corresponding ASM ideal or Schubert determinantal ideal with respect to lexicographic order, 
	    where the variables are ordered reading from left-to-right and top-to-bottom (starting in the northwest corner).
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
	Example
	    diagLexInitNW({1,3,2},CoefficientRing=>ZZ/3001)
            diagLexInitNW(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///


doc ///
    Key 
        diagRevLexInit
    	(diagRevLexInit, List)
	(diagRevLexInit, Matrix)
	[diagRevLexInit, CoefficientRing]
	[diagRevLexInit, Variable]
	CoefficientRing
	Variable
    Headline
        Diagonal initial ideal of an ASM ideal with respect to revlex, ordering variables from NW corner 
    Usage
    	diagRevLexInit w
	diagRevLexInit A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
	CoefficientRing => Ring
	Variable => Symbol    
    Outputs
    	:MonomialIdeal 
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, return the diagonal initial ideal of the corresponding ASM ideal or Schubert determinantal ideal with respect to reverse lexicographic order, 
	    where the variables are ordered reading from left-to-right and top-to-bottom (starting in the northwest corner).
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
	Example
	    diagRevLexInit({1,3,2},CoefficientRing=>ZZ/3001)
            diagRevLexInit(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        subwordComplex    
        (subwordComplex, List)
    Headline
        to find the subword complex associated to w
    Usage
        subwordComplex w
    Inputs
        w:List
    Outputs
    	:SimplicialComplex
    Description
        Text
            Given a permutation in 1-line notation, compute the subword complex associated to w (that is, the Stanley-Reisner complex of antiDiagInit).
        Example
            subwordComplex({2,5,4,1,3})
/// 

doc ///
    Key
        entrywiseMinRankTable
        (entrywiseMinRankTable, List)
    Headline
        compute the entrywise minimum rank table of a list of ASMs
    Usage
        entrywiseMinRankTable(L)
    Inputs
        L:List
            of ASMs of equal size
    Outputs
        :Matrix
    Description
        Text
            Computes the rank tables of a list of ASMs, then returns the entrywise minimum.
        Example
            L = {{4,3,1,2},{2,4,3,1}} / permToMatrix;
            entrywiseMinRankTable L
/// 

doc ///
    Key
        entrywiseMaxRankTable
        (entrywiseMaxRankTable, List)
    Headline
        compute the entrywise maximum rank table of a list of ASMs
    Usage
        entrywiseMaxRankTable L
    Inputs
        L:List
            of ASMs of equal size
    Outputs
        :Matrix
    Description
        Text
            Computes the rank tables of a list of ASMs, then returns the entrywise maximum.
        Example
            L = {{4,3,1,2},{2,4,3,1}} / permToMatrix;
            entrywiseMaxRankTable L
///         


doc ///
    Key
        schubDecomposition
        (schubDecomposition, Ideal)
	(schubDecomposition, Matrix)
    Headline
        finds the decomposition of an ASM ideal into Schubert determinantal ideals
    Usage
        schubDecomposition I
    Inputs
        I:Ideal
	    or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given an ASM ideal $I_A$, it can be decomposed into Schubert determinantal ideals as $I_A = I_{w_1} \cap ... \cap I_{w_k}$, where the $w_i$ are permutations.
            As output, each element in the list is the permutation associated to a prime component in the Schubert decomposition of the antidiagonal initial ideal of $I$.
	Example
	    A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
	    J = schubDetIdeal A;
	    netList schubDecomposition J
	Text
	    If the ASM ideal for an ASM $A$ has not het been computed, one may also give the ASM $A$ as input.
	Example
	    A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}};
	    netList schubDecomposition A	    
///

doc ///
    Key
        permOverASM    
        (permOverASM, Matrix)
    Headline
        finds the permutation set of an alternating sign matrix
    Usage
        permOverASM A
    Inputs
        A:Matrix
    Outputs
    	:List    	
    Description
        Text
            Given an alternating sign matrix $A$, this routine computes Perm$(A) = \{w \in S_n \mid A \leq w$, and $v \in S_n$ with $A \leq v \leq w$ implies $ v=w\}$ (where $\leq$ is in (strong) Bruhat order).  This computation is performed by taking the antidiagonal initial ideal determined by $A$ and extracting the permutations indexing its components via schubDecomposition.
	Example 
	    A = matrix{{0,1,0,0},{0,0,1,0},{1,-1,0,1},{0,1,0,0}}
	    permOverASM A
///

doc ///
    Key
        isIntersectionSchubIdeals
        (isIntersectionSchubIdeals, Ideal)
    Headline
        whether an ideal is the intersection of Schubert determinantal ideals
    Usage
        isIntersectionSchubIdeals I
    Inputs
        I:Ideal
    Outputs
    	:Boolean
    Description
        Text
            Checks if the input ideal $I$ can be written as $I = I_{w_1} \cap ... \cap I_{w_k}$,
            where each $I_{w_i}$ is a Schubert determinantal ideal. 
	    
	    This function computes the antidiagonal initial ideal of $I$ (using the default term order in Macaulay2, which is antidiagonal), finds the primes in the decomposition of $I$, reads a permutation from each such prime, and checks if $I$ is the intersection of the Schubert determinantal ideals of those permutations.
	    
	    The following theorems combine to guarantee that, if $I$ can be written as the intersection of Schubert determinantal ideals, it is exactly the intersection of the Schubert determinantal ideals found by the algorithm described above.
	    
	     @UL {{"[KM05, Theorem B]: Knutson and Miller, Gröbner geometry of Schubert polynomials (see ", arXiv "0110058", ")."},}@ 
            
            @UL {{"[Wei, Proposition 5.4]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
            
            @UL {{"[BB93, Theorem 3.7]: Bergeron and Billey, RC-graphs and Schubert polynomials."},}@
	    
	Example
	    A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
	    J = schubDetIdeal A;
	    isIntersectionSchubIdeals J
///

doc ///
    Key
        isASMIdeal
        (isASMIdeal, Ideal)
    Headline
        whether an ideal is an ASM ideal
    Usage
        isASMIdeal I
    Inputs
        I:Ideal
    Outputs
    	:Boolean
    Description
        Text
            Every ASM ideal can be written as an intersection of Schubert determinantal ideals.  Given an ideal $I$, this function first uses schubDecomposition to find the set of permutations that must index the minimal primes of $I$ if indeed $I$ is an ASM ideal. Then $I$ is an ASM ideal if and only if $I=I_A$ for the ASM $A$ whose rank table is the determined by taking entrywise maxima (using entrywiseMaxRankTable) in the rank tables of the permutations found by schubDecomposition.
	
	    When this function returns true, it also stores the ASM $A$ so that $I=I_A$.  The matrix $A$ can then be accessed using getASM.
	Example
    	   I1=schubDetIdeal {3,4,1,2};
	   I2=sub(schubDetIdeal {3,2,4,1},ring I1);
	   I = intersect(I1,I2);
	   isASMIdeal I
///

doc ///
    Key
        isASMUnion
        (isASMUnion, List)
    Headline
        whether the union of matrix Schubert varieties is an ASM variety 
    Usage
        isASMUnion L
    Inputs
        L:List
    Outputs
    	:Boolean    
    Description
        Text
            Given a list of permutations in 1-line notation, check whether the union of their matrix schubert varieties is an ASM variety. This function uses entrywiseMaxRankTable to construct the rank table that is the entrywise maximum of the rank tables of the input permutations. It then constructs an ASM $A$ from that rank table and uses permOverASM to check if the permutation set of $A$ is equal to the input list of permutations.
        
	    If the union of the matrix Schubert varieties of the input list of permutations is an ASM variety, it must be the ASM variety considered by this algorithm.
	Example 
            isASMUnion {{2,1,3,4},{4,2,3,1}} -- false 
            isASMUnion {{4,1,3,2},{3,4,1,2},{2,4,3,1}} -- true
///

doc ///
    Key
        getASM
        (getASM, Ideal)
    Headline
        get the ASM of an ideal (if it exists)
    Usage
        getASM I
    Inputs
        I:Ideal
    Outputs
    	:Matrix
    Description
        Text
            Given an ideal $I$, gets the ASM $A$ so that $I=I_A$, if such an ASM $A$ exists.
            If the ASM $A$ has already been computed and stored in the cache of I using isASMIdeal, then this function produces $A$ immediately. Otherwise, an attempt will be made to compute the ASM. Once the ASM is computed, it is stored in the cache of $I$.
	Example
	    A = matrix{{0,0,1,0,0},{1,0,0,0,0},{0,1,-1,1,0},{0,0,0,0,1},{0,0,1,0,0}};
	    I = schubDetIdeal A;
	    getASM I
///

doc ///
    Key
        isMinRankTable
        (isMinRankTable, Matrix)
    Headline
        whether a matrix is the canonical rank table of some partial ASM
    Usage
        isMinRankTable T
    Inputs
        T:Matrix
    Outputs
        :Boolean
    Description
        Text
            Checks whether {\tt T} is a the canonical rank table of some partial ASM.  These are the rank tables that are constructed in Section 1 of [Wei] and from which the partial ASM can be determined using equation (21) of the same paper.
        
	    @UL {{"[Wei]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
	    
	Example
            T = matrix {{0,1,1},{1,1,2},{1,2,3}}
            isMinRankTable T

            U = matrix {{1,1,1,1,1},{1,2,2,2,2},{1,2,2,2,3},{1,2,2,3,3}}
            isMinRankTable U
///

doc ///
    Key
        rankTableToASM
        (rankTableToASM, Matrix)
    Headline
        to find the ASM associated to a given rank table
    Usage
        rankTableToASM T
    Inputs
        T:Matrix
    Outputs
    	:Matrix  
    Description
        Text
            Given a square matrix that is a valid minimal rank table, returns the unique ASM associated to it. This algorithm follows
	    
	    @UL {{"[Wei, Equation (21)]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
	    
        Example
            T = matrix {{0,0,1,1},{0,1,1,2},{1,2,2,3},{1,2,3,4}}
            rankTableToASM(T)
            U = matrix {{0,0,1,1,1},{1,1,1,2,2},{1,2,2,3,3},{1,2,3,4,4},{1,2,3,4,5}}
            rankTableToASM(U)
///           

doc ///
    Key
        rankTableFromMatrix
        (rankTableFromMatrix, Matrix)
    Headline
        to find the minimal rank table from an arbitrary square integer matrix
    Usage
        rankTableFromMatrix M
    Inputs
        M:Matrix
    Outputs
    	:Matrix
    Description
        Text
            Given a square integer matrix (viewed as rank conditions to be imposed on a generic matrix), returns the unique integer matrix that both defines the same ideal of minors and also is the minimal rank table of some ASM.  See Section 1 and Equation (21) of 
	    
	    @UL {{"[Wei]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
	    
        Example
            M = matrix {{1,0,0},{0,23,24},{23,24,25}};
            rankTableFromMatrix M
///   

doc /// 
    Key
        schubIntersect
        (schubIntersect, List)
    Headline
        compute the intersection of ASM ideals
    Usage 
        schubIntersect L
    Inputs 
        L:List 
            of ASMs or permutations in 1-line notation
    Outputs
        :Ideal 
    Description
        Text
            Given a list of ASMs or permutations in 1-line notation, compute the intersection of the corresponding Schubert determinantal ideals.
        Example
            schubIntersect {{3,2,1,4}, {2,1,4,3}}
            schubIntersect {matrix {{0,1,0},{1,-1,1},{0,1,0}}, {3,2,1}}
///

doc ///
    Key
        schubAdd
        (schubAdd, List)
    Headline
        compute the sum of ASM ideals
    Usage 
        schubAdd L
    Inputs 
        L:List 
            of ASMs or permutations in 1-line notation
    Outputs
        :Ideal 
    Description
        Text
            Given a list of ASMs or permutations in 1-line notation, compute the sum of the corresponding Schubert determinantal (or ASM) ideals
	    
	    An arbitrary (finite) sum of partial ASM ideals is again a partial ASM ideal.  See 
	    
	     @UL {{"[Wei, Section 3.5]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
	    
	     @UL {{"[KW, Lemma 2.6]: Klein and Weigant, Bumpless pipe dreams encode Gröbner geometry of Schubert polynomials (see ", arXiv "2108.08370", ")."},}@
	    
	    The canonical rank table of the sum will be entrywise minimum of the rank tables of the summands.  This function computes the sum of ASM ideals by computing the individual rank tables, using entrywiseMinRankTable to find the entrywise minimum, and then constructing the partial ASM from that rank table.
	    
        Example
            schubAdd {{3,2,1,4}, {2,1,4,3}}
            schubAdd {matrix {{0,1,0},{1,-1,1},{0,1,0}}, {3,2,1}}
///
