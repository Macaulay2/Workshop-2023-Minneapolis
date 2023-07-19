
---------------------------------
---------------------------------
-- **DOCUMENTATION SECTION** --
---------------------------------
---------------------------------

doc ///
    Key
        MatrixSchubert
    Headline
        a package for investigating matrix Schubert varieties and ASM varieties
    Description
        Text
            This package provides functions for constructing and investigating matrix Schubert varieties. Many of the functions in this package can take as input either a permutation matrix in 1-line notation, or an alternating sign matrix.
        Text
            @UL {
            {"Allen Knutson and Ezra Miller, ",
            HREF("https://arxiv.org/abs/math/0110058", EM "Grobner geometry of Schubert polynomials"),
            " , Annals of Mathematics (2005): 1245-1318."},
            {"Oliver Pechenik, David Speyer, and Anna Weigandt, ",
            HREF("https://arxiv.org/abs/2111.10681", EM "Castelnuovo-Mumford regularity of matrix Schubert varieties"),
            " , arxiv preprint 2111.10681."},
            {"Anna Weigandt, ",
            HREF("https://arxiv.org/abs/1708.07236", EM "Prism tableaux for alternating sign matrix varieties"),
            " , arXiv preprint 1708.07236."}
            }@
        Example
            w = {1,5,3,4,2};
            essentialSet w   
            netList fultonGens w
        Text
            This package also contains functions for studying homological properties of ASM varieties.
        Example
            grothendieckPoly w
            betti res antiDiagInit w	  
        Text
            @SUBSECTION "Contributors"@
        Text
            The following people have generously contributed code, improved existing code, or enhanced the documentation:
            @HREF("https://sites.google.com/illinois.edu/shiliang-gao", "Shiliang Gao")@,
            @HREF("https://www.math.tamu.edu/directory/graduate.html", "Pooja Joshi")@, and
            @HREF("https://www.clemson.edu/science/academics/departments/mathstat/about/profiles/arakoto", "Antsa Tantely Fandresena Rakotondrafara")@.
///

doc ///
    Key
        (isPartialASM, Matrix)
        isPartialASM
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
        (partialASMToASM, Matrix)
        partialASMToASM
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
        (antiDiagInit, List)
        (antiDiagInit, Matrix)
        antiDiagInit
    Headline
        compute the (unique) antidiagonal initial ideal of an ASM ideal
    Usage
        antiDiagInit(w)
        antiDiagInit(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
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
            antiDiagInit(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        (rankMatrix, List)
        (rankMatrix, Matrix)
        rankMatrix
    Headline
        compute a matrix of rank conditions that determines a Schubert determinantal ideal or, more generally, an alternating sign matrix ideal.
    Usage
        rankMatrix(w)
        rankMatrix(A)
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
	:Matrix
    Description
        Text
            Given an alternating sign matrix or a permutation in 1-line notation, outputs the matrix of rank condition associated to that alternating sign matrix or permutation.
        Example
            rankMatrix({1,3,2})
            rankMatrix(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        (rotheDiagram, List)
        (rotheDiagram, Matrix)
        rotheDiagram
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
        (augmentedRotheDiagram, List)
        (augmentedRotheDiagram, Matrix)
        augmentedRotheDiagram
    Headline
        find the Rothe diagram of a partial alternating sign matrix together with the rank conditions determining the alternating sign matrix variety
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
        (essentialSet, List)
        (essentialSet, Matrix)
        essentialSet
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
        (augmentedEssentialSet, List)
        (augmentedEssentialSet, Matrix)
        augmentedEssentialSet
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
        (schubDetIdeal, List)
        (schubDetIdeal, Matrix)
        schubDetIdeal
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
            Given a permutation in 1-line notation or, more generally, a partial alternating sign matrix, outputs the associated alternating sign matrix ideal (which is called a Schubert determinantal ideal in the case of a permutation).  (The convention throughout this package is that the permuation matrix of a pemutation w has 1's in positions (i,w(i)).)
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
        Example
            schubDetIdeal({1,3,2},CoefficientRing=>ZZ/3001)
            schubDetIdeal(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///


doc ///
    Key 
        (fultonGens, List)
        (fultonGens, Matrix)
        fultonGens
    Headline
        compute the Fulton generators of an ASM ideal (for example, a Schubert determinantal ideal)
    Usage
        fultonGens w
        fultongens A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:List
    Description
        Text
            Given a partial alternating sign matrix or permutation in 1-line notation, returns the list of Fulton generators for the corresponding Schubert determinantal ideal or, more generally, ASM ideal.
        Example 
            netList fultonGens {2,5,4,1,3}
            netList fultonGens matrix{{0,1,0},{1,-1,1},{0,1,0}}
///


doc ///
    Key 
        (diagLexInitNW, List)
	(diagLexInitNW, Matrix)
	diagLexInitNW
    Headline
        Diagonal initial ideal of an ASM ideal with respect to lex, starting from NW corner
    Usage
    	diagLexInitNW w
	diagLexInitNW A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@  
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
    	(diagLexInitSE, List)
	(diagLexInitSE, Matrix)
        diagLexInitSE
    Headline
        Diagonal initial ideal of an ASM ideal with respect to lex, starting from SE corner
    Usage
    	diagLexInitSE w
	diagLexInitSE A
    Inputs
        w:List
            or {\tt A} is a @TO Matrix@
    Outputs
    	:MonomialIdeal
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, return the diagonal initial ideal of the corresponding ASM ideal or Schubert determinantal ideal with respect to lexicographic order, 
	    where the variables are ordered reading from right to left and bottom-to-top (starting in the southeast corner).
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
	Example
	    diagLexInitSE({1,3,2},CoefficientRing=>ZZ/3001)
            diagLexInitSE(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key 
    	(diagRevLexInit, List)
	(diagRevLexInit, Matrix)
        diagRevLexInit
    Headline
        Diagonal initial ideal of an ASM ideal with respect to revlex, ordering variables from NW corner 
    Usage
    	diagRevLexInit(w)
	diagRevLexInit(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@    
    Outputs
    	:MonomialIdeal 
    Description
        Text
            Given a partial alternating sign matrix or a permutation in 1-line notation, return the diagonal initial ideal of the corresponding ASM ideal or Schubert determinantal ideal with respect to reverse lexicographic order, 
	    where the variables are ordered reading from left-to-right and top-to-bottom (starting in the northwest corne).
	    
	    This function computes over the coefficient field of rational numbers unless an alternative is specified.
	Example
	    diagRevLexInit({1,3,2},CoefficientRing=>ZZ/3001)
            diagRevLexInit(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        (subwordComplex, List)
        subwordComplex
    Headline
        to find the subword complex associated to w (i.e. SR-ideal of antiDiagInit)
    Usage
        subwordComplex(w)
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
        (schubDecomposition, Ideal)
        schubDecomposition
    Headline
        finds the decomposition of an ASM ideal into Schubert determinantal ideals
    Usage
        schubDecomposition I
    Inputs
        I:Ideal
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
///

doc ///
    Key
        (permOverASM, Matrix)
        permOverASM
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
        (isIntersectionSchubIdeals, Ideal)
        isIntersectionSchubIdeals
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
        (isASMIdeal, Ideal)
        isASMIdeal
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
        (isASMUnion, List)
        isASMUnion
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
        (getASM, Ideal)
        getASM
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
        (isMinRankTable, Matrix)
        isMinRankTable
    Headline
        whether a matrix is a valid rank table
    Usage
        isMinRankTable(A)
    Inputs
        A:Matrix
    Outputs
        :Boolean
    Description
        Text
            Checks whether {\tt A} is a valid rank table.
        Example
            T = matrix {{0,1,1},{1,1,2},{1,2,3}}
            isMinRankTable T

            T = matrix {{1,1,1,1,1},{1,2,2,2,2},{1,2,2,2,3},{1,2,2,3,3},{1,2,3,3,3}}
            isMinRankTable T
///

doc ///
    Key
        (rankTableToASM, Matrix)
        rankTableToASM
    Headline
        to find the ASM associated to a given rank table
    Usage
        rankTableToASM M
    Inputs
        M:Matrix
    Outputs
    	:Matrix  
    Description
        Text
            Given a square matrix which is a valid minimal rank table, returns the unique ASM associated to it.
        Example
            A = matrix {{0,0,1,1},{0,1,1,2},{1,2,2,3},{1,2,3,4}}
            rankTableToASM(A)
            B = matrix {{0,0,1,1,1},{1,1,1,2,2},{1,2,2,3,3},{1,2,3,4,4},{1,2,3,4,5}}
            rankTableToASM(B)
///           

doc ///
    Key
        (rankTableFromMatrix, Matrix)
        rankTableFromMatrix
    Headline
        to find the minimal rank table which is associated to a unique ASM, given a square integer matrix
    Usage
        rankTableFromMatrix M
    Inputs
        M:Matrix
    Outputs
    	:Matrix
    Description
        Text
            Given a square integer matrix, return the associated valid minimal rank table which is associated to a unique ASM.
        Example
            A = matrix {{1,0,0},{0,23,24},{23,24,25}};
            rankTableFromMatrix A
///   

doc ///
    Key
        (getPermFromASM, Matrix)
        getPermFromASM
    Headline
        returns permutation associated to perumation matrix 
    Usage
        getPermFromASM(A)
    Inputs
        A:Matrix
    Outputs
        :List
    Description
        Text
	    When {\tt A} is a permutation matrix, returns the corresponding permutation. Otherwise, returns the empty permutation.
        Example
            T = matrix {{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}
            getPermFromASM T
	    
	    T = matrix{{1,0},{0,0}}
    	    getPermFromASM T
	    
            T = matrix {{0,1,0,0},{1,0,0,0},{0,0,1,0},{0,0,0,1}}
            getPermFromASM T
///

doc ///
    Key
        (entrywiseMaxRankTable, List)
        entrywiseMaxRankTable
    Headline
        compute the entrywise max rank table of a list of ASMs
    Usage
        entrywiseMaxRankTable(L)
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
        (entrywiseMinRankTable, List)
        entrywiseMinRankTable
    Headline
        compute the entrywise min rank table of a list of ASMs
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
        (schubAdd, List)
        schubAdd
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
            Given a list of ASMs or permutations in 1-line notation, compute the sum of the corresponding schubert determinantal ideals.
        Example
            schubAdd {{3,2,1,4}, {2,1,4,3}}
            schubAdd {matrix {{0,1,0},{1,-1,1},{0,1,0}}, {3,2,1}}
///

doc /// 
    Key
        (schubIntersect, List)
        schubIntersect
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
            Given a list of ASMs or permutations in 1-line notation, compute the intersection of the corresponding schubert determinantal ideals.
        Example
            schubIntersect {{3,2,1,4}, {2,1,4,3}}
            schubIntersect {matrix {{0,1,0},{1,-1,1},{0,1,0}}, {3,2,1}}
///
