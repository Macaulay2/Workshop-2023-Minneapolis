
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
            essentialBoxes w   
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
        whether a matrix is a partial alternating sign matrix.
    Usage
        isPartialASM(M)
    Inputs
        M:Matrix
    Description
        Text
            Given an integer matrix, checks that the matrix is a partial alternating sign matrix. A partial alternating sign matrix is a matrix with entries in $\{-1,0,1\}$ such that:

            - The nonzero entries in each row and column alternate in sign,

            - Each row and column sums to $0$ or $1$, and

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
        to extend a partial alternating sign matrix to an alternating sign matrix
    Usage
        partialASMToASM(M)
    Inputs
        M:Matrix
    Description
        Text
            Given a partial alternating sign matrix returns the unique smallest alternating sign matrix by adding necessary number of rows and columns.
        Example
            M = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            partialASMToASM(M)

	    
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
    Description
        Text
            Let $Z = (z_{i,j})$ be a generic matrix and $R=k[Z]$ is a polynomial ring in the entries of Z over the ring k.  We call a term order on R antidiagonal if the lead term of the determinant of each submatrix $Z'$ of Z is the product of terms along the antidiagonal of $Z'$. 

            This method relies on these theorems.  It computes the antidiagonal initial ideal of an alternating sign matrix ideal by directly forming the ideal of the lead terms of the Fulton generators.  

            @UL {{"[KM05]: Knutson and Miller, Gröbner geometry of Schubert polynomials (see ", arXiv "0110058", ")."},}@ 
            
            tells us that the Fulton generators of each Schubert determinantal ideal form a Gröbner basis, a result extended to alternating sign matrix ideals by  
            
            @UL {{"[Knu]: Knutson, Frobenius splitting, point-counting, and degeneration (see ", arXiv "0911.4941", ")."},}@ 
            
            and by 
            
            @UL {{"[Wei]: Weigandt, Prism tableaux for alternating sign matrix varieties (see ", arXiv "1708.07236", ")."},}@
        Example
            antiDiagInit({1,3,2})
            antiDiagInit(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///

doc ///
    Key
        (rankMatrix, List)
        (rankMatrix, Matrix)
        rankMatrix
    Headline
        Computes a matrix of rank conditions that determines a Schubert determinantal ideal or, more generally, an alternating sign matrix ideal.
    Usage
        rankMatrix(w)
        rankMatrix(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
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
        to find the Rothe diagram of a partial alternating sign matrix
    Usage
        rotheDiagram(w)
        rotheDiagram(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a permutation in 1-line notation or a partial alternating sign matrix returns the Rothe diagram.
        Example
            w = {2,5,4,1,3}
            rotheDiagram(w)
            M = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            rotheDiagram(M)

	    
///

doc ///
    Key
        (augmentedRotheDiagram, List)
        (augmentedRotheDiagram, Matrix)
        augmentedRotheDiagram
    Headline
        to find the Rothe diagram of a partial alternating sign matrix together with the rank conditions determining the alternating sign matrix variety
    Usage
        augmentedRotheDiagram(w)
        augmentedRotheDiagram(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a permutation in 1-line notation or a partial alternating sign matrix returns list of entries of Rothe diagram with the ranks of each entry.
        Example
            w = {2,5,4,1,3}
            augmentedRotheDiagram(w)
            M = matrix{{0,1,0},{1,-1,0},{0,0,0}}
            augmentedRotheDiagram(M)
///


doc ///
    Key
        (essentialBoxes, List)
        (essentialBoxes, Matrix)
        essentialBoxes
    Headline
        computes a list of the essential boxes in the Rothe Diagram for an alternating sign matrix M or a permutation in 1-line notation.
    Usage
        essentialBoxes(w)
        essentialBoxes(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given an alternating sign matrix or a permutation in 1-line notation, outputs a list of the essential boxes in the Rothe diagram for that matrix. 
        Example
            essentialBoxes({1,3,2})
            essentialBoxes(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})
///



-*
--TODO: figure out error in isSchubertCM doc node
doc ///
    Key
    	(isSchubertCM, Matrix)
    	(isSchubertCM, List)
    	isSchubertCM
    Headline
    	whether $R/I_A$ is Cohen-Macaulay
    Usage
    	isSchubertCM(A)
    	isSchubertCM(w)
    Inputs
		A:Matrix
			or {\tt w} is a @TO List@
	Outputs
		:Boolean
    Description
    	Text
			Given an alternating sign matrix $A$ (resp. permutation $w$),
			checks whether $R/I_A$ (resp. $R/I_w$) is Cohen-Macaulay.

			If the input is a permutation $w$, the output is always true
			since $I_w$ is a Schubert determinantal ideal, and a theorem
			of Fulton says $R/I_w$ is always Cohen-Macaulay.
        Example
            A = matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}}
			isSchubertCM(A)

			w = {1,3,2}
            isSchubertCM(w)
///
*-

doc ///
    Key
        (schubDetIdeal, List)
        (schubDetIdeal, Matrix)
        schubDetIdeal
    Headline
        Computes an alternating sign matrix ideal
    Usage
        schubDetIdeal(w)
        schubDetIdeal(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a permutation in 1-line notation or, more generally, a (partial) alternating sign matrix, outputs the associated alternating sign matrix ideal (which is called a Schubert determinantal ideal in the case of a permutation).  (The convention throughout this package is that the permuation matrix of a pemutation w has 1's in positions (i,w(i)).)
        Example
            schubDetIdeal({1,3,2})
            schubDetIdeal(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///


doc ///
    Key 
        (fultonGens, List)
        (fultonGens, Matrix)
        fultonGens
    Headline
        compute the Fulton generators for a Schubert determinantal ideal 
    Usage
        fultonGens(w)
        fultongens(M)
    Inputs
        w:List
            or {\tt M} is a @TO Matrix@
    Description
        Text
            Given a partial alternating sign matrix or permutation in 1 line notation, return the list of fulton generators for the corresponding Schubert determinantal ideal.
        Example 
            netList fultonGens {2,5,4,1,3}
            netList fultonGens matrix{{0,1,0},{1,-1,1},{0,1,0}}
///


doc ///
    Key 
        diagLexInit
    Headline
        tmp 
    Description
        Text
            This is a stub
///

doc ///
    Key 
        diagRevLexInit
    Headline
        tmp 
    Description
        Text
            This is a stub
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
    Description
        Text
            Given a permutation in 1-line notation, compute the subword complex associated to w (i.e. SR-ideal of antiDiagInit).
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
        schubDecomposition(I)
    Inputs
        I:Ideal
    Description
        Text
            Given an ASM ideal, it can be decomposed into Schubert determinantal ideals as I = I_{w_1} \cap ... \cap I_{w_k}, where the w_i are permutations.
            As output, each element in the list is the permutation associated to a prime component in the Schubert decomposition of the antidiagonal initial ideal of I.
///

doc ///
    Key
        (isIntersectionSchubIdeals, Ideal)
        isIntersectionSchubIdeals
    Headline
        whether an ideal is the intersection of Schubert determinantal ideals
    Usage
        isIntersectionSchubIdeals(I)
    Inputs
        I:Ideal
    Description
        Text
            Checks if the input ideal I can be written as $I = I_{w_1} \cap ... \cap I_{w_k}$,
            where each $I_{w_i}$ is a Schubert determinantal ideal.
///

doc ///
    Key
        (isASMIdeal, Ideal)
        isASMIdeal
    Headline
        whether an ideal is ASM
    Usage
        isASMIdeal(I)
    Inputs
        I:Ideal
    Description
        Text
            An ideal is ASM if it is radical and can be written as the intersection of Schubert determinantal ideals.
///

doc ///
    Key
        (getASM, Ideal)
        getASM
    Headline
        gets the ASM of an ideal (if it exists)
    Usage
        getASM(I)
    Inputs
        I:Ideal
    Description
        Text
            Gets the alternating sign matrix (ASM) of an ideal $I$, if it exists.
            If the ASM has not been computed yet, then an attempt will be made to compute the ASM. Once the ASM is computed, it is stored in the cache of $I$.
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
    Description
        Text
            Given a square integer matrix, return the associated valid minimal rank table which is associated to a unique ASM.
        Example
            A = matrix {{1,0,0},{0,23,24},{23,24,25}};
            rankTableFromMatrix A
///   
