
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
       This package provides functions for constructing and investigating matrix Schubert varieties.
       Many of the functions in this package can take as input either a permutation matrix in 1-line notation,
       or an alternating sign matrix.
      Example
	   w = {1,5,3,4,2};
	   essentialBoxes w   
	   netList fultonGens w
      Text
       	   This package also contains functions for studying homological properties of ASM varieties.
      Example
      	  grothendieckPoly w
	  betti res antiDiagInit w	   
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
	 Given an integer matrix, checks that the matrix is a partial alternating sign matrix.
	 A partial alternating sign matrix is a matrix with entries in $\{-1,0,1\}$ such that:
	 
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
    	(schubertDetIdeal, Matrix)
    	schubertDetIdeal
    Headline
    	Computes Schubert determinantal ideal for a given permutation.
    Usage
    	schubertDetIdeal(M)
    Inputs
    	M:Matrix
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs the Schubert determinantal ideal associated to that matrix.
	Example
	 schubertDetIdeal({1,3,2})
	 schubertDetIdeal(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///

doc ///
    Key
        (composePerms, List, List)
        composePerms
    Headline
        computes the composition of two permutations
    Usage
        composePerms(u,v)
    Inputs
        u:List
        v:List
    Description
        Text
            Computes the composition of two permutations, u and v, as u*v.
            Note that the permutations must be written as a list in one-line notation.
        Example
            u = {2,3,4,1}
            v = {4,3,2,1}
            composePerms(u,v)

            u = {1,2,3,4,5}
            v = {3,5,2,1,4}
            composePerms(u,v)

            u = {3,5,2,1,4}
            v = {1,2,3,4,5}
            composePerms(u,v)
///


doc ///
    Key
        (schubertDecomposition, Ideal)
        schubertDecomposition
    Headline
        finds the decomposition of an ASM ideal into Schubert determinantal ideals
    Usage
        schubertDecomposition(I)
    Inputs
        I:Ideal
    Description
        Text
            Given an ASM ideal, it can be decomposed into Schubert determinantal ideals
            as I = I_{w_1} \cap ... \cap I_{w_k}, where the w_i are permutations.
            As output, each element in the list is the permutation associated 
            to a prime component in the Schubert decomposition of the antidiagonal 
            initial ideal of I.
///

doc ///
    Key
        (isIntersectionSchubIdeals, Ideal)
        isIntersectionSchubIdeals
    Headline
        whether an ideal is ASM
    Usage
        isIntersectionSchubIdeals(I)
    Inputs
        I:Ideal
    Description
        Text
            An ideal I is ASM if I is radical and I = I_{w_1} \cap ... \cap I_{w_k},
            where the I_{w_i} are Schubert determinantal ideals.
///

doc ///
    Key
    	rankMatrix
	(rankMatrix, Matrix)
	(rankMatrix, List)
    Headline
    	Computes a matrix of rank conditions that determines a Schubert determinantal ideal or, more generally, an alternating sign matrix ideal.
    Usage
    	rankMatrix(M)
	rankMatrix(w)
    Inputs
    	M:Matrix
	w:List
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs the matrix of rank condition associated to that alternating sign matrix or permutation.
	Example
	 rankMatrix({1,3,2})
	 rankMatrix(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///

doc ///
    Key
    	essentialBoxes
	(essentialBoxes, Matrix)
	(essentialBoxes, List)
    Headline
    	Computes a list of the essential boxes in the Rothe Diagram for an alternating sign matrix M or a permutation in 1-line notation.
    Usage
    	essentialBoxes(M)
	essentialBoxes(w)
    Inputs
    	M:Matrix
	w:List
    Description
    	Text
	 Given an alternating sign matrix or a permutation in 1-line notation, 
	 outputs a list of the essential boxes in the Rothe diagram for that matrix. 
	Example
	 essentialBoxes({1,3,2})
	 essentialBoxes(matrix{{0,0,0,1},{0,1,0,0},{1,-1,1,0},{0,1,0,0}})

///
