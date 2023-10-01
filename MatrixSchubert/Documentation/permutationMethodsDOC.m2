doc ///
    Key
        (isPerm, List)
        isPerm
    Headline
        whether a list is a permutation in 1-line notation
    Usage
        isPerm w
    Inputs
        w:List
    Description
        Text
            Given a list of length $n$, checks if the entries of the permutation are the integers from 1 to $n$.
        Example
            w = {5,3,4,6,1,2}
            isPerm(w)

            v = {4,3,2}
            isPerm v
///

doc ///
    Key
        (permToMatrix, List)
        permToMatrix
    Headline
       converts a permutation in 1-line notation into a permutation matrix
    Usage
        permToMatrix w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, produces the permutation matrix with 1's in location $(i,w_i)$.
        Example
            w = {7,2,5,8,1,3,6,4}
            permToMatrix w

            v = {1,6,9,2,4,7,3,5,8}
            permToMatrix v
///

doc ///
    Key
        (lastDescent, List)
        lastDescent
    Headline
       finds the location of the last descent of a permutation
    Usage
        lastDescent w
    Inputs
        w:List
    Description
        Text
            Given a non-identity permutation in 1-line notation, finds the location of its last descent, i.e., the greatest $i$ so that $w_(i+1)<w_i$.
        Example
            w = {7,2,5,8,1,3,6,4}
            lastDescent w

            v = {1,6,9,2,4,7,3,5,8}
            lastDescent v
///

doc ///
    Key
        (firstDescent, List)
        firstDescent
    Headline
       finds the location of the first descent of a permutation
    Usage
        firstDescent w
    Inputs
        w:List
    Description
        Text
            Given a non-identity permutation in 1-line notation, finds the location of its first descent, i.e., the least $i$ so that $w_(i+1)<w_i$.
        Example
            w = {7,2,5,8,1,3,6,4}
            lastDescent(w)

            v = {1,6,9,2,4,7,3,5,8}
            lastDescent(v)
///

doc ///
    Key
	(permLength, List)
        permLength
    Headline
    	to find the length of a permutation in 1-line notation.
    Usage
        permLength w
    Inputs
    	w:List
    Description
    	Text
	 Given a permutation in 1-line notation returns the Coxeter length of the permutation.
	Example
    	    w = {2,5,4,1,3}
	    permLength w

	    
///

doc ///
    Key
	(inverseOf, List)
        inverseOf
    Headline
    	to return the inverse of a permutation in 1-line notation.
    Usage
        inverseOf w
    Inputs
    	w:List
    Description
    	Text
	 Given a permutation in 1-line notation returns the inverse of the permutation in 1-line notation.
	Example
    	    w = {2,5,4,1,3}
	    inverseOf w
///

doc ///
    Key
	(longestPerm, ZZ)
        longestPerm
    Headline
    	to return the longest permutation of length n
    Usage
        longestPerm n
    Inputs
    	n:ZZ
    Description
    	Text
	    Given an integer $n$, returns the permutation $\{n,n-1,...2,1\}$.
	Example
    	    longestPerm 7
///

doc ///
    Key
	(getOneReducedWord, List)
        getOneReducedWord
    Headline
        given a permutation in 1-line notation, finds one reduced word
    Usage
        getOneReducedWord w
    Inputs
    	w:List
    Description
    	Text
	    This is a stub.
///

doc ///
    Key
        (toOneLineNotation, List, ZZ)
        toOneLineNotation
    Headline
    	rewrites a transposition in 1-line notation
    Usage
        toOneLineNotation(perm, maxIdx)
    Inputs
    	perm:List
        maxIdx:ZZ
    Outputs
        :List
    Description
    	Text
    	    Converts a transposition $(a,b)$ to 1-line notation.
            {\tt maxIdx} is the $n$ for which to regard {\tt perm} as an 
            element of $S_n$, the symmetric group on $n$ letters.
        Example
            perm = {2,4}
            maxIdx = 5
            toOneLineNotation(perm, maxIdx)
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
    Outputs
        :List
    Description
        Text
            Computes the composition of two permutations, $u$ and $v$, as $u*v$.
            Note that the permutations must be written as a list in 1-line notation and must both be permutations of (the same) $n$ letters.
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
        (isPatternAvoiding, List, List)
        isPatternAvoiding
    Headline
        whether a permutation avoids certain patterns, e.g. $2143$-avoiding or $312$- and $231$-avoiding
    Usage
        isPatternAvoiding(w, pattern)
    Inputs
        w:List
        pattern:List
    Outputs
        :Boolean
    Description
        Text
            Given a permutation, checks if the permutation is pattern-avoiding, e.g. $2143$-avoiding or $1432$-avoiding.
            For example, a permutation $w$ is $2143$-avoiding if there does not exist indices $i < j < k < l$
            such that $w_j < w_i < w_l < w_k$.
        Example
            w = {7,2,5,8,1,3,6,4}
            pattern2143 = {2,1,4,3}
            isPatternAvoiding(w, pattern2143)

            v = {2,3,7,1,5,8,4,6}
            pattern1432 = {1,4,3,2}
            isPatternAvoiding(v, pattern1432)
///


doc ///
    Key
        (isVexillary, List)
        isVexillary
    Headline
        whether a permutation is vexillary, i.e. 2143-avoiding
    Usage
        isVexillary w 
    Inputs
        w:List
    Outputs
        :Boolean
    Description
        Text
            Given a permutation in 1-line notation, checks if the permutation is vexillary, i.e. $2143$-avoiding.
            A permutation $w$ is $2143$-avoiding if there do not exist indices $i < j < k < l$
            such that $w_j < w_i < w_l < w_k$.
        Example
            w = {7,2,5,8,1,3,6,4}
            isVexillary w

            v = {1,6,9,2,4,7,3,5,8}
            isVexillary v
///

doc ///
    Key
        (avoidsAllPatterns, List, List)
        avoidsAllPatterns
    Headline
        whether a permutation avoids all of the given patterns
    Usage
        avoidsAllPatterns(perm, patterns)
    Inputs
        perm:List
        patterns:List
    Outputs
        :Boolean
    Description
        Text
            Given a permutation in one-line notation, and a list of patterns checks if the permutation avoids every pattern.
	Example 
	    w = {7,2,5,8,1,3,6,4}
	    patterns = {{2,1,4,3},{1,4,3,2}}
	    avoidsAllPatterns(w,patterns)
///

doc ///
    Key
        (isCartwrightSturmfels, List)
        isCartwrightSturmfels
    Headline
        whether a permutation is Cartwright-Sturmfels
    Usage
        isCartwrightSturmfels w
    Inputs
        w:List
    Outputs
        :Boolean
    Description
        Text
            Given a permutation in 1-line notation, checks if the permutation is Cartwright-Sturmfels.  By [CDG22], the matrix
	    Schubert variety $X_w$ is Cartwright-Sturmfels if and only if $w$ avoid all of the patterns 
	    $\{12543, 13254, 13524, 13542, 21543, 125364, 125634, 215364, 215634, 315264, 315624, 315642\}$.
	    
	     @UL {
            {"[CDG22] A. Conca, E. De Negri, and E. Gorla, ",
            HREF("https://arxiv.org/abs/2108.10115", EM "Radical generic initial ideals"),
            ", Vietnam J. Math. 50 (2022), no. 3, 807-827."}
            }@
	    
        Example
            w = {7,2,5,8,1,3,6,4}
            isCartwrightSturmfels w

            v = {1,6,9,2,4,7,3,5,8}
            isCartwrightSturmfels v
///

doc ///
    Key
        (isCDG, List)
        isCDG
    Headline
        whether a permutation is CDG
    Usage
        isCDG(perm)
    Inputs
        perm:List
    Outputs
        :Boolean
    Description
        Text
            Given a permutation in 1-line notation, checks if the permutation is CDG.  We say that a permutation $w$ is CDG 
	    if a certain modification (see [Kle23] for precise description) of the Fulton generators of the Schubert determinantal
	    ideal $I_w$ form a diagonal Grobner basis.  By [Kle23], $w$ CDG if and only if $w$ avoid all of the patterns 
	    $\{13254, 21543, 214635, 215364, 215634, 241635, 315264, 4261735\}$.
	    
	     @UL {
            {"[Kle23] P. Klein, ",
            HREF("https://arxiv.org/abs/2008.01717", EM "Diagonal degenerations of matrix Schubert varieties"),
            ", Algebr. Comb. 6 (2023), no. 4, 1073-1094."}
            }@
	    
        Example
            w = {7,2,5,8,1,3,6,4}
            isCDG w

            v = {1,6,9,2,4,7,3,5,8}
            isCDG v
///

doc ///
    Key
        (rajCode, List)
        rajCode
    Headline
      finds the Rajchgot code of a permutation
    Usage
        rajCode w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, finds its Rajchgot code, as defined in [PSW].
	    
	    @UL {{"[PWS]: O. Pechenik, A. Weigandt, and D. Speyer, \"Castlenuovo--Mumford regularity of matrix Schubert varieties\" (see ", arXiv "2111.10681", ")."},}@
	    
        Example
            w = {7,2,5,8,1,3,6,4}
	    rajCode w
	    
            v = {1,6,9,2,4,7,3,5,8}
            rajCode v
///

doc ///
    Key
        (rajIndex, List)
        rajIndex
    Headline
      finds the Rajchgot index of a permutation
    Usage
        rajIndex w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, finds its Rajchgot index, as defined in [PSW].
	    
	    @UL {{"[PWS]: O. Pechenik, A. Weigandt, and D. Speyer, \"Castlenuovo--Mumford regularity of matrix Schubert varieties\" (see ", arXiv "2111.10681", ")."},}@
	    
        Example
            w = {7,2,5,8,1,3,6,4}
	    rajIndex w
	    
            v = {1,6,9,2,4,7,3,5,8}
            rajIndex v
///





doc ///
    Key 
        grothendieckPoly
    Headline
        computes the Grothendieck polynomial of a permutation 
    Usage
        grothendieckPoly w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, finds its Grothenieck polynomial.  Two algorithms are impliemented: DividedDifference (which is the default) and PipeDream.
	    
	Example
	    w = {2,1,4,3}
	    grothendieckPoly w
	    grothendieckPoly (w,Algorithm=>"PipeDream")
	    
///

doc ///
    Key 
        schubertPoly
    Headline
        computes the (singe) Schubert polynomial of a permutation 
    Usage
        schubertPoly w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, finds its (single) Schubert polynomial.  Two algorithms are impliemented: DividedDifference (which is the default) and Transition
	    (which makes use of the transition equations for Schubert polynomials).
	 Example 
	    w = {2,1,5,4,3}
	    schubertPoly w
	    schubertPoly (w,Algorithm=>"Transition")
///

doc ///
    Key 
        doubleSchubertPoly
    Headline
        computes the double Schubert polynomial of a permutation 
    Usage
        doubleSchubertPoly w
    Inputs
        w:List
    Description
        Text
            Given a permutation in 1-line notation, finds its double Schubert polynomial.  Two algorithms are impliemented: DividedDifference (which is the default) and Transition
	    (which makes use of the transition equations for double Schubert polynomials).
	 Example 
	    w = {2,1,5,4,3}
	    doubleSchubertPoly w
	    doubleSchubertPoly (w,Algorithm=>"Transition")
///

doc ///
    Key 
        dividedDifference
    Headline
        tmp 
    Description
        Text
            This is a stub
///

doc ///
    Key 
        pipeDreams
    Headline
        tmp 
    Description
        Text
            This is a stub
///

doc ///
    Key 
        pipeDreamsNonReduced
    Headline
        tmp 
    Description
        Text
            This is a stub
///

doc ///
    Key 
        ASMToMonotoneTriangle
    Headline
        converts an ASM to a monotone triangle
    Usage
        ASMToMonotoneTriangle A
    Inputs
        A:Matrix
    Outputs
        :List
    Description
        Text
            Converts an alternating sign matrix (ASM) to a monotone triangle according to the bijection described in [HR].
            More precisely, suppose $A$ is an ASM.
            The unique monotone triangle $T=(T_0,\hdots,T_n)$ corresponding to $A$ is given by $T_m = \sum_{i=1}^m A_m $, where $A_m$ denotes the $m$th row of $A$.
            See [HR] for more details.

            @UL {{"[HR]: Z. Hamaker and V. Reiner, \"Weak Order and Descents for Monotone Triangles\" (see ", arXiv "1809.10571", ")."},}@
        Example
            A = matrix{{0,1,0,0,0,0},{0,0,0,1,0,0},{1,-1,1,-1,0,1},{0,0,0,1,0,0},{0,1,0,-1,1,0},{0,0,0,1,0,0}}
            ASMToMonotoneTriangle A
///

doc ///
    Key 
        MonotoneTriangleToASM
    Headline
        converts a monotone triangle to an ASM
    Usage
        MonotoneTriangleToASM M
    Inputs
        M:List
    Outputs
        :Matrix
    Description
        Text
            Converts an monotone triangle to an alternating sign matrix (ASM) according to the bijection described in [HR].
            More precisely, suppose $T=(T_0,\hdots,T_n)$ is an ASM.
            The unique ASM $A$ corresponding to $T$ is given by $A_m = \mathbb{1}_{T_m} - \mathbb{1}_{T_{m-1}}$, where $A_m$ denotes the $m$th row of $A$.
            See [HR] for more details.

            @UL {{"[HR]: Z. Hamaker and V. Reiner, \"Weak Order and Descents for Monotone Triangles\" (see ", arXiv "1809.10571", ")."},}@
        Example
            M = {{}, {2}, {2, 4}, {1, 3, 6}, {1, 3, 4, 6}, {1, 2, 3, 5, 6}, {1, 2, 3, 4, 5, 6}}
            MonotoneTriangleToASM M
///