
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
        (isPatternAvoiding, List, List)
        isPatternAvoiding
    Headline
        whether a permutation avoids certain pattersn, e.g. 2143-avoiding or 312- and 231-avoiding
    Usage
        isPatternAvoiding(w, pattern)
    Inputs
        w:List
        pattern:List
    Description
        Text
            Given a permutation, checks if the permutation is pattern-avoiding, e.g. 2143-avoiding or 1432-avoiding.
            For example, a permutation w is 2143-avoiding if there does not exist indices i < j < k < l
            such that w_j < w_i < w_l < w_k.
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
    Description
        Text
            Given a permutation in 1-line notation, checks if the permutation is vexillary, i.e. 2143-avoiding.
            A permutation w is 2143-avoiding if there do not exist indices i < j < k < l
            such that w_j < w_i < w_l < w_k.
        Example
            w = {7,2,5,8,1,3,6,4}
            isVexillary w

            v = {1,6,9,2,4,7,3,5,8}
            isVexillary v
///

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
            Given a list of length n, checks if the entries of the permutation are the integers from 1 to n.
        Example
            w = {5,3,4,6,1,2}
            isPerm(w)

            v = {4,3,2}
            isPerm v
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
            Given a non-identity permutation in 1-line notation, finds the location of its last descent, i.e., the greatest i so that w_(i+1)<w_i
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
            Given a non-identity permutation in 1-line notation, finds the location of its first descent, i.e., the least i so that w_(i+1)<w_i
        Example
            w = {7,2,5,8,1,3,6,4}
            lastDescent(w)

            v = {1,6,9,2,4,7,3,5,8}
            lastDescent(v)
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
            Given a permutation in 1-line notation, produces the permutation matrix with 1's in location (i,w_i).
        Example
            w = {7,2,5,8,1,3,6,4}
            permToMatrix w

            v = {1,6,9,2,4,7,3,5,8}
            permToMatrix v
///

doc ///
    Key
        (composePerms, List, List)
        composePerms
    Headline
       composes 2 permutations in 1-line notation
    Usage
        composePerms(w,v)
    Inputs
        w:List
	v:List
    Description
        Text
            Given 2 permutations w and v of the same length in 1-line notation, computes the product wv
        Example
            w = {7,2,5,8,1,3,6,4,9}
            v = {1,6,9,2,4,7,3,5,8}
            composePerms(w,v)
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
        (avoidsAllPatterns, List, List)
        avoidsAllPatterns
    Headline
        whether a permutation avoids all of the given patterns
    Usage
        avoidsAllPatterns(perm, patterns)
    Inputs
        perm:List
        patterns:List
    Description
        Text
            This is filler text.
///

doc ///
    Key
        (isCartwrightSturmfels, List)
        isCartwrightSturmfels
    Headline
        whether a permutation is CartwrightSturmfels
    Usage
        isCartwrightSturmfels(perm)
    Inputs
        perm:List
    Description
        Text
            This is filler text.
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
    Description
        Text
            This is filler text.
///


doc ///
    Key
        (rajCode, List)
        rajCode
    Headline
        to return the rajCode of a permutation in 1-line notation.
    Usage
        rajCode(perm)
    Inputs
        perm:List
    Description
        Text
            Given a permutation, returns the rajCode of the permutation.
            For the definition of rajCode see CASTELNUOVO–MUMFORD REGULARITY OF MATRIX SCHUBERT VARIETIES by OLIVER PECHENIK, DAVID E SPEYER, AND ANNA WEIGANDT, https://arxiv.org/pdf/2111.10681.
        Example
            w = {2,5,4,1,3}
            rajCode(w)
///



doc ///
    Key
        (rajIndex, List)
        rajIndex
    Headline
        to return the length of a permutation in 1-line notation.
    Usage
        rajCode(perm)
    Inputs
        perm:List
    Description
        Text
            Given a permutation, returns the rajIndex of the permutation.
            For the definition of rajIndex see CASTELNUOVO–MUMFORD REGULARITY OF MATRIX SCHUBERT VARIETIES by OLIVER PECHENIK, DAVID E SPEYER, AND ANNA WEIGANDT, https://arxiv.org/pdf/2111.10681.
        Example
            w = {2,5,4,1,3}
            rajIndex(w)
///
