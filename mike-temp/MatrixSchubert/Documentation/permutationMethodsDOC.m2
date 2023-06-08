
doc ///
    Key
	(permLength, List)
        permLength
    Headline
    	to find the length of a permutation in 1-line notation.
    Usage
        permLength(w)
    Inputs
    	w:List
    Description
    	Text
	 Given a permutation in 1-line notation returns the Coxeter length of the permutation.
	Example
    	    w = {2,5,4,1,3}
	    permLength(w)

	    
///


doc ///
    Key
        (isPatternAvoiding, List, List)
        isPatternAvoiding
    Headline
        whether a permutation is pattern-avoiding, e.g. 2143-avoiding or 1432-avoiding
    Usage
        isPatternAvoiding(perm, pattern)
    Inputs
        perm:List
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
        isVexillary(perm)
    Inputs
        perm:List
    Description
        Text
            Given a permutation, checks if the permutation is vexillary, i.e. 2143-avoiding.
            A permutation w is 2143-avoiding if there do not exist indices i < j < k < l
            such that w_j < w_i < w_l < w_k.
        Example
            w = {7,2,5,8,1,3,6,4}
            isVexillary(w)

            v = {1,6,9,2,4,7,3,5,8}
            isVexillary(v)
///
