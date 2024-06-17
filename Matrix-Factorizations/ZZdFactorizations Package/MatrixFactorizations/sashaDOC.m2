doc ///
    Key
    	isdFactorization
    Headline
    	check if a ZZdFactorization object is a factorization of some polynomial
    Usage
    	isdFactorization(X)
    Inputs
    	X: ZZdFactorization
    Outputs
    	A sequence (x, f), where x is true or false. If x is true, f is a ring element. If x is false, f is the string "no potential".
    Description
    	Text
	    A ZZdFactorization object X of period p may not have (X.dd)^p equal to a scalar multiple of the identity map. 
	    This function checks if this is true, and returns the polynomial if it is true.
	    
	    Here is a simple factorization of period 3.
	Example
	    S = ZZ/101[x, y, z];
	    X = ZZdfactorization{x, y, z};
	    isdFactorization X
	Text
	    Here is a ZZdFactorization which is not a factorization. We do not expect the direct sum of two different polynomials to be a d-factorization.
	Example
	    Y = ZZdfactorization{x, x, x};
	    dsum = X ++ Y;
	    (dsum.dd)^(period dsum)
	    isdFactorization(X++Y)
	Text
	    The endomorphisms of a matrix factorization form a factorization of 0.
	Example
	    Z = ZZdfactorization{x, y};
	    isdFactorization Z
	    E = End(Z)
	    isdFactorization E
    SeeAlso
    	"isdComplex"
///

doc ///
    Key
    	isZZdComplex
    Headline
    	check if a ZZdFactorization has maps composing to 0
    Usage
    	isZZdComplex(X)
    Inputs
    	X: ZZdFactorization
    Outputs
    	true or false
    Description
    	Text
	    For a ZZdFactorization X of period p, check if d^p = 0 for its differential d.
	    In the example below, isZZdComplex returns false, since X is a factorization of xyz.
        Example
	    S = ZZ/101[x,y,z];
	    X = ZZdfactorization{x, y, z};
	    isZZdComplex X
	Text
	    Taking its endomorphisms will give a factorization of 0.
	Example
	    E = End(adjoinRoot(X, t))
	    isZZdComplex E
	SeeAlso
	    "adjoinRoot"
	
///

doc ///
    Key
    	dShift
    Headline
    	apply a shift functor to a ZZ/d-factorization of period larger than 2
    Usage
    	dShift(n, X, t), dShift(n, X, w)
    Inputs
    	X, a ZZdFactorization of period larger than 2
	n:ZZ
	    specifying how far to shift by
	t:RingElement
	    which is a pth root of unity in {\tt ring X}, where p is the period of X
	w:Symbol
	    if a root of unity needs to be adjoined to {\tt ring X}
    Outputs
    	a ZZdFactorization (does something happen to the differential?)
    	
