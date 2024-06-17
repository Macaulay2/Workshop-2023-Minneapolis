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
    Headline
    Description
        Text
        Example
	Text
	Example
	Text
	Example
	Text
	SeeAlso
///
