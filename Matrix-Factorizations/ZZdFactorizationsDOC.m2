doc ///
    Key
        ZZdFactorizations
    Headline
        A package for creating and computing objects in the category of ZZ/d-graded factorizations, such as matrix factorizations
    Description
        Text
            A ZZ/d-graded factorization F  of a ring element f is a ZZ/d-graded complex of free R-modules equipped with a degree -1 (mod d) endomorphism d^F such that (d^F)^d = f * id_F. 
	    In practice, a ZZdFactorization may be visualized as a sequence of R-module maps:
	    
	    $F_0 \leftarrow F_1 \leftarrow \cdots \leftarrow F_{d-1}$
	    
	    with the caveat that $d^F_0 : F_0 \to F_{d-1}$, since one should count degree modulo d. Any 2-periodic complex may be reinterpretted as a ZZ/2-graded factorization of 0, and 
	    likewise a matrix factorization of a ring element f is equivalently a ZZ/2-graded factorization of f. 
        Example
            S = simplicialModule(ZZ^2,3,Degeneracy => true) --the integer 3 specifies a top degree, the option Degeneracy specifies whether or not to compute degeneracy maps
	    S.dd
	    S.ss
	Text
	    In general, simplicial objects are infinite objects. Because of this, the user can specify a top degree for the resulting simplicial object. If no top degree is specified, then
	    the default top degree is given by the length of the input (viewed as a chain complex). The category of simplicial R-modules is equivalent to the category of nonegatively graded
	    chain complexes via an equivalence known as the Dold-Kan correspondence. 
	    
	    This means that there is a functor that converts a chain complex into a simplicial object, known as the Dold-Kan functor. In practice, the image of the Dold-Kan functor is
	    highly nontrivial to compute by hand. However, this package uses an algorithm of Sakuranath and Kock to compute Dold-Kan images quite efficiently by using the simplicialModule
	    command:
	Example
	    R = ZZ/101[x_1..x_3];
	    K = koszulComplex vars Q
	    simplicialModule(K) --defaults to top degree 3
	    simplicialModule(K,6) --specify top degree 6
	Text
	    The other piece of the Dold-Kan correspondence is the functor that converts a simplicial module into a chain complex. This functor is often referred to as the normalization functor,
	    and is also implemented as the command normalize: 
	Example
	    Kn = normalize simplicialModule(K,4)
	    Kn.dd
	    K == Kn
	Text
	    As the above example shows, applying the normalization to the Dold-Kan image recovers the original complex.
	Text
	    @SUBSECTION "Canonical Extensions of Functors to the Category of Chain Complexes"
	Text
	    One of the main utilities of the Dold-Kan correspondence from the perspective of a commutative algebraist is as a method of canonically extending any
	    endofunctor of R-modules to an endofunctor of nonnegatively chain complexes. This is because endofunctors of R-modules naturally extend to endofunctors
	    of simplicial R-modules by just applying the functor degree-wise. Since any nonnegatively graded chain complex can be converted into a simplicial R-module
	    (and vice versa) this yields a canonical way to extend endofunctors of R-modules to chain complexes in a way that preserves homotopy equivalence. 
	    
	    This extension is often called the Dold-Puppe extension, and was introduced in CITE, where the authors did pioneering work on the theory of derived nonlinear functors.
	    In general, the Dold-Puppe extension of a functor looks quite different from the classically defined functor, and in some cases may not be quasi-isomorphic!
	Example
	    Q = ZZ/101[x_1,x_2];
	    K1 = complex {matrix{{x_1}}};
	    K2 = complex {matrix{{x_2}}};
	    T1 = K1**K2
	    T1.dd
	    T2 = simplicialTensor({K1,K2})
	    T2.dd
	    phi1 = extend(T1,T2,id_(T1_0));
	    phi2 = extend(T2,T1,id_(T1_0));
	    phi1*phi2 == id_T1
	    isNullHomotopic(phi2*phi1 - id_T2)
	Text
	    The above example allows us to directly verify that the classically defined tensor product is homotopy equivalent to the Dold-Kan version. This is true in general,
	    but there are other functors such as the symmetric/exterior power functors that have classical definitions for chain complexes that are not necessarily
	    homotopy equivalent to the Dold-Puppe extensions:
	Example
	    needsPackage "ChainComplexOperations"
	    Q = ZZ/2[x_1..x_3];
	    K = koszulComplex vars Q;
	    W1 = complex wedge2(chainComplex K) --the classically defined exterior power functor
	    W2 = ext(2,K) --the simplicial version of the exterior power functor
	    prune HH W1
	    prune HH W2
	Text
	    Notice that the two definitions of exterior power yield complexes that are not even quasi-isomorphic, and hence certainly not homotopy equivalent. Moreover, the "naive"
	    definition yields a complex that does not even have finite length homology, while the Dold-Puppe extension applied to a complex with finite lenth homology is guaranteed
	    to have finite length homology. This preservation of "good behavior" for Dold-Puppe extensions of functors was the main motivation of the authors of CITE for using
	    simplicial techniques to prove the Total Rank Conjecture in characteristic 2. The ext command along with Schur functors in general have also been implemented in a 
	    functorial way and can take morphisms of chain complexes as inputs:
	Example
	    Q=ZZ/101[x_1,x_2];
	    K = koszulComplex vars Q;
	    F = complex res ((ideal vars Q)^2);
	    phi = extend(K,F,id_(K_0))
	    ext(2,phi)
	    ephi = ext(3,phi,TopDegree => 4);
	    ephi_1
	    ephi_2
	    schurMap({2},phi) --the second symmetric power
	    sphi = schurMap({2,1},phi,TopDegree => 3);  --schur functor corresponding to partition {2,1} applied to phi
	    sphi_1
	Text    
	    This package thus lays the groundwork for computing with Dold-Puppe extensions of nonlinear functors, deriving nonlinear functors, and also allows any endofunctor of R-modules
	    implemented in Macaulay2 to immediately be canonically extended to an endofunctor at the level of chain complexes. 
	Text
	    @SUBSECTION "Contributors"
	Text
	    The following people have generously contributed code or improved existing code:
///
