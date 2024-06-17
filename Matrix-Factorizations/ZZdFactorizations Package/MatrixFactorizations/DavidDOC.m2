
doc ///
	Key
		(KoszulMF, List, RingElement)
		(KoszulMF, ideal, RingElement)
	Headline
                Create a Koszul factorization
	Usage
		koszulMF(L,f)
		koszulMF(I,f)
	Inputs
		L:List
		f:RingElement
		I:Ideal
	Outputs
		:ZZdfactorization
	Description
		Text
		Given a list  $L = \{f_1, ..., f_n\}$ and a polynomial $f$, chooses
		$f = \sum a_if_i$ and outputs the corresponding Koszul factorization $\{a_1, ... a_n\}, \{f_1, ..., f_n\}$.
		Given an ideal, it takes generators of that ideal as a list.
		Example
			S = QQ[x,y,z]
			koszulMF({x^2,y^2,z}, x^3+y^4+z^5)
		Example
			S = QQ[x,y,z]
			koszulMF(ideal(x^2,y^2,z), x^3+y^4+z^5)		
	SeeAlso



doc ///
	Key
		(EulerMF, RingElement)
	Headline
                Create the Euler factorization of a polynomial
	Usage
		eulerMF(f)
	Inputs
		f:RingElement
	Outputs
		:ZZdfactorization
	Description
		Text  Builds the koszul factorization whose input is the Jacobian ideal.  
        
		Example
			S = QQ[x,y,z]
			eulerMF(x^3+y^4+z^5)
	SeeAlso
	

doc ///
	Key
		(randomTailMF, RingElement, ZZ, ZZ, ZZ)
	Headline
                Create a random MF from a resolution
	Usage
		randomTailMF(f,m,n,b)
	Inputs
		f:RingElement
		m:integer (rank of source)
		n:integer (rank of target)
		b:integer (randomness bound)
	Outputs
		:ZZdfactorization
	Description
		Text  Takes the matrix corresponding to the cokernel of a random matrix  $m \times n$ with values in $S/f$ where $S$ is the ring where $f$ lives.
        
		Example
			S = QQ[x,y,z]
			randomTailMF(x^3+y^4+z^5, 5,7,4)
	SeeAlso

doc ///
	Key
		(randomLinearMF, ZZ, Ring)
		(randomLinearMF, ZZ, Ring, RingElement)
		(randomLinearMF, ZZ,Ring,Symbol)
	Headline
                Create a random Koszul factorization using linearMF
	Usage
		randomLinearMF(d, S)
		randomLinearMF(d, S, f)
		randomLinearMF(d, S, t)
	Inputs=
		d:integer (randomness bound)(rank of source)
		S: Ring
		f:RingElement
		t:symbol
	Outputs
		:ZZdfactorization
	Description
		Text  Takes the matrix corresponding to the cokernel of a random matrix  $m \times n$ with values in $S/f$ where $S$ is the ring where $f$ lives.
        
		Example
			S = QQ[x,y,z]
			randomLinearMF(5, S)
			
		Example
			S = QQ[x,y,z]
			randomLinearMF(5, S, x^3+y^4+z^5)
	      	Example
			S = QQ[x,y,z]
			randomLinearMF(5, S, t)
	SeeAlso
