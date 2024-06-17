
doc ///
	Key
		(permute, ZZdFactorization, ZZ)
	Headline
		cyclically permute the maps in a matrix factorization
	Usage
		permute(X,k)
	Inputs
		X:ZZdFactorization
		k:ZZ
	Outputs
		:ZZdfactorization
	Description
		Text
			Applies a cyclic permutation to the maps in a matrix factorization $k$ times. 
		Example
			S = QQ[x,y,z,w]
			X = ZZdfactorization{x,y,z,w}
			permute(X,2)
	SeeAlso


doc ///
	Key
		(fullCollapse, ZZdFactorization, ZZ, ZZ)
	Headline
		Converts a factorization of arbitrary period into a factorization of period $2$.  
	Usage
		fullCollapse(X,n,k)
	Inputs
		X:ZZdFactorization
			of arbitrary period d
		n:ZZ
		k:ZZ
			with $k<d$
	Outputs
		:ZZdfactorization
			of period 2
	Description
		Text
			Create a period $2$ factorization by composing k differentials ending with the $n$th differential.	
		Example
		  	S = QQ[x_1 .. x_5]
		  	X = ZZdfactorization{x_1,x_2,x_3,x_4,x_5}
		  	fullCollapse(X,2,3)
		Example
			S = QQ[x,y,t]/(t^2+t+1)
			X = ZZdfactorization{x,x,x,x}
			Y = ZZdfactorization{y,y,y,y}
			Z = dTensor(X,Y,t)
			fullCollapse(Z,1,1)
	SeeAlso


doc ///
	Key
		(collapseMF, ZZdFactorization, ZZ, ZZ)
	Headline
		reduces the period by 1  
	Usage
		collapseMF(X,k)
	Inputs
		X:ZZdFactorization
		k:ZZ
	Outputs
		:ZZdfactorization
			of period X.period-1
	Description
		Text
			Compose two maps in a period d factorization to obtain a period d-1 factorization. The integer k determines which two maps get composed, namely, the $k$th and the $(k+1)$st differential.
		Example
		  	S = QQ[x_1 .. x_5]
		  	X = ZZdfactorization{x_1,x_2,x_3,x_4,x_5}
		  	collapseMF(X,3)
	SeeAlso


doc ///
	Key
		(tailMF,Module)
	Headline
		generate a factorization of period $2$ from a module
	Usage
		tailMF(M)
	Inputs
		M:Module
	Outputs
		:ZZdfactorization
	Description
		Text

		Example
			S = QQ[x,y,t]/(t^2+t+1)
			R = S/ideal(x^3+y^3)
			M = R^1/ideal(x,y)
			tailMF(M)
	SeeAlso
///