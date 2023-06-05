--examplePackage.m2
--thanks to Thomas Yahl for this tutorial
newPackage(
    "examplePackage",
    Version=>"1.0",
    Date=>"Oct 9, 2020",
    Authors=>{
	{Name=>"Thomas Yahl",
	 Email=>"Thomasjyahl@math.tamu.edu",
	 HomePage=>"https://math.tamu.edu/~thomasjyahl"}
    },
    Headline=>"Small package to illustrate how to make packages",
    --PackageImports is for using functions from other packages in your package.
    --PackageExports is for using functions from other packages, but also
    ----loading the package into the current M2 session.
    PackageImports=>{},
    PackageExports=>{},
    --DebuggingMode set to true runs the debugger if any error occurs.
    DebuggingMode=>true
    )


--These are the functions that the users can actually use. Generally not everything
----should be exported.
export{
    "myFunction",
    "anotherFunction",
    "notAnotherFunction"
    }


McNugget = method()
McNugget(ZZ,ZZ) := ZZ => (n,m)->(
    n*m - n - m
    )


myFunction = method()
myFunction(ZZ,ZZ) := ZZ => (n,m)->(
    print("This is myFunction");
    McNugget(n,m)
    )


anotherFunction = method()
anotherFunction(Matrix) := Number => A->(
    det A
    )


notAnotherFunction = method()
notAnotherFunction(Ring) := ZZ => R->(
    numgens R
    )

notAnotherFunction(Matrix) := ZZ => A->(
    numgens target A
    )


--This is one of many was to make documentation.
beginDocumentation()
doc///
    	Key
	    	myFunction
		    (myFunction,ZZ,ZZ)
	Headline
	    	Computes the McNugget number corresponding to two integers.
	Description
	    	Text
		    	The McNugger number corresponding to two integers n and m is the largest
			positive integer that cannot be written as a positive linear
			combination of n and m.
///

doc///
    	Key
	    	anotherFunction
		    (anotherFunction,Matrix)
	Headline
	    	Compute the determinant of a matrix
	Description
	    	Text
		    	The determinant is the unique, multilinear functional
			on the rows of a matrix normalized so that the determinant
			of the identity matrix is 1.
///

doc///
    	Key 
	    	notAnotherFunction
		    (notAnotherFunction,Ring)
		    (notAnotherFunction,Matrix)
	Headline
	    	Another function
	Description
	    	Text
		    	Determines the number of generators of a ring or the number
			of rows of a matrix.
///
		
		
		
