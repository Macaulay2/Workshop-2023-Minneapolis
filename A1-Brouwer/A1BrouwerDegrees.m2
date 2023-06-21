--A1BrouwerDegrees.m2
newPackage(
    "A1BrouwerDegrees",
    Version=>"1.0",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/nikita-borisov"},
        {Name=>"Thomas Brazelton",
	 Email=>"tbraz@math.upenn.edu",
	 HomePage=>"https://www2.math.upenn.edu/~tbraz/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/frenly-espino"},
         {Name=>"Tom Hagedorn",
	 Email=>"hagedorn@tcnj.edu",
	 HomePage=>"https://hagedorn.pages.tcnj.edu/"},
        {Name=>"Zhaobo Han",
	 Email=>"zbtomhan@sas.upenn.edu",
	 HomePage=>"https://www.linkedin.com/in/zhaobo-han-77b1301a2/"},
     	{Name=>"Jordy Lopez Garcia",
	 Email=>"jordy.lopez@tamu.edu",
	 HomePage=>"https://jordylopez27.github.io/"},
        {Name=>"Joel Louwsma",
	 Email=>"jlouwsma@niagara.edu",
	 HomePage=>"https://www.joellouwsma.com/"},
        {Name=>"Andrew Tawfeek",
	 Email=>"atawfeek@uw.edu",
	 HomePage=>"https://www.atawfeek.com/"},
        {Name=>"Wern Juin Gabriel Ong",
	 Email=>"gong@bowdoin.edu",
	 HomePage=>"https://wgabrielong.github.io/"}
	},
    Headline=>"TODO",
    PackageImports=>{},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    -- types
    "GrothendieckWittClass",
    -- methods
    "baseField",
    "gwClass",
    "gwAdd",
    "gwMultiply"
    }


---------------
-- Types
---------------

-- We define GrothendieckWittClass to be a new type,
--    meant to represent the isomorphism class of a symmetric bilinear form
--    over a base field.
GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"

gwClass = method()

-- A class in GW can be constructed from a representing matrix
gwClass (Matrix) := GrothendieckWittClass => M -> (
  new GrothendieckWittClass from {
      symbol matrix => M,
      symbol cache => new CacheTable
      }
)

-- This allows us to extract the matrix from a class
matrix GrothendieckWittClass := Matrix => beta -> beta.matrix

-- Check if a constructed class is well-defined
isWellDefined GrothendieckWittClass := Boolean => beta -> (
    
    -- Returns false if a matrix isn't symmetric
    --	  Note that this will also return false if a matrix isn't square,
    --	      so we don't need another check for that.
    if not (transpose(beta.matrix) == beta.matrix) then(
	<< "-- Defining matrix is not symmetric" << endl;
	return false
	);
    
    -- Returns false if a matrix is defined over a field of characteristic two
    if char ring beta.matrix == 2 then(
	<< "-- Package does not support base fields of characteristic two" <<endl;
	return false
	);
    
    -- Returns false if the matrix has determinant zero
    if det beta.matrix == 0 then(
	<< "-- Matrix defines a degenerate bilinear form" <<endl;
	return false
        );
    
    -- Returns false if the defining matrix isn't defined over a field
    if not isField ring beta.matrix then(
	<< "-- Matrix is not defined over a field" <<endl;
	return false
	);
    
    true);

-- Method for returning the ring the matrix is defined over
baseField = method()
baseField GrothendieckWittClass := Ring => beta -> (
    if(isWellDefined(beta)===true) then (ring beta.matrix)
)

-- Method for taking a direct sum of two Grothendieck-Witt classes
gwAdd = method()

gwAdd(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    
    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b++g)
    )

-- Method for taking a tensor product of two Grothendieck-Witt classes
gwMultiply = method()

gwMultiply(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (

    -- Returns an error if the underlying fields of the two inputted symmetric bilinear forms are different
    if not ring beta.matrix === ring gamma.matrix then error "Error: these classes have different underlying fields";

    b := beta.matrix;
    g := gamma.matrix;
    
    return gwClass(b**g)
    )




-- New type of hash table called "GrothendieckWittType"
-- GW-type.m2


---------------
-- Operations with matrices
---------------

-- Diagonalization of symmetric bilinear forms
-- diagonalize.m2


-- basic booleans about matrices
-- matrixBooleans.m2


-- Checks if a Gram matrix can easily be seen to be upper left triangular
-- easyUpperTriangular.m2

---------------
-- Operations with k-algebras
---------------

-- localAlgebraBasis.m2


---------------
-- Diagonalization over QQ
---------------

-- Takes in a rational number or integer and outputs the smallest magnitude integer in its square class
-- squarefreepart.m2

-- Given diagonal matrix, split off any <a>+<-a> and return number of times we can do this as well as smaller matrix with none of these
-- splitoffobvioushyperbolics.m2

-- Takes in symmetric matrix over QQ and diagonalizes, removes squares from entries, and splits off hyperbolic forms that immediately appear as <a> + <-a>
-- rationalsimplify.m2





-- Checks if 
-- easyIsomorphicGW
















--load();








TEST ///
    assert(1==1);
    

///


end
