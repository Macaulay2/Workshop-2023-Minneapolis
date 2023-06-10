--A1BrouwerDegrees.m2
newPackage(
    "A1BrouwerDegrees",
    Version=>"0.1",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>""},
        {Name=>"Thomas Brazelton",
	 Email=>"tbraz@math.upenn.edu",
	 HomePage=>"https://www2.math.upenn.edu/~tbraz/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>""},
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
	 HomePage=>"https://wgabrielong.github.io/"},
	},
    Headline=>"",
    PackageImports=>{},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    -- methods
    
    }


---------------
-- Types
---------------


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
















load();








TEST ///
    assert(1==1);
    

///


end
