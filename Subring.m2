-- -*- coding: utf-8 -*-
newPackage(
    "Subring",
    Version => "1.0",
    Date => "June 6, 2023",
    Authors => {
	{Name => "Casey Hill", Email => "caseybhill2@gmail.com"},
        {Name => "Trevor K. Karn", Email => "karnx018@umn.edu"},
        {Name => "Miranda Moore", Email => "moor2340@umn.edu"},
        {Name => "Christopher O'Neill", Email => "cdoneill@sdsu.edu"}},
    Headline => "a package for subrings",
    Keywords => {"Documentation"},
    DebuggingMode => true
    )

export {"Subring",
    "subring"}

protect symbol ambientRing --Take this out in the future

veryLocalVar:=3;

Subring = new Type of HashTable

-- a method to create subrings
subring = method()
subring Matrix := genMatrix -> (
    localBar:=3;
    new Subring from {
	generators => genMatrix,
	ambientRing => ring genMatrix,
	"otherKey" => 3,
	-- presentation ring: one variable for each generator
	-- presentation map: presentation ring --> ambient ring 
	-- presentation ideal?? -- probably compute into cache 
	cache => new CacheTable
	}
    )

subring List := genList -> (
    subring matrix {genList}
    )

-- beginDocumentation()
-*
doc ///
 Node
  Key
   Subring
  Headline
     an example Macaulay2 package
  Description
   Text
    {\em FirstPackage} is a basic package to be used as an example.
  Caveat
    Still trying to figure this out.
  Subnodes
    firstFunction
 Node
  Key
   (firstFunction,ZZ)
   firstFunction
  Headline
   a silly first function
  Usage
   firstFunction n
  Inputs
   n:
  Outputs
   :
    a silly string, depending on the value of {\tt n}
  Description
   Text
    Here we show an example.
   Example
    firstFunction 1
    firstFunction 0
///

TEST ///
    assert ( firstFunction 2 == "D'oh!" )
///
*-

end--

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
