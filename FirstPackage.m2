-- -*- coding: utf-8 -*-
newPackage(
    "Subrings",
    Version => "1.0",
    Date => "June 6, 2023",
    Authors => {
	{Name => "Casey Hill", Email => "caseybhill2@gmail.com"},
        {Name => "Miranda Moore", Email => "moor2340@umn.edu"},
        {Name => "Jane Doe", Email => "doe@math.uiuc.edu", HomePage => "http://www.math.uiuc.edu/~doe/"}},
    Headline => "an example Macaulay2 package",
    Keywords => {"Documentation"},
    DebuggingMode => false
    )

export {"firstFunction"}

firstFunction = method(TypicalValue => String)
firstFunction ZZ := String => n -> if n == 1 then "Hello World!" else "D'oh!"

beginDocumentation()

doc ///
 Node
  Key
   FirstPackage
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

end--

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
