newPackage(
         "RandomIdealsPlus",
         Version => "0.1",
         Date => "",
         Headline => "",
         Authors => {{ Name => "", Email => "", HomePage => ""}},
         AuxiliaryFiles => false,
         DebuggingMode => false
         )

     export {}

     -* Code section *-
mixedIdeal = (S, BList, MList) -> (randomPureBinomialIdeal(BList, S)+randomMonomialIdeal(MList,S))   	

     -* Documentation section *-
     beginDocumentation()

     doc ///
     Key
       RandomIdealsPlus
     Headline
     Description
       Text
       Tree
       Example
       CannedExample
     Acknowledgement
     Contributors
     References
     Caveat
     SeeAlso
     Subnodes
     ///

     doc ///
     Key
     Headline
     Usage
     Inputs
     Outputs
     Consequences
       Item
     Description
       Text
       Example
       CannedExample
       Code
       Pre
     ExampleFiles
     Contributors
     References
     Caveat
     SeeAlso
     ///

     -* Test section *-
     TEST /// -* [insert short title for this test] *-
     -- test code and assertions here
     -- may have as many TEST sections as needed
     ///

     end--

     -* Development section *-
     restart
     debug needsPackage "RandomIdealsPlus"
     check "RandomIdealsPlus"

     uninstallPackage "RandomIdealsPlus"
     restart
     installPackage "RandomIdealsPlus"
     viewHelp "RandomIdealsPlus"

	"randomMonomial",
	"randomArtinIdeal",
	"randomArtinDegreeIdeal",
	"randomArtinMonomialIdeal",
	"randomArtinDegreeMonomialIdeal",
	"randomArtinRing",
	"randomArtinDegreeRing",
	"randomArtinMonomialRing",
	"randomArtinDegreeMonomialRing",
	"random toric ideal"




