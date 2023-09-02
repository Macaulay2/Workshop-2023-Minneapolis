newPackage(
         "ResolutionIdeals",
         Version => "0.1",
         Date => "May 2023",
         Headline => "Minors in Infinite Resolutions",
         Authors => {{ Name => "David Eisenbud", Email => "de@msri.or"}},
         AuxiliaryFiles => false,
         DebuggingMode => true
         )

     export {"resolutionIdeals",
	     "Size" -- option for resolutionIdeals
	     }

     -* Code section *-
      resolutionIdeals = method(Options => {Size => 1})
      resolutionIdeals (ChainComplex, ZZ) := o-> (F,r) -> (
	  --sums of r consecutive ideals of Size x Size minors of 
	  --differentials the complex F
      for i from min F to max F - r list (
	  trim sum(r, j -> minors(o.Size, F.dd_(j+1+i)))
	      )
	  )

     -* Documentation section *-
     beginDocumentation()

     doc ///
     Key
      resolutionIdeals
      (resolutionIdeals, ChainComplex, ZZ)
      [resolutionIdeals, Size]
     Headline
      list sums of ideals of minors of r consecutive differentials
     Usage
      L = resolutionIdeals(F, r)
      L = resolutionIdeals(F, r, Size => n)      
     Inputs
      F: ChainComplex
      r: ZZ
       number of ideals to sum
      Size => ZZ
       size of minors
     Outputs
      L: List
       of ideals
     Description
       Text
        Experiments by Dao and Eisenbud suggested that ideals of
	minors in consecutive differentials of a resolution often stabilize
	or become periodic of period 2, and this was proven by by 
	Brown, Dao and Sridhar in the case of resolutions of 
	modules over complete intersections and Golod rings.
        Examples of Gasharov and Peeva show 
	that this is not universal, but allow the possibility 
	that sums of 2 consecutive ideals of this type do stabilize.
	
	The following is the example of Gasharov-Peeva.
       Example
        kk = GF(13)
	
        S = kk[x_1..x_5]
        a = 2_kk -- primitive element
	apply(12, i -> a^i)
        I = ideal(
	   a*x_1*x_3+x_2*x_3,
	   (x_1+x_2)*x_4, 
	   x_3^2-x_2*x_5+a*x_1*x_5,
	   x_4^2-x_2*x_5+x_1*x_5, 
	   x_1^2,
	   x_2^2,
	   x_3*x_4,
	   x_3*x_5,
	   x_4*x_5,
	   x_5^2
	   )
        betti res(I, LengthLimit => 8)
	R = S/I
    	M = coker map(R^2, R^{2:-1}, 
	    matrix{{x_1, a*x_3+x_4},
		   {0,   x_2}})
	F = res(M, LengthLimit => 8)
	betti F
	#unique resolutionIdeals(F, 1)
	#unique resolutionIdeals (F,2)
        --example of Gasharov+Peeva
     ///

     -* Test section *-
     TEST /// -* [insert short title for this test] *-
     -- test code and assertions here
     -- may have as many TEST sections as needed
     ///

     end--

     -* Development section *-
     restart
     loadPackage ("ResolutionIdeals", Reload =>true)
     debug needsPackage "ResolutionIdeals"
     check "ResolutionIdeals"

     uninstallPackage "ResolutionIdeals"
     restart
     installPackage "ResolutionIdeals"
     viewHelp "ResolutionIdeals"

