-*
   Copyright 2020, Luigi Ferraro, Federico Galetto,
   Francesca Gandini, Hang Huang, Matthew Mastroeni, Xianglong Ni.

   You may redistribute this file under the terms of the GNU General Public
   License as published by the Free Software Foundation, either version 2 of
   the License, or any later version.
*-

newPackage(
        "InvariantRing",
        Version => "2.0", 
        Date => "October 16, 2020",
        Authors => {
	    {Name => "Luigi Ferraro", 
		 Email => "lferraro@ttu.edu", 
		 HomePage => "http://www.math.ttu.edu/~lferraro/"
		 },
             {Name => "Federico Galetto", 
		 Email => "f.galetto@csuohio.edu", 
		 HomePage => "https://math.galetto.org"
		 },
             {Name => "Francesca Gandini", 
		 Email => "fra.gandi.phd@gmail.com", 
		 HomePage => "https://sites.google.com/a/umich.edu/gandini/home"
		 },
	     {Name => "Hang Huang", 
		 Email => "hhuang235@tamu.edu", 
		 HomePage => "https://math.tamu.edu/~hhuang"
		 },
	     {Name => "Thomas Hawes", Email => "thomas.hawes@maths.ox.ac.uk"},
	     {Name => "Matthew Mastroeni", 
		 Email => "mmastro@okstate.edu", 
		 HomePage => "https://mnmastro.github.io/"
		 },
             {Name => "Xianglong Ni", 
		 Email => "xlni@berkeley.edu", 
		 HomePage => "https://math.berkeley.edu/~xlni/"}
             },
        Headline => "invariants of group actions",
	AuxiliaryFiles => true,
        DebuggingMode => true
        )



export {
    "GroupAction",    	      	  
    "finiteAction",    	       	  
    "FiniteGroupAction",    	  
    "group",	    	    	  
    "isAbelian",    	    	  
    "permutationMatrix",          
    "schreierGraph",	    	  
    "words",    	       	  
    "cyclicFactors",	    	  
    "DiagonalAction",	     	  
    "diagonalAction",	     	  
    "equivariantHilbert",    	  
    "equivariantHilbertSeries",   
    "weights",	      	      	  
    "actionMatrix",    	       	  
    "groupIdeal",    	     	  
    "hilbertIdeal",    	       	  
    "linearlyReductiveAction",	  
    "LinearlyReductiveAction",	  
    "action",	     	     	  
    "definingIdeal",              
    "DegreeBound",    	      	  
    "invariants",    	     	  
    "invariantRing",	    	  
    "isInvariant",    	      	  
    "reynoldsOperator",	       	  
    "UseLinearAlgebra",     	  
    "RingOfInvariants",	       	  
    "UseCoefficientRing",    	  
    "UseNormaliz",    	      	  
    "UsePolyhedra",    	      	  
    "hironakaDecomposition",   	  
    "molienSeries",    	       	  
    "primaryInvariants",    	  
    "secondaryInvariants",    	  
    "Dade",    	       	       	  
    "DegreeVector",    	       	  
    "PrintDegreePolynomial"    	  
    }


needsPackage("Elimination")
needsPackage("Normaliz")
needsPackage("Polyhedra")

GroupAction = new Type of HashTable


-------------------------------------------
--- GroupAction methods -------------------
-------------------------------------------

ring GroupAction := Ring => G -> G.ring

dim GroupAction := ZZ => G -> dim ring G


-------------------------------------------
--- LOAD AUXILIARY FILES ------------------
-------------------------------------------

load "./InvariantRing/FiniteGroups.m2"

load "./InvariantRing/AbelianGroups.m2"

load "./InvariantRing/LinearlyReductiveGroups.m2"

load "./InvariantRing/Invariants.m2"

load "./InvariantRing/Hawes.m2"

beginDocumentation()

load "./InvariantRing/InvariantRingDoc.m2"

load "./InvariantRing/FiniteGroupsDoc.m2"

load "./InvariantRing/AbelianGroupsDoc.m2"

load "./InvariantRing/LinearlyReductiveGroupsDoc.m2"

load "./InvariantRing/InvariantsDoc.m2"

load "./InvariantRing/HawesDoc.m2"

load "./InvariantRing/Tests.m2"


end


