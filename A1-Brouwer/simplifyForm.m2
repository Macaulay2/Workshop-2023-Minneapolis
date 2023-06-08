-- Inputs a form in GW(k)
--     outputs a simple form for its Gram matrix,
--     and something of the form aH+b

load "./GW-type.m2"
load "./diagonalForm.m2"


-- Legendre boolean, returnts true if an element of a finite field is a square and false otherwise
legendreBoolean = method()
legendreBoolean (RingElement) := (Boolean) => a -> (
    if not instance(ring(a),GaloisField) then error "Error: this works only for Galois fields";
    q := (ring(a)).order;
    a^((q-1)//2) == 1 
);


-- simplifyForm method

simplifyForm = method()
simplifyForm (GrothendieckWittClass) := (GrothendieckWittClass, String) => beta -> (
    print("TEST");
    
    -- Check if the diagonalForm has already been computed, if so recall it from the cache    
    gamma := diagonalForm(beta);
    -- Get the matrix, base field, and rank of the form
    A := gamma.matrix;
    k := baseField(gamma);
    
    n := numRows(gamma.matrix);
    
    -- If the form has no rank, stop here and return the zero form. This gets rid of fringe cases later on
    if n == 0 then return (gamma,"zero form");
    
    -- Define the hyperbolic form over our field
    H := matrix(k,{{1,0},{0,-1}});
        
    -----------------------
    -- Field by field cases
    -----------------------
        
    -- In the case of C:
    if (k === CC or instance(k,ComplexField)) then(
        print("GOT HERE");
        if even n then(
        
            -- The form is (n/2)H
            numHypForms := n/2;
        
            simplifiedForm := H;
        
            for i in 1..sub((numHypForms-1),ZZ) do(
                simplifiedForm = simplifiedForm ++H;
		);
        
            outputString := toString(numHypForms) | "H";
            return (simplifiedForm,outputString)
            );
        if odd n then(
        
            -- The form is floor(n/2)H + <1>
            numHypForms := floor n/2;
            simplifiedForm := H;
        
            for i in 1..sub((numHypForms-1),ZZ) do(
                simplifiedForm = simplifiedForm ++H;
		);
        
            simplifiedForm = simplifiedForm ++ matrix(k,{{1}});
        
            outputString := toString(numHypForms) | "H + <1>";
            return (simplifiedForm,outputString)
            );
	);
    	-- End CC case
    
       
    -- If the field is R, look at the sign of elements along the diagonal
    if (k === RR or instance(k,RealField)) then(
    	print("field is R");
        posEntries := 0; --for loop counts the number of positive diagonal entries of diagA
        negEntries := 0; --for loop counts the number of negative diagonal entries
    	for i from 0 to (n-1) do(
            if A_(i,i)>0 then(
                posEntries=posEntries+1;
            	);
            if A_(i,i)<0 then(
                negEntries=negEntries+1;
            	);
            );
            
    	-- Number of hyperbolic forms is the number of positive and negative entries
        wittIndex := min(posEntries,negEntries);
    	
    	simplifiedForm := H;
            
    	for i in 1..(wittIndex-1) do(
	    simplifiedForm = simplifiedForm ++H;
            );
        
    	-- Look at what's left over
    	posEntries = posEntries - wittIndex;
    	negEntries = negEntries - wittIndex;
        
    	if posEntries > 0 then(
            for i in 1..posEntries do(
        	simplifiedForm = simplifiedForm ++ matrix(k,{{1}});
        	);
            
            outputString := toString(wittIndex) | "H + " | toString(posEntries) | "<1>";
            
	    return (simplifiedForm,outputString)
	    );	        
        
    	if negEntries > 0 then(
            for i in 1..negEntries do(
        	simplifiedForm = simplifiedForm ++ matrix(k,{{-1}});
        	);
        
            outputString := toString(wittIndex) | "H + " | toString(negEntries) | "<-1>";
	    return (simplifiedForm,outputString)
            );           
        );
    	-- End RR case
    
    
    --------
    -- Finite field case
    --------
    if instance(k, GaloisField) then(
	
	-- Take the diagonal entries to a list
	diagEntries := {};
	for i from 0 to sub((n-1),ZZ) do(
	    diagEntries = append(diagEntries, A_(i,i));
	    );
	
	-- Start counting squares and nonsquares
        numSquares = 0;
	numNonSquares = 0;
	
	-- Define a new variable which we will set equal to a nonsquare representative
	nonSquareRepresentative := sub(0,k);
	
	
	-- Check if each elt is a square or not
	for x in diagEntries do(
	    
	    -- If it is a square
	    if legendreBoolean(x) then(
		numSquares = numSquares + 1;
		);
	    
	    -- If it is not a square
	    if not legendreBoolean(x) then(
		--Pick x to be a representative of the non-square class in F_q^x / (F_q^x)^2
		nonSquareRepresentative = x;
		numNonSquares = numNonSquares + 1;
		);
	    );
	
	-- Number of hyperbolic forms
        wittIndex := floor(numSquares/2) + floor(numNonSquares/2);
	print("Witt index");
	print(wittIndex);
	
	
	-- If the wittIndex is positive:
	
	if wittIndex > 0 then(
	    
	    print("I can see witt index positive");
	    outputString := toString(wittIndex) | "H";
	    outputMatrix := H;
	    for i in 0..(wittIndex-1) do(
		outputMatrix = outputMatrix ++ H;
		
		);
	    
	    if odd numSquares then(
		outputMatrix = outputMatrix ++ matrix(k,{{1}});
		outputString = outputString | " + <1>";
		);
	    
	    if odd numNonSquares then(
		print("I see odd num squares");
		outputMatrix = outputMatrix ++ matrix(k,{{nonSquareRepresentative}});
		outputString = outputString | " + <" | toString(nonSquareRepresentative) | ">";
		);
	    
	    return(outputMatrix,outputString)
	    
	    );
	
	-- If the wittIndex is negative
	if wittIndex == 0 then(
	    outputString := "";
	    
	    if odd numSquares then(
		outputMatrix := matrix(k,{{1}});
		outputString = "<1>";
		
		if odd numNonSquares then(
		    outputMatrix = outputMatrix ++ matrix(k,{{nonSquareRepresentative}});
		    outputString = outputString | " + <" | toString(nonSquareRepresentative) | ">";
		);
		
	        return (outputMatrix,outputString)
		
		);
	    
	    if odd numNonSquares then(
		outputMatrix = matrix(k,{{nonSquareRepresentative}});
		outputString = " <" | toString(nonSquareRepresentative) | ">";
		
		
		return (outputMatrix,outputString)
		);	    
	    );    
        );
    	-- End finite field case
    
    
    );




