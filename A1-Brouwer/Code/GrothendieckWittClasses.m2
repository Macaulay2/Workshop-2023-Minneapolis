-- We define GrothendieckWittClass to be a new type,
--    meant to represent the isomorphism class of a symmetric bilinear form
--    over a base field.
GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"

-- A class in GW can be constructed from a representing matrix
gwClass = method()
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
    if not isSquareAndSymmetric(beta.matrix) then(
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
    Kb := baseField(beta);
    Kg := baseField(gamma);
    
    -- Galois field case
    if instance(Kb, GaloisField) and instance(Kg, GaloisField) then (
	-- Returns an error if the underlying fields of the two classes beta and gamma are different
	if not Kb.order == Kg.order  then error "Error: these classes have different underlying fields";
	return gwClass(safeBlockSum(beta.matrix, substitute(gamma.matrix,Kb)))
	);
    
    -- remaining cases
    if not Kb === Kg then error "Error: these classes have different underlying fields";
    	return gwClass(safeBlockSum(beta.matrix, gamma.matrix))
    )

-- Method for taking a tensor product of two Grothendieck-Witt classes
gwMultiply = method()

gwMultiply(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    Kb := baseField(beta);
    Kg := baseField(gamma);
    
    -- Galois field case
    if instance(Kb, GaloisField) and instance(Kg, GaloisField) then (
	-- Returns an error if the underlying fields of the two classes beta and gamma are different
	if not Kb.order == Kg.order  then error "Error: these classes have different underlying fields";
	return gwClass(beta.matrix ** substitute(gamma.matrix,Kb))
	);
    
    -- remaining cases
    if not Kb === Kg then error "Error: these classes have different underlying fields";
    	return gwClass(beta.matrix ** gamma.matrix)
    )