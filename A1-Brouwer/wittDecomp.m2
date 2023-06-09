needsPackage "RationalPoints2"

--load "./GW-type.m2"
--load "./matrixBooleans.m2"
--load "./squarefreePart.m2"
--load "./easyUpperLeftTriangular.m2"
load "./diagonalize.m2"

--Nikita Borisov and Frenly Espino

-- Given a matrix A, wittDecomp decomposes A as a sum nH + Q, where n= number of hyperbolic forms, and Q is (presumed to be) anisotropic.  
-- Input:  A square symmetric matrix A over an exact field (not RR or CC)
-- Output: (n, Q), n=integer giving number of hyperbolic spaces, Q=anisotropic forms.

-- Future Work: Shoudl be able to input a bound.

wittDecomp = method()
wittDecomp (Matrix) := (ZZ,Matrix) => (A) -> (
    k:= ring A;   
  
    -- Add error in case the base field is RR or CC
    if (instance(k,InexactFieldFamily) or instance(k,RealField) or instance(k,ComplexField)) then error "Error: base field is inexact, use wittDecompInexact() instead";
    
    n:=numRows(A); --rank of matrix
    R:=k[x_0..x_(n-1)]; -- local variable ring
    
    -- Create quadratic from f from matrix A
    f:=sum (
        for i from 0 to (n-1) list (
            sum (for j from 0 to (n-1) list (A_(i,j)*x_i*x_j))
            )
        );
 
    -- will use rationalPoints package to see if f has a zero.  If so, A is isotropic.  
    -- rationalPoints package seeks a zero upto variable bound
    
    use k;
    solnPt:=new MutableList;
    solnFound = false;
    for bound from 1 to 10 do (
        solns:= new MutableList from rationalPoints(ideal(f),Bound=>bound);
        if (#solns >1) then(
            if solns#0 == toList(n:0) then (solnPt=solns#1) else (solnPt=solns#0);
            solnFound =true;
            break;
        );
    );


    --if no solutions found we assume the form is anisotropic
    if ( not solnFound) then (return (0,A));
    
    --if solution found for rank 2 form, then the form is purely hyperbolic
    if ((n==2) and  (det(A)==0))then (error "Matrix singular, run wittDecompGeneral instead" );
    if (n==2) then (return (1,matrix(k,{{}})));

    -- if found a solution, record it as row matrix z. Then find vector y such that <z,y> <>0.  
    -- z and y then generate a hyberbolic plane. 
    z:=matrix{toList solnPt}; --z as a row matrix
    zA:=z*A; --z*A
    y :=new MutableMatrix from matrix{toList(n:(0/1))};
    for i from 0 to (n-1) do (
        if (zA_(0,i) != 0) then (y_(0,i)=1; break;);
    );
    
    --now z and y span a copy of |H in the bilinear form
    --we need to find a basis of vectors orthogonal (wrt bilinear form) to x and y 
    orthoComp :=gens kernel((z||matrix(y))*A);   -- a (n-2) x n matrix.  
    
    --now recursively apply wittDecomp to orthoComp^T*A*orthoComp a (n-2)-by-(n-2) Gram matrix
    subComputation := wittDecomp(transpose(orthoComp)*A*orthoComp);
    
    -- subComputation_0 gives number of hyperbolic forms in (n-2) x (n-2) subform.  
    -- 1+ subComputation_0 is the number of hyperbolic forms in A
    -- subComputation_1 is the anisotropic part of A and (also) the subform part.
    
    return (1+subComputation_0, subComputation_1);
)


--wittDecomp method for InexactFieldFamily

-- wittDecompInexact calculates the Witt decomposition for matrix A over the fields kk=RR, CC
-- Function assumes that A has maximal rank. 

wittDecompInexact=method()

wittDecompInexact (Matrix) := (ZZ,Matrix) => (A) -> (
    k := ring A;
    
    -- checks that we have RR or CC as our field
    if not (instance(k,RealField) or instance(k,ComplexField) or instance(k,InexactFieldFamily)) then error "Error: base field is not RR or CC";
    
    n:=numRows(A); --rank of matrix
    
    --if k is the complex numbers, witt decomposition depends only on rank
    if (k===CC or instance(k,ComplexField)) then (
        if (n%2==0) then(return (n//2,matrix(CC,{{}}))) --if rank is even, then matrix decomposes into n/2 hyberbolic forms with no anisotropic parts
        else return (n//2,id_(k^1)); --if rank is odd, matrix decomposes into (n-1)/2 hyperbolic forms with 1-by-1 anisotropic part
        );
    
    --if k is the real numbers, witt decomposition depends on rank and signature
    if (k===RR or instance(k,RealField)) then (
        diagA := diagonalize(A);
        posEntries := 0; --for loop counts the number of positive diagonal entries of diagA
        negEntries := 0; --for loop counts the number of negative diagonal entries
	for i from 0 to (n-1) do(
            if diagA_(i,i)>0 then(
                posEntries=posEntries+1;
            );
	    if diagA_(i,i)<0 then(
                negEntries=negEntries+1;
            );
        );

        if (posEntries + negEntries > n) then (print "A is singular");
        wittIndex := min(posEntries,negEntries); -- witt index is given by how many positive-negative diagonal entry pairs exist
        signature := posEntries-negEntries; 
        if signature == 0 then (return (wittIndex,matrix(RR,{{}})))
        else if signature > 0 then ( return (wittIndex, id_(k^(signature)))) --signature characterizes anisotropic part
        else return (wittIndex, -id_(k^(-signature)));
        );
);


