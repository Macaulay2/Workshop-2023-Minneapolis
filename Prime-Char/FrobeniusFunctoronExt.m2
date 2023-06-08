needsPackage "TestIdeals"
FrobeniusonExt=method()
FrobeniusonExt(Ideal,ZZ) := (I,i) -> (
    -- This function takes an ideal in a polynomial ring and a local cohomology index and returns the frobenius map on the Matlis dual of local Cohomology of the quotient module.
R := ring(I);
d :=dim R;   
M2 := R^1/I;
j := d-i;
M1 := R^1/frobenius(1,I);
phi := map(M2,M1,1);
Ext^j (phi,R^1)
)
