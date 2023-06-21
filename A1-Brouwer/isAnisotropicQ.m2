--Nikita, Tom, Zhaobo
load "getInvariants.m2"

--isAnisotropicQ takes n GWClass over QQ and returns a Boolean based on whether or not 
--the class is anisotropic
--unlike simplifyForm it can say for certain if the form is anisotropic or not since it
--does not use rationalPoints

isAnisotropicQ = method()
isAnisotropicQ (GrothendieckWittClass) := Boolean => (alpha) -> (
    A:= alpha.matrix;
    n:= numRows(A);
    kk:= ring A;

    if (not (kk===QQ)) then (error "GrothendieckWittClass is not over QQ");

    --check if form is degenerate
    signature := (getInvariants(alpha))_1;
    if (signature_1 > 0) then (return false);
    
    --if rank>=5, we can use signature do decide this
    --the non-degenerate form will be anisotropic iff all diagonal entries have same sign
    if (n>= 5) then (
        return ((signature_0 == n) or (signature_2 == n));
    );

    --if rank<=4, we need to take p-adic completions

);