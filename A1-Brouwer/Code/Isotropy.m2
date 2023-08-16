-- Boolean returning if a symmetric bilinear form is anisotropic	
isAnisotropic = method()
isAnisotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    k:=baseField(alpha);
    -- Ensure base field is supported
    if not (k === CC or instance(k,ComplexField) or k === RR or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    A:=alpha.matrix;
    n:= numRows(A);
    return n == anisotropicDimension(alpha)
    )

-- Boolean returning if a symmetric bilinear form is isotropic	
isIsotropic = method()
isIsotropic (GrothendieckWittClass) := (Boolean) => (alpha) -> (
    return (not isAnisotropic(alpha));
    )
