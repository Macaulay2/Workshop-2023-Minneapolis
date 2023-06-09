-- Block sum "++" of a zero matrix with something else outputs the wrong thing
safeBlockSum = method()
safeBlockSum (Matrix, Matrix) := Matrix => (A,B) -> (
    if numColumns A == 0 then return B;
    if numColumns B == 0 then return A;    
    return A ++ B
    
)
