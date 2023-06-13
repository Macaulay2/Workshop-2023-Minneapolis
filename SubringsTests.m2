TEST ///
    assert(true == true)
///

TEST ///    
    R := QQ[x,y];
    P := QQ[p_1..p_3];
    S := subring {x^2, x*y, y^2};
    assert(isSubringElement(x^2, S))
    assert(not isSubringElement(x, S))
///
