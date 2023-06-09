-- Note from Tommy: goal is to work over rationals for any prime
-- long term goal is to implement for any local field

-- This is already implemented in Sage! See source here:
-- https://github.com/sagemath/sage/blob/2f426a11f4c38e8e73ccd3f39374a74aa08e91ab/src/sage/arith/misc.py#L4743
-- Will adapt the python implementation in particular here


-- Being worked on by Andrew






-- Input:
-- a,b integers
-- p integer, prime or (TODO???) -1 which represents the archimedian place

-- Output:
-- integer 0, -1, or 1


hilbertSymbol = method()

hilbertSymbol(ZZ, ZZ, ZZ) := ZZ => (a,b,p) -> (
    if (a == 0 or b == 0) then (
        return 0
    );

    if p != 1 then (
        -- TODO, working on
    )
)