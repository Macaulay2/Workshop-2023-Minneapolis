-- discForm computes the product of the entries of a list;
-- used to compute the discriminant of a diagonal form, where form is given by list of diagonal entries

discForm = method ()

discForm (List):=(ZZ) => (f) -> (
    -- Input: f = list of diagonal elements of quadratic form
    -- Output: Product of the entries
    disc:=1;
    --for loop counts the number of negative diagonal entries
    for i from 0 to (#f -1) do(
	 disc = disc* f_i;
	 );
    return disc;
      	     
    );
