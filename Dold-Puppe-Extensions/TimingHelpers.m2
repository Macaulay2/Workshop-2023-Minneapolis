race = (fun1,fun2,args,numIterations) -> (
    time M = (fun1 args);
    time N = (fun2 args);
    print(M==N);
    L = for i to numIterations list (
	    M = timing fun1 args;
	    N = timing fun2 args;
	    M#0>N#0
	);
    T = tally L;
    Tcat = tally (L|{true,false});
    -- Concattenating with true false so that 
    if T#false>T#true then print(toString(fun1) | " was slower than " | toString(fun2)) else print(toString(fun2) | " was slower than " | toString(fun1));
    T
    )
	    
